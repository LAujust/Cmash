import requests
import numpy as np
import ligo.skymap
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.postprocess import find_greedy_credible_levels
from astropy.table import Table, join, join_skycoord, vstack
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
from ligo.skymap.postprocess import crossmatch
import astropy.units as u
from astropy.time import Time
import os, sys



def crossmatch_GW_AGN(skymap_dir,save_dir,wise_dir,milliquas_dir):
    """
    skymap_dir [str]:   skymap file dir or url
    save_dir[str]:      saved dir of 
    """

    #Load Cat
    #wise_agn_table = Table.read("/Users/liangrunduo/Desktop/Aujust/NAOC/EP/Crossmatch/catalogs/WISE_AGN.csv", format="csv")
    #milliquas_table = Table.read("/Users/liangrunduo/Desktop/Aujust/NAOC/EP/Crossmatch/catalogs/Milliquas.csv", format="csv")
    wise_agn_table = Table.read(wise_dir, format="csv")
    milliquas_table = Table.read(milliquas_dir, format="csv")

    skymap = read_sky_map(skymap_dir,moc=True)
    #skymap = read_sky_map('/Users/liangrunduo/EP/GW/S241102br_skymap.fits',moc=True)
    milliquas_table_valid = milliquas_table[milliquas_table['Z']>0]
    dist = Planck18.luminosity_distance(milliquas_table_valid['Z'])
    coordinates = SkyCoord(milliquas_table_valid['RA']*u.deg, milliquas_table_valid['DEC']*u.deg, dist)
    result = crossmatch(skymap, coordinates)
    #print(milliquas_table_valid[result.searched_prob_vol < 0.9])
    matched_milliquas = milliquas_table_valid[result.searched_prob_vol < 0.95]

    wise_agn_table_valid = wise_agn_table[wise_agn_table['z']>0]
    dist = Planck18.luminosity_distance(wise_agn_table_valid['z'])
    coordinates = SkyCoord(wise_agn_table_valid['_RAJ2000']*u.deg, wise_agn_table_valid['_DEJ2000']*u.deg, dist)
    result = crossmatch(skymap, coordinates)
    #print(wise_agn_table_valid[result.searched_prob_vol < 0.9])
    matched_wise = wise_agn_table_valid[result.searched_prob_vol < 0.95]
    matched_wise.rename_column('HMQ','NAME')
    for i in range(len(matched_wise)):
        if type(matched_wise['NAME'][i]) is not np.str_:
            matched_wise['NAME'][i] = matched_wise['WISEA'][i]


    # matched_milliquas.write('/Users/liangrunduo/EP/GW/crossmatch/S240413p_milliquas.csv',format='csv',overwrite=True)
    # matched_wise.write('/Users/liangrunduo/EP/GW/crossmatch/S240413p_wise.csv',format='csv',overwrite=True)

    if len(matched_milliquas) * len(matched_wise) > 0:
        matched_all = join(matched_milliquas, matched_wise, keys='NAME',join_type='outer')
    else:
        matched_all = vstack([matched_milliquas,matched_wise])
    print(matched_all)
    
    for i in range(len(matched_all)):
        if not isinstance(matched_all['RA'][i], np.floating):
            matched_all['RA'][i] = matched_all['_RAJ2000'][i]
            matched_all['DEC'][i] = matched_all['_DEJ2000'][i]

        if not isinstance(matched_all['Z'][i], np.floating):
            matched_all['Z'][i] = matched_all['z'][i]

    matched_all.write(save_dir,format='csv',overwrite=True)



def match_cat(source_cat,cat,radius,nthneighbor=1,seperation=False):
    idx, sep, _ = source_cat.match_to_catalog_sky(cat,nthneighbor=nthneighbor)
    filtered_id = sep < radius
    cat_matched_idx, cat_matched_sep = idx[filtered_id], sep[filtered_id]
    source_matched_idx = filtered_id
    if seperation:
        return source_matched_idx, cat_matched_idx, cat_matched_sep
    else:
        return source_matched_idx, cat_matched_idx
    
    
def toverlap(tstart, tend, obs_start, obs_end):
    """
    Args:
        tstart (astropy.time): window start
        tend (astropy.time): window end
        obs_start (astropy.time): obs start
        obs_end (astropy.time): obs end

    Returns:
        float: overlap time in sec
    """
    tt = min(Time(tend), Time(obs_end)) - max(Time(tstart), Time(obs_start))
    return max(0, tt.to_value('sec'))

def on_rm_error(func, path, exc_info):
    # Try to change the file permission and retry
    os.chmod(path, 0o777)
    func(path)
    
def clean_folder(folder_path):
    for root, dirs, files in os.walk(folder_path, topdown=False):
        # Remove files
        for file in files:
            file_path = os.path.join(root, file)
            try:
                os.remove(file_path)
                print(f"Removed file: {file_path}")
            except Exception as e:
                print(f"Failed to remove {file_path}: {e}")
        
        # Remove empty directories
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            try:
                os.rmdir(dir_path)
                print(f"Removed dir: {dir_path}")
            except Exception as e:
                print(f"Failed to remove dir {dir_path}: {e}")
    os.rmdir(folder_path)