from .utils import *
from .data import *
from .products import *
import healpy as hp
from tqdm import tqdm
from lal import gpstime

def get_event_info(event):
    t_0_gps = event['t_0']
    t_0_utc = gpstime.gps_to_utc(t_0_gps)
    file_url = event['links']['files']
    file_dict = requests.get(file_url).json()
    flink = None
    for key, link in file_dict.items():
        if ('offline0.fits' in key or 'offline1.fits' in key) and 'Bilby' in key:
            flink = link
            print(flink)
            break
        elif not 'multiorder.fits' in key or 'Bilby' in key:
            if 'bayestar' in key:
                flink = link
                print(flink)
                break
            
    if flink is None:
        print('%s Cannot find skymap link' % event['superevent_id'])
                
    return (t_0_utc, t_0_gps, flink)

def EP_coverage(skymap_url,pointings):
    """
    Args:
        skymap_url (str): url of skymap file in gracedb api
        pointings (astropy.table): EP Pointings includes columns of ra,dec
    """
    wxt_radius = np.sqrt(75)/1.5
    skymap, meta = read_sky_map(skymap_url,nest=False,distances=True)
    EP_pnts = SkyCoord(np.array(pointings['ra']),np.array(pointings['dec']),unit=u.deg)
    prob_norm = skymap[0]/max(skymap[0])
    npix = len(skymap[0])
    nside = hp.npix2nside(npix)
    credible_levels = find_greedy_credible_levels(skymap[0])

    #Loop all sky
    N = 8000
    i = 0
    source90, sourceEP = 0, 0
    sample = np.random.choice(np.arange(0,len(prob_norm),1),N,replace=True,p=skymap[0])
    with tqdm(total=N) as pbar:
        for i,ipix in enumerate(sample):
            ra, dec = hp.pix2ang(nside,ipix,lonlat=True)
            CL = credible_levels[ipix]
            if CL < 0.9:
                source90 += 1
                pbar.update(1)
                c = SkyCoord(ra*u.deg,dec*u.deg)
                sep = c.separation(EP_pnts)
                if min(sep.value) < wxt_radius:
                    #In EP FoV
                    sourceEP += 1
                else:
                    continue
            else:
                continue
    return sourceEP/source90


def GW_EP_data_search(gw_info):
    pass


class Search(object):
    
    def __init__(self):
        self.one_orbit = 3600
        self.pointings = None
    
    def from_gw(self,gw_info,event_name,window=[-50,50],l23=None,l23_dir='/mnt/rdliang/AGN/BBH_AGN/data/L23/',save_pnt=False,pnt_dir='/mnt/rdliang/AGN/BBH_AGN/data/WXT_pnt',plot=False):
        #gw info
        self.t_0_utc, self.t_0_gps, self.flink = gw_info
        self.window = window
        #self.skymap, meta = read_sky_map(self.flink,distances=True)
        
        #output dir
        self.l23_dir = os.path.join(l23_dir,'%s_%s'%(event_name,l23))
        if not os.path.exists(self.l23_dir) and l23:
            os.mkdir(self.l23_dir)
        
        self.t_0_utc = Time(self.t_0_utc)
        self.tsearch_start = self.t_0_utc + self.window[0] * u.second - self.one_orbit * u.second
        self.tsearch_end = self.t_0_utc + self.window[1] * u.second + self.one_orbit * u.second
        self.tstart = self.t_0_utc + self.window[0] * u.second
        self.tend = self.t_0_utc + self.window[1] * u.second
        
        epdata = retrieve_epdata(start=self.tsearch_start.iso.replace(' ','T'),
                                 end=self.tsearch_end.iso.replace(' ','T'))
        
        maxoverlap, cov = 0., 0.
        
        if len(epdata) > 1 or isinstance(epdata,list):
            table_data = []
            for j in range(len(epdata)):
                obs_start = Time(epdata[j]['obs_start'])
                obs_end = Time(epdata[j]['obs_end'])
                obsid = epdata[j]['obs_id']
                cmosid = epdata[j]['detnam'][4:]
                maxoverlap = max(toverlap(self.tstart,self.tend,obs_start,obs_end),maxoverlap)
                
                if toverlap(self.tstart,self.tend,obs_start,obs_end) > 0:
                    
                    #Extract Level23 EVT
                    obsid = epdata[j]['obs_id']
                    cmosid = epdata[j]['detnam'][4:]
                    clevtfile = '/mnt/epdata_pipeline/L23/obs/%s/ep%swxt%spo_cl.evt'%(obsid,obsid,cmosid)
                    expfile = '/mnt/epdata_pipeline/L23/obs/%s/ep%swxt%s.exp'%(obsid,obsid,cmosid)
                    stemout = '%s_%sCMOS%s'%(event_name,obsid,cmosid)
                    stemout = stemout.replace(' ','')
                    table_data.append([obs_start.iso,epdata[j]['pnt_ra'],epdata[j]['pnt_dec'],clevtfile,expfile,stemout])
                    
                    #L23 peoducts
                    if l23:
                        #create sub dir
                        if not os.path.exists(os.path.join(self.l23_dir,'%sCMOS%s'%(obsid,cmosid))):
                            os.mkdir(os.path.join(self.l23_dir,'%sCMOS%s'%(obsid,cmosid)))
                            
                        if l23 == 'bb':
                            bb_products(obsid,cmosid,os.path.join(self.l23_dir,'%sCMOS%s'%(obsid,cmosid))+'/test.json')
                        elif l23 == 'default':
                            extract_products(clevtfile,expfile,os.path.join(self.l23_dir,'%sCMOS%s'%(obsid,cmosid)),
                                stemout=stemout,
                                ra=epdata[j]['pnt_ra'],dec=epdata[j]['pnt_dec'])
                        
            self.pnt_dir = os.path.join(pnt_dir,'%s_%s.fits'%(event_name,window[1]-window[0]))    
            self.pointings = Table(np.array(table_data),names=['t','ra','dec','evt','exp','stemout'],dtype=['U30']+['f8']*2+['U75']*3)
            
            if maxoverlap > 0.:
                cov = EP_coverage(self.flink,self.pointings)
                print('Time overlap %.1f s & Sky coverage %.2f for %s'%(maxoverlap,cov,event_name))
             
        if save_pnt:
            self.pointings.write(self.pnt_dir,format='fits',overwrite=True)
        return maxoverlap, cov, self.pointings
    
    
    
    
    def from_fermi(self):
        pass
    
    
    def from_tns(self):
        pass