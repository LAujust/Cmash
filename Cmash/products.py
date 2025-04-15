from .utils import *
import subprocess

def extract_products(evtname,out_dir,stemout,ra,dec,expname,regionfile=None,binsize=1.):
    cmd = ['wxtproducts',
            'infile=%s'%evtname,
            'outdir=%s'%out_dir,
            'stemout=%s'%stemout,
            'bkgextract=yes'
            'imagefile=NONE',
            'expofile=%s'%expname,
            'binszie=%.1f'%binsize,
            'chatter=1',]
    if regionfile:
        cmd.append('regionfile=%s'%regionfile)
    else:
        cmd.append('ra=%.4f'%ra)
        cmd.append('dec=%.4f'%dec)
        cmd.append('radius=%s'%67)  #3 arcmin
    cmd.append('clobber=yes')
    print(cmd)
    subprocess.call(cmd)