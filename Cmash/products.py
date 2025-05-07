from .utils import *
import subprocess

def extract_products(evtname,expname,out_dir,stemout,ra,dec,radius=67,regionfile=None,gtifile=None,binsize=1.,img=None):
    cmd = ['wxtproducts',
            'infile=%s'%evtname,
            'outdir=%s'%out_dir,
            'stemout=%s'%stemout,
            'bkgextract=yes',
            'imagefile=NONE',
            'expofile=%s'%expname,
            'binsize=%.1f'%binsize,
            'chatter=2',]
    if regionfile:
        cmd.append('regionfile=%s'%regionfile)
    else:
        cmd.append('ra=%.4f'%ra)
        cmd.append('dec=%.4f'%dec)
        cmd.append('radius=%s'%radius)  #3 arcmin
        
    if gtifile:
        cmd.append('gtifile=%s'%gtifile)

    if not img:
        cmd.append('imagefile=NONE')
    cmd.append('clobber=yes')
    print(cmd)
    subprocess.call(cmd)
    
    
def bb_products(obsid,cmosid,outdir): #Short-time transient search
    
    cmd = ['python',
           '/mnt/rdliang/wxt_analysis-master/data/epcodes/wxt_chain.py',
           obsid,
           cmosid,
           '/mnt/epdata_pipeline/L1/obs/%s'%obsid,
           '/mnt/epdata_pipeline/L23/obs/%s'%obsid,
           outdir] #outdir suggested to endwith **filtered1/test.json
    print(cmd)
    subprocess.call(cmd)