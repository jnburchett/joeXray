from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.cosmology import Planck18 as cosmo

import numpy as np
import heasoftpy as hsp
from joeXray import events,getData

import subprocess

def center_crop(infile,outfile,coords,boxdims=(30.*u.arcmin,30.*u.arcmin),
                binning=80,event_file=True,emin=0.2*u.keV,emax=10.0*u.keV,
                 clobber=True):
    if event_file:
        ff = fits.open(infile)
        w= WCS(naxis=2)
        w.wcs.crpix = [ff[1].header['TCRPX4'],ff[1].header['TCRPX5']]
        w.wcs.cdelt = [ff[1].header['TCDLT4'],ff[1].header['TCDLT5']]
        w.wcs.ctype = [ff[1].header['TCTYP4'],ff[1].header['TCTYP5']]
        w.wcs.crval = [ff[1].header['TCRVL4'],ff[1].header['TCRVL5']]
    else:
        w = WCS()

    coox,cooy = w.world_to_pixel(coords)
    xdelta = boxdims[0].to(u.deg).value / w.wcs.cdelt[0]
    ydelta = boxdims[1].to(u.deg).value / w.wcs.cdelt[1]
    x1,x2 = (coox-xdelta,coox+xdelta)
    print(xdelta,ydelta)
    y1,y2 = (cooy-ydelta,cooy+ydelta)
    spatialpart='[bin x = {}:{}:{}, y = {}:{}:{}]'.format(int(x1),int(x2),binning,int(y1),int(y2),binning)

    if event_file:
        eventpart='[EVENTS]'
    else: 
        eventpart=''
    
    energypart = events.filter_energy(emin,emax)

    if clobber: clob='yes'
    else: clob='no'

    filt = '{}{}{}'.format(eventpart,energypart,spatialpart)
    infile4copy = '{}{}'.format(infile,filt)
    print(infile4copy)

    hsp.ftcopy(infile=infile4copy,outfile=outfile,clobber=clob)

def image_extent(coords,z,extent,outname,**kwargs):
    scale = cosmo.kpc_comoving_per_arcmin(z)
    skysep = extent / scale
    print(skysep)
    evf = getData.grab_evt(coords)
    evf = evf.split('/')[-1]
    #import pdb; pdb.set_trace()
    imout = center_crop('events/'+evf,outname,coords,
                               boxdims=[skysep,skysep],**kwargs)
    
def smooth_image(inimage,outname,sigma,exact=False):
    '''if exact:
        exact = '-exact'
    else:
        exact = ''
    args_imsmo = ['imsmo',inimage]
    import pdb; pdb.set_trace()
    subprocess.run(args_imsmo)'''
    from scipy import ndimage
    ff = fits.open(inimage)

    
    


          
    