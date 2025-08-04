from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import heasoftpy as hsp
import os, subprocess


def extract_image(infile,outfile,emin=0.2*u.keV,emax=10.0*u.keV,binning=80,
                  events=True,clobber=False):
    energypart = filter_energy(emin,emax)
    if events:
        eventpart = '[EVENTS]'
    if binning is not None:
        spatialpart='[bin (x,y)={}]'.format(binning)
    else:
        spatialpart=''

    filt = '{}{}{}'.format(eventpart,energypart,spatialpart)
    infile4copy = '{}{}'.format(infile,filt)
    print(infile4copy)
    if clobber: clob='yes'
    else: clob='no'
    hsp.ftcopy(infile=infile4copy,outfile=outfile,clobber=clob)

def filter_energy(emin=0.2*u.keV,emax=10.0*u.keV):
    if emin is None:
        eminpart = ''
    else:
        eminpart = 'PI>{}'.format(int(emin.to(u.eV).value))
    if emax is None:
        emaxpart = ''
    else:
        emaxpart = 'PI<{}'.format(int(emax.to(u.eV).value))
    if (len(eminpart)>0)&(len(emaxpart)>0):
        energypart = '[({}&&{})]'.format(eminpart,emaxpart)
    else:
        energypart = '[{}{}]'.format(eminpart,emaxpart)
    
    return energypart

def extract_spectrum(evtfile,name,coords,output_dir='srctoolProducts',
                     apertype='circle',aperrad=120*u.arcsec,backtype='annulus',
                     backrad=[150,180]*u.arcsec):
    """Produce spectrum from events file

    Parameters
    ----------
    evtfile : str
        Name of event file 
    name : str
        Object name (used as a prefix for output files)
    coords : SkyCoord
        Center coordinates of extraction region
    output_dir : str, optional
        Name of directory for output files by default 'srctoolProducts'
    apertype : str, optional
        Type of extraction region, by default 'circle'
    aperrad : Quantity, optional
        Radius of circle or radii of annulus if tuple, by default 120*u.arcsec
    backtype : str, optional
        Type extraction region for background, by default 'annulus'
    backrad : list, optional
        Radius of circle or radii of annulus if tuple, by default [150,180]
    """    
    if output_dir[-1] != '/':
        output_dir=output_dir+'/'
    if apertype == 'circle':
        aperrad = int(aperrad.to(u.arcsec).value)
        srcreg = 'srcreg=fk5;circle * * {:n}"'.format(aperrad)
    if backtype == 'annulus':    
        backrad1 = int(backrad[0].to(u.arcsec).value)
        backrad2 = int(backrad[1].to(u.arcsec).value)
        backreg =  'backreg=fk5;annulus * * {:n}" {:n}"'.format(backrad1, backrad2)
    cmd = ['srctool',
      'eventfiles={}'.format(evtfile),
          'prefix={}_'.format(os.path.join(output_dir,name)),
          'suffix='.format(''),
          'srccoord=fk5; {:.5f} {:.5f}'.format(coords.ra.deg,coords.dec.deg),
           srcreg,
           backreg,
           'clobber=yes',
           'writeinsts=8'
          ]
    print(cmd)
    if os.path.isdir(output_dir)==False:
        os.mkdir(output_dir)
    subprocess.check_call(cmd)
    
