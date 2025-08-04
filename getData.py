from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
import glob, shutil
import os, os.path, subprocess
from urllib import request

erodatadir = '/Users/jnburchett/Research/data/eROSITA/eRASS/ero_archive/'
skytilefile = erodatadir+'SKYMAPS_052022_MPE.fits'
skytiledat = Table.read(skytilefile)

def get_skytile(coords):
    rasel=(coords.ra.deg<skytiledat['RA_MAX'])&(coords.ra.deg>skytiledat['RA_MIN'])
    decsel=(coords.dec.deg<skytiledat['DE_MAX'])&(coords.dec.deg>skytiledat['DE_MIN'])
    matches = np.where(rasel&decsel)[0]
    if len(matches) == 1:
        tileid = skytiledat['SRVMAP'][rasel&decsel]
    elif len(matches) > 1:
        raise ValueError('Coordinates found in multiple tiles.')
    elif len(matches) <=0:
        raise ValueError('Coordinates not found in any tile:')
    tileid=tileid[0]
    return tileid

def grab_evt(coords,outdir='./events'):
    tileid = str(get_skytile(coords))
    if (coords.ra.deg < 10.)&(len(tileid)==4):
            tileid_rapart='00'+str(tileid[0])
    elif (coords.ra.deg < 15.)&(len(tileid)==4):
            tileid_rapart='00'+str(tileid[0])
    elif coords.ra.deg <100:
            tileid_rapart='0'+str(tileid)[:2]
    else:
            tileid_rapart = str(tileid[:3])
    '''
    effdec = 90.-coords.dec.deg
    if (coords.dec.deg < 10.)&(coords.dec.deg > 0.):
        tileid_decpart='00'+str(tileid[-1])
    elif (coords.dec.deg<0.)&((np.abs(coords.dec.deg)+90)<100):   
        #if coords.dec.deg 
        tileid_decpart='0'+str(tileid[-1])  
    else: 
        tileid_decpart=str(tileid[-3:])
'''
    tileid_decpart=str(tileid[-3:])
    tid=tileid_rapart+tileid_decpart
    evtpath = erodatadir+tid[-3:]+'/'+str(tid[:3])+'/EXP_010/'

    try:
        evt_fn = glob.glob(evtpath+'/*')[0]
    except:
        if not os.path.exists(evtpath):
             print('File not found, attempting download from eROdat archive.')
             if not os.path.exists(evtpath[:-8]):
                os.mkdir(evtpath[:-8])
             if not os.path.exists(evtpath):
                os.mkdir(evtpath)
             filepart = 'em01_{}{}_020_EventList_c010.fits.gz'.format(tileid_rapart,tileid_decpart)
             url2get = 'https://erosita.mpe.mpg.de/dr1/erodat/data/download/' + \
             '{}/{}/EXP_010/{}'.format(tileid_decpart,tileid_rapart,filepart)
             outfile = evtpath+filepart
             request.urlretrieve(url2get,outfile)
        evt_fn = glob.glob(evtpath+'/*')[0] 
        #import pdb; pdb.set_trace()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    try:
        shutil.copy2(evt_fn,outdir)
    except:
        import pdb; pdb.set_trace()
    return evt_fn