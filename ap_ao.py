#!/usr/bin/env python3

"""`ap_ao.py`

Python aperture photometry script, designed to run on HST WFPC2 "c0m" fits files (https://mast.stsci.edu/search)

Will run over all HST WFPC2 c0m files in its directory.

Parameters:
-------    

The data will be dsiplayed with ds9 and user will be prompted to estimate 
the center of the extended object, and estimate the galaxy type
after executing the script (elliptical has the most overlap with HST bands).

Returns:
--------

-Figure displaying apertures and flux [.png]

-Figure displaying the profiles of the object with the Moffat fit and FWHM [.png]

-Object name, aperture summed flux, FWHM, converted Johnson R band magnitude [.txt]

TBD
-----
- make useful for PCI and WPC3 [ ]

"""

__version__ = '20221121'
__author__ = 'bcooper, astephens'

import warnings
warnings.filterwarnings("ignore")
import argparse
from astropy.io import fits as fits
from astropy.modeling import models, fitting
from astropy.nddata import CCDData
from astropy import units
import ccdproc as ccdp
import glob
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import os
import photutils
import photutils.aperture as pap
from photutils.aperture import CircularAperture,CircularAnnulus
from scipy.interpolate import interp1d

# Ensuring XPA methods match for pyds9 import, and pointing to stsynphot config files
os.environ['XPA_METHOD']='localhost' 
os.environ['PYSYN_CDBS']='/Users/bcooper/Documents/PT_AO/trds'
import pyds9
from synphot import SourceSpectrum, Observation
import stsynphot


def ds9_disp(img, deletef='no', newf='no', tile_yn='no'):
    """Looks for and connects to open DS9 Instance, starts one otherwise. 
    Displays object for user to estimate object's 
    central pixel coodinates. Inputs are image, desired frame and tile(yes or no). 
    Returns ds9 instance.
    """
    ds9tgts = pyds9.ds9_targets()
    if ds9tgts:
        d = pyds9.DS9()
    else:
        d = pyds9.DS9(start=True)
    d.set('iconify no')
    if deletef == 'yes':
        d.set('frame delete all')
    if newf == 'yes':
        d.set('frame new')
    d.set(f'tile {tile_yn}')
    d.set_np2arr(img)
    #d.set('scale mode minmax')
    d.set(f'scale limits {np.percentile(img,1)} {np.percentile(img,99.999)}')
    d.set('scale lock limits')
    d.set('scale log')
    return d


def ap_disp(flux,x,y,objName,name):
    """Display aperture and sky annulus overtop of the user-provided target. 
    Inputs are the estimated central coordinates x,y, target name, and filename. 
    displays a plot and saves the file as a png.
    """
    ap    = CircularAperture((x,y), r=11.0)
    skyap = CircularAnnulus((x,y), 11.0, 12.1) 
    plt.figure()
    crop_f = flux[int(y-20):int(y+20), int(x-20):int(x+20)]
    # Keep positive vals for log normalization on colour bar
    plt.imshow(flux, norm = LogNorm(vmin = np.abs(np.nanmin(crop_f)), vmax = np.abs(np.nanmax(crop_f)))) 
    plt.title(f"{objName} Flux (erg/cm2/s/A) and Apertures", fontsize=15)
    plt.colorbar()
    ann_patches = skyap.plot(color='white', lw=0.5, label = 'Annulus (Width = 0.05")')
    ap_patches = ap.plot(color='red', lw=0.5,label='Aperture (Radius = 0.5")')
    # Plot close to object to see difference in apertures
    plt.xlim(x-50,x+50)
    plt.ylim(y-50,y+50)
    handles = (ap_patches[0], ann_patches[0])
    plt.legend(facecolor='#458989', labelcolor='white', handles=handles, prop={'weight': 'bold', 'size': 11})
    plt.xlabel(f"{name} Pixels")
    plt.ylabel("Pixels")
    plt.savefig(f'{objName}_aperture.png', dpi=1200, format='png')
    plt.show(block=False)
        

def get_centre(ds9,flux,objName,name):
    """Prompt user to enter central coordinates, and asks for confirmation. 
    Input is current ds9 instance from ds9_disp, image numpy arry of flux, target and filename. 
    Depends on ap_disp(flux,x,y,objName,name).
    Returns coordinates as a list of integers.
    """
    while True:
        print('\nPlease click on approximate center of the object...')
        coords1=ds9.get('iexam coordinate image ')
        coords2=coords1.rsplit()
        coords=(int(float(coords2[0])),int(float(coords2[1])))
        ap_disp(flux,coords[0],coords[1],objName,name)
        happy=input('\n Does this aperture look centred? (y/n): ')
        plt.close()
        if happy == 'y':
            print('\n Yay!')
            break
    return coords  


def lin_interp(x, y, i, half):
    """Linear Interpolation to find the location of the 
    fwhm crossings on the Moffat fit.
    """
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))


def half_max_x(x, y):
    """Determine Half Max. Uses lin_interp, 
    and requires x and y numpy arrays.
    """
    half    = ((max(y)-min(y))/2.0)+min(y)
    signs   = np.sign(np.add(y, -half))
    roots   = (signs[0:-2] != signs[1:-1])
    roots_i = np.where(roots)[0]
    return [lin_interp(x, y, roots_i[0], half),lin_interp(x, y, roots_i[1], half)]


def magT(objName, mfilt, filt, mfilt_u=None): 
    """Converts magnitude to Johnson R magnitude using galaxy profiles and 
    the corresponding Hubble filter. Returns R Mag.
    """
    fdictionary={
            'F300W' : stsynphot.band('wfpc2, f300w'),
            'F439W' : stsynphot.band('wfpc2, f439w'),
            'F450W' : stsynphot.band('wfpc2, f450w'),
            'F547M' : stsynphot.band('wfpc2, f547m'),
            'F555W' : stsynphot.band('wfpc2, f555w'),
            'F569W' : stsynphot.band('wfpc2, f569w'),
            'F606W' : stsynphot.band('wfpc2, f606w'),
            'F622W' : stsynphot.band('wfpc2, f622w'),
            'F656N' : stsynphot.band('wfpc2, f656n'),
            'F658N' : stsynphot.band('wfpc2, f658n'),
            'F675W' : stsynphot.band('wfpc2, f675w'),
            'F702W' : stsynphot.band('wfpc2, f702w'),
            'F791W' : stsynphot.band('wfpc2, f791w'),
            'F814W' : stsynphot.band('wfpc2, f814w'),
    }
    johnson_r = stsynphot.band('johnson, r') 
    gdictionary={
        'bulge'      : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'bulge_template.fits')),
        'elliptical' : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'elliptical_template.fits')),
        'sc_galaxy'  : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'sc_template.fits')),
        'qso'        : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'qso_template.fits')),
        'seyfert1'   : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert1_template.fits')),
        'seyfert2'   : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert2_template.fits')),
        'b'          : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'bulge_template.fits')),
        'e'          : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'elliptical_template.fits')),
        'sc'         : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'kc96', 'sc_template.fits')),
        'q'          : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'qso_template.fits')),
        's1'         : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert1_template.fits')),
        's2'         : SourceSpectrum.from_file(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert2_template.fits')),
    }
    while True:
        gkine = input('\nPlease input galaxy type. \nbulge(b), elliptical(e), sc_galaxy(sc), qso(q), seyfert1(s1), seyfert2(s2): ')
        if gkine in gdictionary.keys():
            break
        else:
            gkine = input('\nOops! Please input one of these types: \nbulge(b), elliptical(e), sc_galaxy(sc), qso(q), seyfert1(s1), seyfert2(s2): ')
            break
    mag_in  = mfilt * units.STmag    
    spec    = gdictionary[gkine].normalize(mag_in, band=fdictionary[filt])
    obs     = Observation(spec, band=johnson_r)
    mag_out = obs.effstim(flux_unit='vegamag', vegaspec=stsynphot.Vega)
    print(f'\n{objName}: {filt} {mag_in:.3f} {gkine} galaxy -> Johnson R {mag_out:.3f}')
    
    if mfilt_u:
        mag_in    = mfilt_u * units.STmag    
        spec      = gdictionary[gkine].normalize(mag_in, band=fdictionary[filt])
        obs       = Observation(spec, band=johnson_r)
        mag_out_u = obs.effstim(flux_unit='vegamag', vegaspec=stsynphot.Vega)
        
    print(f'\n{objName} (No Background Subtraction): {filt} {mag_in:.3f} {gkine} galaxy -> Johnson R {mag_out_u:.3f}')
    return mag_out, mag_out_u



def main(args): 
    
    ## Open and run on all fits files in current dir
    fn_list = glob.glob( '*.fits')
    
    ## Create file for results output
    if args.filename:
        fout = args.filename 
    else:
        fout ='results.txt'
    f      = open(f'{fout}', 'w')
    header = ['OBJECT ', 'FWHM(") ', 'AP_FLUX(erg/cm2/s/A) ', 'AP_FLUX_UNSUB(erg/cm2/s/A) ', 'SKY_MEAN(erg/cm2/s/A) ', 'HST_FILT ', 'FILT_MAG ', 'R_MAG ', 'R_MAG_UNCERT ', 'R_MAG_UNSUB ']
    f.writelines(header)
    f.close()
    
    for num, fn in enumerate(fn_list):
        
        ## Get useful keywords from header 
        print(f'\nCurrent Image: {fn}\n')
        name = fn
        file = fits.open(name)
        file.verify('fix')
        imgs     = file[1].data  #fits data
        objName  = file[0].header['TARGNAME']
        photflam = file[1].header['photflam']
        exptime  = file[0].header['exptime']
        photzpt  = file[1].header['photzpt']
        filt     = file[0].header['FILTNAM1']
        
        ## Run cosmic ray rejection algorithm TBD: display result and prompt user if target was affected
        data       = CCDData.read(name, hdu=1, unit='adu', format='fits')
        clean_data = ccdp.cosmicray_lacosmic(data, readnoise=10, sigclip=1, verbose=False)
        imgs       = np.asarray(clean_data)
        ds9        = ds9_disp(imgs, deletef='yes', newf='yes')
        imgs       = np.asarray(data)
        ds9        = ds9_disp(imgs, newf='yes', tile_yn='yes')
        cr_yn      = input('\nDoes the target centre (left img) look corrupted by the cosmic ray rejection algorithm ? (y/n): ')
        
        if cr_yn == 'y':
            ds9 = ds9_disp(imgs,newf='yes', tile_yn='no')
            print('\n Got it, using the raw image instead.')
        else:
            imgs = np.asarray(clean_data)
            ds9  = ds9_disp(imgs,newf='yes', tile_yn='no')
            print('\n Okay, sticking with the cosmic ray rejection image.')
        
        ## Convert data to flux in units of erg/cm2/s/A
        flux = imgs * photflam / exptime
        
        ## Get estimate for central coordinates of obj
        coords    = get_centre(ds9,flux,objName,name)
        dsf       = np.ones((1,2))
        dsf[0][0] = coords[0]
        dsf[0][1] = coords[1]#getting the proper format of coordinates for DAOstarfinder

        print("\nPlease Hold...")
        
        ## Centroid Detection using initial guess as starting point
        try:
            galaxy       = photutils.detection.DAOStarFinder(fwhm=12, threshold=5, brightest=1, xycoords=dsf).find_stars(data=imgs)
            x            = float(galaxy['xcentroid'].data)
            y            = float(galaxy['ycentroid'].data)
            final_coords = (x,y)
            print(f'\n Centroid Coords: {x:.2f} {y:.2f}')
            ap_disp(flux,final_coords[0],final_coords[1],objName,name)
        except:
            print(f'\nDAO StarFinder failed. \nContinuing with initial guess {coords}, as central coordinates.')
            x            = coords[0]
            y            = coords[1]
            final_coords = coords 
            
        ## Aperture on centroid + sky ap (radius = 11 px or 0.5", annulus width of 0.05" or 1.1 px)
        ap             = CircularAperture(final_coords, r=11.0)
        skyap          = CircularAnnulus(final_coords, 11.0, 12.1) # annulus
        skymask        = skyap.to_mask(method='center')
        sky            = skymask.get_values(imgs/exptime) #[ADU/s]
        bkg            = np.mean(sky) #[ADU/s]
        for i in range(np.shape(imgs)[0]):
            for j in range(np.shape(imgs)[1]):
                if imgs[i][j]<0:
                    imgs[i][j]=0
        error          = np.sqrt(imgs/exptime) #[ADU/s]
        phot_table_sub = pap.aperture_photometry(imgs/exptime-bkg, ap, error=error) # with background subtraction [ADU/s]
        phot_table     = pap.aperture_photometry(imgs/exptime, ap, error=error) # without background subtraction [ADU/s]
        bgsf           = flux - bkg * photflam #background subtracted flux [erg/cm2/s/A]
        uncert         = 2.5 * np.log10(1 + np.sqrt((phot_table_sub['aperture_sum_err'].data)**2 + (np.std(sky)/(np.sqrt(np.shape(imgs)[0] * np.shape(imgs)[1])))**2))
        print(f"\nSummed Aperture Flux: {float(phot_table_sub['aperture_sum'].data * photflam)} erg/cm2/s/A")
        
        
        ## Radial Profile and FWHM - Moffat 2D Fit to Background Subtracted Flux
        rad = 12 
        k   = bgsf[int(y-rad):int(y+rad),int(x)]
        j   = bgsf[int(y),int(x-rad):int(x+rad)]
        y1  = (j+k)/2
        
        for val in range(0,len(y1)):
            if y1[val ] < 0:
                y1[val] = 0
        x1     = np.linspace(1,rad*2,rad*2)
        xplot  = np.linspace(-rad,rad,rad*2)
        g_init = models.Moffat2D(amplitude=1, x_0=rad,y_0=rad, gamma=1, alpha=1) #Moffat 2D   
        fitg   = fitting.LevMarLSQFitter()
        g      = fitg(g_init,x1,j,k,)
        gi     = interp1d(x1,g(x1,y1),kind='cubic') #Interpolated Moffat
        xi     = np.linspace(1,rad*2,rad*8)
        xiplot = np.linspace(-rad+0.5,rad-0.5,rad*8)
        halfg  = ((max(g(x1,y1))-min(g(x1,y1)))/2.0)+min(g(x1,y1))
        lab    = 'Moffat 2D Fit'
        
        # If the 2D Moffat failed (crossing points cannot be located,then perform Moffat 1D instead)
        # Find the two crossing points
        try:
            hmx = half_max_x(xi,gi(xi)) 
        except:
            print('Moffat 2D failed, attempting Moffat 1D instead.')
            g_init = models.Moffat1D(amplitude=1, x_0=rad, gamma=1, alpha=1) #Moffat 1D
            fit_g  = fitting.LevMarLSQFitter()
            g      = fit_g(g_init, x1, y1)
            gi     = interp1d(x1,g(x1),kind='cubic') #Interpolated Moffat
            hmx    = half_max_x(xi,gi(xi))
            halfg  = ((max(g(x1))-min(g(x1)))/2.0)+min(g(x1))
            lab    = 'Moffat 1D Fit'
        fwhm = hmx[1] - hmx[0]
        print(f'\nMoffat FWHM:{fwhm:.1f} pixels')
        
        ## Profile, Moffat, FWHM Plot and Save Fig
        plt.figure()
        plt.scatter(xplot,k, c='b', label='Vertical Profile')
        plt.scatter(xplot,j, c='r',  label='Horizontal Profile')
        plt.plot(xplot,y1, label='Mean Profile')
        plt.title(f'{objName} Radial Profile', fontsize=15)
        plt.xlabel(f'{name} Pixels')
        plt.ylabel('Background Subtracted Flux (erg/cm2/s/A)')
        plt.plot(xiplot,gi(xi), label=lab)
        plt.plot([hmx[0]-rad-0.5,hmx[1]-rad-0.5], [halfg, halfg],'g--', label=f'FWHM:{fwhm/11:.2f}(")')
        plt.legend(facecolor='#458989', labelcolor='white', prop={'weight': 'bold', 'size': 11})
        plt.ylim(0,max(np.nanmax(j), np.nanmax(k), np.nanmax(y1), np.nanmax(gi(xi))))
        plt.xlim(-rad+0.5,rad-0.5)
        plt.savefig(f'{objName}_profile.png', dpi=1200, format='png')
        plt.show(block=False)
        
        # Magnitude Computation and conversion
        mfilt = -2.5 * np.log10(float(phot_table_sub['aperture_sum'].data * photflam)) + photzpt
        mfilt_usub = -2.5 * np.log10(float(phot_table['aperture_sum'].data * photflam)) + photzpt
        print(f'\n{filt} Magnitude: {mfilt:.3f}')
        rmag, rmag_usub = magT(objName=objName, mfilt=mfilt, filt=filt, mfilt_u=mfilt_usub)
        
        
        # Append results to file
        f = open(f'{fout}', 'a')
        vals=[f'{objName} ', f'{fwhm/11:.3f} ', f"{float(phot_table_sub['aperture_sum'].data * photflam)} ", f"{float(phot_table['aperture_sum'].data * photflam)} " , f'{bkg * photflam} ' , f'{filt} ' ,f'{mfilt:.3f} ', f'{rmag:.3f} ', f'{uncert[0]:.3f} ', f'{rmag_usub:.3f} ']
        f.write('\n')
        f.writelines(vals)
        f.close()
        
        ds9.set('exit')
        plt.close('all')
            
    return
        

        
def create_parser():
    """Argument parser, will display help and take in arguments."""
    parser = argparse.ArgumentParser(description="Gemini North Altair target aperture flux assessment. Runs over HST WFPC2 c0m files in current directory, and the user will be prompted for inputs.", epilog='Version: ' + __version__)
    parser.add_argument('-f', '--filename', nargs='?', default=None, help='Optional: customize filename for output results' )
    return parser



if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
