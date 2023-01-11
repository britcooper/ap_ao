# ap_ao

Ap_AO estimates the feasbility of potential extended object observations with Gemini
 North's Laser Guide Star (LGS) and adaptive optics unit [(Altair)](https://www.gemini.edu/instrumentation/altair) . 


## Usage

    usage: ap_ao.py [-h] [-f [FILENAME]]
    
    Gemini North Altair target aperture flux assessment. To be run over HST WFPC2 c0m
    files in current directory, and the user will be prompted for inputs.
    
    optional arguments:
      -h, --help            show this help message and exit
      -f [FILENAME], --filename [FILENAME]
                        Optional: customize filename for output results

To assess feasibilty, Ap_AO performs a background subtraction, sums the flux, and estimates the full width half maximum (FWHM) within Altair's 
LGS tip/tilt star aperture using an image of the target captured above the atmosphere.

Currently the image should only take the form of an HST WFPC2 "c0m" fits file, available from the [HST archive](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html) . 

Ap_AO will run through the ``.fits`` files in its working directory, one by one:

1. Runs a cosmic ray detection and removal algorithm using the LA-Cosmic method (user is prompted to confirm it did not corrupt the target of interest)
1. Displays the image with DS9 and prompts the user to click on the center of the object
1. Runs a DAO Star Finder centroiding algorithm with the inital guess from the user (if it fails, prompts user again and runs with user estimate)
1. Converts image pixel values to flux (erg/cm2/s/A)
1. Subtracts the background, sums the aperture flux, determines the uncertainty
1. Approximates the aperture FWHM using a 2D (or 1D if 2D fails) Moffat Fit
1. Converts the aperture flux to HST fiter magnitude, then to an aproximate Johnson R magnitude
1. Outputs a ``.png`` of the target with the aperture and sky annulus overlaid, and a ``.png`` of the target and interpolated Moffat profile
1. Outputs a text file (headers = OBJECT FWHM(") AP_FLUX(erg/cm2/s/A) AP_FLUX_UNSUB(erg/cm2/s/A) SKY_MEAN(erg/cm2/s/A) HST_FILT FILT_MAG R_MAG R_MAG_UNCERT R_MAG_UNSUB) 

## Dependencies

* python 3.7:
    * astropy
    * ccdproc
    * matplotlib
    * numpy
    * photutils
    * pyds9
    * scipy
    * synphot
    * stsynphot

* ds9

### Notes: 
1. Gemini DRAGONS python 3.7 conda environment with a couple additional `conda install` of libraries is an easy way to install and configure the python and ds9 requirements.
Install instructions [here](https://www.gemini.edu/observing/phase-iii/understanding-and-processing-data/data-processing-software/download-latest#dragons) . 
2. Stsynphot requires calibration files downloaded locally and an environment variable that must be set (line 54). Info [here](https://stsynphot.readthedocs.io/en/latest/index.html) .  
3. Pyds9 and ds9 connection issues have been a [problem in the past](https://github.com/astroconda/astroconda/issues/86) when using a conda environment. 
To thwart this,`XPA_METHOD` value must be set to `localhost` in both the DS9 environemnt and the python environment (line 53).

## Future Improvements (TBD)

Some useful updates to be added in the future (TBD):

* Make useful for HST PCI and WPC3 to expand range of targets 


## Support

Contact brittney.cooper@noirlab.edu with questions, problems, suggestions.


