# ap_ao

Ap_AO estimates the feasbility of potential extended object observations with Gemini North's Laser Guide Star (LGS) and adaptive optics unit (Altair). 


## Usage

To assess feasibilty, Ap_AO performs background subtraction and sums the flux within Altair's 
LGS tip/tilt aperture using an image of the potential target observed above the atmosphere.

Currently, the image can only take the form of a HST WFPC2 "c0m" fits file, available on the [HST archive](https://mast.stsci.edu/search) .

Ap_AO will run through all ``.fits`` files in its working directory:

1. Runs a cosmic ray detection and removal algorithm using the LA-Cosmic method
1. Displays the image with DS9 and prompts the user to estimate the center of the object
1. Runs a DAO Star Finder centroiding algorithm with the inital guess from the user
1. Converts the image array values to flux (erg/cm2/s/A)
1. Subtracts the background and sums the aperture flux
1. Determines the aperture full width half maximum (FWHM) using a 2D (or 1D if 2D fails) Moffat Fit
1. Converts the aperture flux to ST magnitude, then to aproximate Johnson R magnitude
1. Outputs a ``.png`` of the target displaying the aperture and background annulus, and a ``.png`` of the Moffat profile
1. Outputs a text file with the object name, aperture FWHM, aperture flux (erg/cm2/s/A), HST filter magnitude and Johnson R magnitude

An optional bash wrapper script `ap_ao.sh` has been created to activate an appropriate environment in which ap_ao.py can be run, 
and sets the `XPA_METHOD` key to dissuade pyds9 connection issues. More into in **Dependencies** and **Potential Issues**.


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

* Optional: A Gemini DRAGONS `geminiconda3` conda environment in py3.7 with a couple additional `conda install`s will cover the python and ds9 requirements.
Install instructions [here](https://www.gemini.edu/observing/phase-iii/understanding-and-processing-data/data-processing-software/download-latest#dragons) . 
The optional wrapper script `ap_ao.sh` activates this conda environment before calling `ap_ao.py`.

## Potential Bugs

Pyds9 and ds9 connection issues have been a [problem in the past](https://github.com/astroconda/astroconda/issues/86) when using a conda environment. 

To thwart this, the `XPA_METHOD` key value can be set to `localhost` 
in both the DS9 environemnt and the python environment. 

It is necessary for `ap_ao.py` to set this variable before `import pyds9` for it to be effective.

The optional bash wrapper also sets this key value in the ds9 environment before calling `ap_ao.sh`. 
If `ap_ao.sh` is not used or modified for your local env, and you are using a conda install of ds9, this line will likely need to be run after your environment is activated.


## Future Improvements (TBD)

Some useful updates to be added in the future (TBD):

* Improve for CR detection (i.e. show and ask user if they want it or if it inadvertantly the target of interest) 
* Make useful for HST PCI and WPC3 to expand range of targets 

## Support

Contact brittney.cooper@noirlab.edu with questions, problems, suggestions.


