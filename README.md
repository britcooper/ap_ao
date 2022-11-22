# ap_ao

Ap_AO estimates the feasbility of potential extended object observations with Gemini North's Laser Guide Star (LGS) and adaptive optics unit (Altair). 


## Usage

To assess feasibilty, Ap_AO performs background subtraction and sums the flux within Altair's 
LGS tip/tilt aperture using an image of the potential target observed above the atmosphere.

Currently, the image takes the form of a HST WFPC2 "c0m" fits file, available on the [HST archive](https://mast.stsci.edu/search) .

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


## Dependencies

TBD

## Potential Issues

TBD

## TBD

TBD

## Support

Contact brittney.cooper@noirlab.edu with questions or problems.


