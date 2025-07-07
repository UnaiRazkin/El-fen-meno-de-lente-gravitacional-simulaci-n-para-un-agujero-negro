http://archive.stsci.edu/prepds/frontier/lensmodels/
http://archive.stsci.edu/pub/hlsp/frontier/abell370/models/
hlsp_frontier_model_abell370_v1_readme.txt

This directory contains gravitational lensing models of the Hubble Frontier Fields (HFF) cluster Abell 370 based on pre-HFF imaging.  The models were produced by five map making teams who shared observational constraints but then worked independently to derive these models using robust, established methodologies.

These models may also be interactively displayed:
http://archive.stsci.edu/prepds/frontier/lensmodels/

And magnification estimates at a given position and redshift (RA, Dec, z) may be extracted using this web tool:
http://archive.stsci.edu/prepds/frontier/lensmodels/webtool/magnif.html

Eight lens models were submitted for Abell 370.  In addition to models from the Bradac, CATS, Sharon, and Williams teams, the Merten-Zitrin team submitted four models.  More information about the methodologies may be found in the README files in each team's model subdirectory as well as at this webpage:
http://www.stsci.edu/hst/campaigns/frontier-fields/Lensing-Models

The following lensing primer helps explain the contents of the model directories:
http://archive.stsci.edu/prepds/frontier/lensmodels/webtool/hlsp_frontier_model_lensing_primer.pdf

Each team's model subdirectory contains FITS images with the following data:

- z01-magnif: magnification for a lensed galaxy at z = 1
- z02-magnif: magnification for a lensed galaxy at z = 2
- z04-magnif: magnification for a lensed galaxy at z = 4
- z09-magnif: magnification for a lensed galaxy at z = 9

and the following scaled to DLS / DS = 1 (see lensing primer):

- kappa: mass surface density 
- gamma: weak lensing shear

Based on these mass and shear maps, magnification maps may be calculated for any redshift using the Python script provided in this directory, or as described in the lensing primer.

Some teams also provided the following (again scaled to DLS / DS = 1):

- gamma1: weak lensing shear component 1
- gamma2: weak lensing shear component 2

- x-arcsec-deflect: image deflection (from source to image plane) along the x-axis [in arcseconds]
- y-arcsec-deflect: image deflection (from source to image plane) along the y-axis [in arcseconds]
- x-pixels-deflect: image deflection (from source to image plane) along the x-axis [in pixels]
- y-pixels-deflect: image deflection (from source to image plane) along the y-axis [in pixels]

Some teams also provided files containing model parameters and/or object catalogs as described in their README files.

All teams also provided a range of possible models as constrained by pre-HFF lensing observables.  Mass and shear maps are given as numbered files in subdirectories under each team:

- range/*map###*.fits

From this range of models, we may calculate a range of possible magnification maps at any redshift (as described above) according to each method.

To download all of the maps in these directories, you may use wget commands such as:
  wget -nH --cut-dirs=7 -r -l0 -c -N -np -R 'index*' -erobots=off http://archive.stsci.edu/pub/hlsp/frontier/abell370/models/cats/range/

To view these images in ds9, we recommend scalings of:
- magnif: log 1 100
- kappa: linear 0 3
- shear: linear 0 1

The models cover the following areas with the following resolutions as given in the FITS WCS headers:

model   |  center RA & Dec (J2000)   |  width  | resolution | pixels on a side
CATS      02:39:53.129  -01:34:26.58    4.003'     0.2002"     1200
Sharon    02:39:52.948  -01:34:36.75    3.368'     0.0300"     6736
Zitrin    02:39:52.906  -01:34:37.82    2.500'     0.0500"     3000
Williams  02:39:53.027  -01:34:36.90    2.328'     0.2794"      500
Bradac    02:39:54.697  -01:34:19.10    3.000'     0.0439"     4096
Merten    02:39:53.075  -01:34:56.10   25.002'     6.2500"      240

If you have any questions about these data products, please e-mail Dan Coe <DCoe@STScI.edu>.
