# BXA-Plasma
Log-norm temp distribution of APEC plasmas fitted with BXA.

BXA connects the X-ray spectral analysis environments Xspec/Sherpa to the nested sampling algorithm 'UltraNest' for Bayesian Parameter Estimation and Model comparison. 
It was designed to deal with systematically analysing a large data-set, comparing multiple models and analysing low counts data-set with realistic models.

The BXA-Plasma combine the BXA framework embedded with the PCA-based background model to connect the Astrophysical Plasma Emission Code (APEC, Smith et al. 2001) or VAPEC which allow the variant of the abundances. 
## Method/how it works
See how BXA works here: http://johannesbuchner.github.io/BXA/index.html

## Install 
You will need to install 'UltraNest' see documentation here: https://johannesbuchner.github.io/UltraNest/installation.html

## How to run
For this work, I worked in Sherpa environment and Chandra data. 

To load BXA in a session with sherpa:

import bxa.sherpa as bxa

from bxa.sherpa.background.pca import auto_background

from bxa.sherpa.background.models import ChandraBackground

from bxa.sherpa.background.fitters import SingleFitter

## How to report issues

## How to cite
If you make use of these scripts, please cite this work: https://arxiv.org/abs/2401.17303

## License

see LICENCE file

