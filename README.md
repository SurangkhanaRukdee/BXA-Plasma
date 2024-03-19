# BXA-Plasma
Log-norm temp distribution of APEC plasmas fitted with Bayesian X-ray Analysis (BXA).

BXA connects the X-ray spectral analysis environments Xspec/Sherpa to the nested sampling algorithm 'UltraNest' for Bayesian Parameter Estimation and Model comparison. 
It was designed to deal with systematically analysing a large data-set, comparing multiple models and analysing low counts data-set with realistic models.

The BXA-Plasma combine the BXA framework embedded with the PCA-based background model to connect the Astrophysical Plasma Emission Code (APEC, Smith et al. 2001) or VAPEC which allow the variant of the abundances. 

Here we introduced multi temperature distributions instead of using one or two temperature (grid). The temperature distributions resemble a log-Gaussian during quiescent phases, and also during flares, which
typically exhibit a broader and hotter temperature distribution (Robrade & Schmitt 2005). Therefore, we model the spectra with a plasma featuring a log-Gaussian temperature distribution. This approach is expected to capture the behavior of the plasma temperature better than a single point (kT1 or kT2). In practice, this is achieved by an ensemble of 10 APEC components on a logarithmic temperature grid, where the normalizations follow a bell curve as part of the modeling process in BXA.

![multiT_model](https://github.com/SurangkhanaRukdee/BXA-Plasma/assets/9215336/f440c9f4-e038-4e7a-b062-90c814cce3bb)

The left image is the cartoon model for this multi temperature distribution idea while on the right image is the implementation example of multi temperature distributions (in various colors) during the quasi-quiescence stage. The average of these distributions are plotted in black color with the gray area of uncertainty.

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

