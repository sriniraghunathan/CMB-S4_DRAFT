

# CMB-S4: Dark Radiation Anisotropy Flowdown Team (DRAFT) tool

Overview:
===============================================
* Optimally combines data from different bands using noise + foreground signals in different bands.
  * Supports standard / constrained / partial internal linear combinations.
* Computes lensing noise curves using residual noise.
* Combines delensed CMB spectra and lensing spectra to forecast cosmological parameter constraints using Fisher formalism.
  * Delensing reference: Green, Meyers, van Engelen 2016, arXiv: [1609.08143](https://arxiv.org/abs/1609.08143).
  * Delensing + Fisher code reference: Hotinli, Meyers, Trendafilova, Green, van Engelen 2021, arXiv: [2111.15036](https://arxiv.org/abs/2111.15036).
    *  [Code repo](https://github.com/ctrendafilova/FisherLens).
  * [Simple Fisher code repo without delensing](https://github.com/sriniraghunathan/cmbs4_fisher_forecasting).
* Estimates biases in cosmological parameters due to residual foregrounds also using Fisher formalism **(yet to be integrated into this repo)**.
  * References: Huterer & Takada 2004, arXiv: [0412142](https://arxiv.org/abs/astro-ph/0412142); Loverde, Hui, & Gaztanaga 2006, arXiv: [0611539](https://arxiv.org/abs/astro-ph/0611539); Amara & Réfrégier 2007, arXiv: [0710.5171](https://arxiv.org/abs/0710.5171).
---------

Sub-modules and Dependencies:
===============================================
1. CLASS_delens
https://github.com/selimhotinli/class_delens (9944e0ec5f0b734b178436a8652b186046f63fdb)
 
  This code uses a wrapper for the CLASS_delens code to facilitate Fisher forecasting of cosmological parameter constraints from CMB spectra.

  ### Authors: 
  * Selim C. Hotinli, Joel Meyers, Cynthia Trendafilova, Daniel Green, Alexander van Engelen
---------

Steps to run the code:
===============================================
 ### Cloning:
  1. Clone the repo first: ```git clone git@github.com:sriniraghunathan/CMB-S4_DRAFT.git```
  2. list the submodules: ```git submodule```. This will show the ```CLASS_delens``` and the specific commit ```9944e0ec5f0b734b178436a8652b186046f63fdb```
  3. pull the submodule: ```git submodule update --init --recursive```
 ### CLASS installation
  1. Install CLASS: ```cd CLASS_delens; make class```
  2. Ensure CLASS is working properly by running the following document within the ```CLASS_delens``` folder: ```./class explanatory.ini```
---------

Contributors:
=============================================== 
[_Joel_ **Meyers**](https://joelmeyers.github.io/), [_Cynthia_ **Trendafilova**](https://github.com/ctrendafilova), and [_Benjamin_ **Wallisch**](https://www.ias.edu/scholars/benjamin-wallisch).

### Delensing and Fisher formalism contributions: 
* Selim C. Hotinli, Joel Meyers, Cynthia Trendafilova, Daniel Green, Alexander van Engelen
---------

CMB-S4 instrument/noise specs:
===============================================
* Chilean LAT: [S4-Wide](https://cmb-s4.uchicago.edu/wiki/index.php/Expected_Survey_Performance_for_Science_Forecasting#Instrument_Definition)
* Delensing SouthPole LAT: [S4-Ultra deep](https://cmb-s4.uchicago.edu/wiki/index.php/Delensing_sensitivity_-_updated_sensitivities,_beams,_TT_noise)
---------

Foreground modelling:
===============================================
* **Extragalactic foregrounds**: Radio, CIB, tSZ and  kSZ power spectra from SPT measurements (George et al. 2015, arXiv: [1408.3161](https://arxiv.org/abs/1408.3161) and Reichardt et al. 2020, arXiv: [2002.06197](https://arxiv.org/abs/2002.06197)).
  * Assumed polarisation fractions: CIB = 2%; Radio = 3%; tSZ/kSZ = 0. But these are configurable. Look into [params.ini](https://github.com/sriniraghunathan/DRAFT/blob/master/scripts/notebooks/params.ini).
* **Galactic foregrounds**: Dust and Synchrotron power spectra obtained from [pySM3](https://github.com/CMB-S4/s4mapbasedsims/tree/master/202002_foregrounds_extragalactic_cmb_tophat) simulations.
---------

Results:
===============================================
* **ILC curves**:
  * Look into this [link](https://github.com/sriniraghunathan/DRAFT/tree/master/results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/). 
     * Look into *read_file.py* script to read ILC curves. 
     * Standard ILC curves for S4-Wide (PBDR configuration) will be under this [link](https://github.com/sriniraghunathan/DRAFT/tree/master/results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4like_mask_v2/TT-EE/baseline). 
     * Constrained ILC curves (for galactc dust) are also in the above folder. 
---------
