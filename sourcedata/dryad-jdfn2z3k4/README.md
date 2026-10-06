# Data from: Neural interactions in the human frontal cortex dissociate reward and punishment learning

[https://doi.org/10.5061/dryad.jdfn2z3k4](https://doi.org/10.5061/dryad.jdfn2z3k4)

Preprocessed data containing intracranial EEG high-gamma activity during a probabilistic learning task

## Data organisation

The data are organized into three folders :

\- **power :** netcdf files of gamma power of individual subjects. Each matrix has a shape of `(n_{trials}, n_{roi}, n_{times})`. When loading this file using Xarray, each dimension of the matrix contains the labels (i.e. whether it's a rewarding or punishing trial, the channel name and time points).

\- **behavior :** behavioral xlsx of length `(n_{trials},)`. This table contains the values of the prediction error during the rewarding and punishing conditions.

\- **anatomy :** xlsx file of length `(n_{roi},)` containing anatomical information (i.e. the name of the brain region) of each recording contact.

To load the gamma power using Python we recommend using the Xarray python :

\```python

import xarray as xr

data = xr.load_dataarray(".../power_subject-X.nc")

### Code/Software

The Python scripts and notebooks to reproduce the results presented here are hosted on Github

([https://github.com/brainets/papercode/tree/main/combrisson_pblt_2024](https://github.com/brainets/papercode/tree/main/combrisson_pblt_2024)).
