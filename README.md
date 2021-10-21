# TPPSS (TopopoSunsun)

Computes sunrise / sunset times taking into account local topography provided by the user as a DEM (it can also simply compute an horizon based on the DEM). Provides a library and a command line script to compute sunrise / sunset for a location for a single day or a whole year.

# Install

*Not yet on PyPI*

TODO mention rasterio + conda

The tool requires Python 3.7+.

To install, launch :

```console
pip install tppss
```

The command above will install the `tppss` Python library and its dependencies. The library includes a command-line script, also named `tppss`, whose functionality is described below.

Library documentation is avaiable (here)[http://example.com]. (*not yet*)

# Instructions

# Acknowledgements

The horizon computation is based on this Matlab code (adapted for Python + Numpy):

Benjamin Pillot (2021). DEM-based topography horizon model 
(https://www.mathworks.com/matlabcentral/fileexchange/59421-dem-based-topography-horizon-model),
MATLAB Central File Exchange. Retrieved October 19, 2021. 


The sun position computation is adapted from this code:

John Clark Craig. Python Sun Position for Solar Energy and Research
(https://levelup.gitconnected.com/python-sun-position-for-solar-energy-and-research-7a4ead801777)


# TODO 

- Add CLI + CSV export for year
- Projected DEM
- Example that draws horizon + sun course through the sky
- Generate doc
- Sample: Document how the input DEM is obtained
