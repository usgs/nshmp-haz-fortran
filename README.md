## National Seismic Hazard Model (NSHM) Fortran Codes

This repository houses the codes used to generate the [2008](http://pubs.usgs.gov/of/2008/1128/) and [2014](http://pubs.usgs.gov/of/2014/1091/) updates to the NSHM for the conterminous US. The 2008 codes are tagged as [nshm2008r3](https://github.com/usgs/nshmp-haz-fortran/tree/nshm2008r3). The 2014 codes are tagged as [nshm2014r1](https://github.com/usgs/nshmp-haz-fortran/tree/nshm2014r1). These codes include all necessary configuration and data files and were compiled and run using [GFortran](http://gcc.gnu.org/fortran/) to produce the maps and data currently available on the [USGS website](http://earthquake.usgs.gov/hazards/products/conterminous/).

*__Note:__ At this time, these codes are provided with limited documentation and support. We anticipate releasing a unified code base with more thorough support later this fall. See notes on 2014 update, below.*

#### Notes on the 2014 update (IMPORTANT)

The California component of the 2014 update was based on the 3rd version of the Uniform California Earthquake Rupture Forecast ([UCERF3](http://pubs.usgs.gov/of/2013/1165/)). This model was developed using [OpenSHA](http://www.opensha.org). The contributions to hazard from UCERF3 sources were also calculated using OpenSHA; these were then added to the contributions from other Western US sources calculated using the 2014 Fortran codes. Replicating the hazard curves for a site such as Reno, NV is difficult at this time as it requires running both codes and summing the results.

For sites in California, the UCERF3 model may be accessed via the nightly builds of OpenSHA applications (e.g. the [Hazard Curve Calculator](http://www.opensha.org/apps-HazardCurveLocal)) or by checking out the repository and running the code directly from a development environment such as Eclipse.

The UCERF3 source model comes in a variety of flavors, all of which are collectively referred to as [Fault System Solutions](http://opensha.usc.edu/trac/wiki/UCERF3FaultSystemSolutions) (see [this page](http://opensha.usc.edu/trac/wiki/UCERF3FaultSystemSolutions) for detailed information). Hazard curves for the 2014 NSHM were computed using the "True Mean" solution. Although OpenSHA provides some facility to downsample (via averaging of rupture parameters) this model, be aware that it takes a large amount of memory (generally > 8GB) to load and "Out of Memory" errors are not uncommon for users without the necessary hardware.

