# deep-convection-initiation_code
This repository highlights some of the algorithms I developed for my master's project on the sensitivity of deep convection initiation (DCI) to horizontal grid resolution

**Background:**\
Moist convection is fundamental in regulating a variety of atmospheric phenomena ranging from global-scale circulations to precipitating cumulus clouds with sub-kilometre horizontal length scales. Numerical simlation of this process exhibits strong sensitivity to model grid resolution, primarily because turbulence and microphysical processes, which are integral to convective clouds, are inadequately resolved on model grids, leading to substantial errors and biases in DCI forecasts. The focus of this project is placed on addressing two important objectives: (i) to thoroughly evaluate the sensitivity of simulated deep-convection initiation over a mesoscale convergence to horizontal grid spacing and (ii) to physically interpret the cloud-layer and subcloud processes regulating this sensitivity.

**Model setup:**\
Numerical simulations are conducted using Bryan Cloud Model (CM1) version 19, with domain size of 120 km, 180 km, and 20 km in the x, y, and z directions, respectively. While the vertical grid size increases with height, the horizontal grid size is varied from 125 m, to 250 m, to 500 m, and to 1 km in different experiments. A diurnally time-varying and Gaussian-shaped surface heating function with peak amplitude located at the center of the domain (x=0 km) is imposed to induce mesoscale convergence in the x-z plane. The heating function is symmetric in x and uniform in y. The former generates locally strong horizontal convergence and a subcloud updraft that breaches into the cloud layer to initiate cumulus convection, while the latter allows statisitcal sampling of the numerically simulated clouds.

**Analyses:**\
There are two scripts presented for each analysis, one for analyzing numerical output and the other for making plots using the processed data. Model output is written to file (.nc) every simulated minute (60s) or 90s, depending on cases. We used Matlab to read variables from the NetCDF files. 

(1) Cloud-top height:
creating a time series that shows the maximum height of cloud over the course of the model integration time (12hrs), with different grid spacings used.
scripts are 'cloudth.m' and 'cloudth_1_plot.m'

(2) Surface convergence:
creating a time series which shows the magnitude of surface convergence, averaged in y along the convergence line, and in x over the central 1 km of the domain.
scripts are 'conver_2.m' and 'subplot_4fig_1.m'

(3) Mass flux profile:








**Ref.**\
[Details of CM1 model](https://www2.mmm.ucar.edu/people/bryan/cm1/)
