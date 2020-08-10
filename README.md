# deep-convection-initiation_code
This repository highlights some of the algorithms I developed for my master's project on the sensitivity of deep convection initiation (DCI) to horizontal grid resolution

**Background:**\
Moist convection is fundamental in regulating a variety of atmospheric phenomena ranging from global-scale circulations to precipitating cumulus clouds with sub-kilometre horizontal length scales. Numerical simlation of this process exhibits strong sensitivity to model grid resolution, primarily because turbulence and microphysical processes, which are integral to convective clouds, are inadequately resolved on model grids, leading to substantial errors and biases in DCI forecasts. The focus of this project is placed on addressing two important objectives: (i) to thoroughly evaluate the sensitivity of simulated deep-convection initiation over a mesoscale convergence to horizontal grid spacing and (ii) to physically interpret the cloud-layer and subcloud processes regulating this sensitivity.

**Model setup:**\
Numerical simulations are conducted using Bryan Cloud Model ([CM1](https://www2.mmm.ucar.edu/people/bryan/cm1/)) version 19, with domain size of 120 km, 180 km, and 20 km in the x, y, and z directions, respectively. While the vertical grid size increases with height, the horizontal grid size is systematically varied from 125 m, to 250 m, to 500 m, and to 1 km in different experiments. A diurnally time-varying and Gaussian-shaped surface heating function with peak amplitude located at the center of the domain (x=0 km) is imposed to induce mesoscale convergence in the x-z plane. The heating function is symmetric in x and uniform in y. The former generates locally strong horizontal convergence and a subcloud updraft that breaches into the cloud layer to initiate cumulus convection, while the latter allows statisitcal sampling of the numerically simulated clouds.

**Analyses:**\
There are generally two scripts presented for each analysis, one for analyzing numerical output and the other for making plots using the processed data. Model output is written to file (.nc) every simulated minute (60s) or 90s, depending on cases. Variables from the NetCDF files are read using Matlab. Each analysis compares the behaviour or amplitude of the stated quantity among the different grid spacings used. 


(I) Cloud-layer analysis
  1. Cloud-top height:
  * a time series that shows the maximum height of cloud over the course of the model integration time (12hrs), with different grid spacings used.  <br/>
  * related scripts: 'cloudth.m' and 'cloudth_1_plot.m'
  
  2. Surface convergence:
  * a time series which shows the magnitude of surface convergence, averaged in y along the convergence line, and in x over the central 1 km of the domain.  <br/>
  * related scripts: 'conver_2.m' and 'subplot_4fig_1.m'
  
  3. Mass flux profile
  * creates profiles of (a) total mass flux, (b) total area, (c)average vertical velocity w, average buoyancy b, and average hydrometero mixing ratio qh, of convective cores over the whole the domain during each simulated hour.
  * related scripts: 'ascmsfx_wbpqcqi_hydrometeor.m' and 'ascmf_wbpqcqi_hydrometeor_plot.m'
  
  4. Fractional entrainment and detrainment profiles:
  * uses ice-liquid water potential temperature as conserved variable when calculating moist static energy (see [Bryan and Fritsch (2004)](https://doi.org/10.1175/1520-0493(2004)132%3C2421:AROIWP%3E2.0.CO;2) for details)
  * create profiles of (a) fractional entrianment rate and (b) fractional detrainment rate, averaged over all convective clouds during each hour.
  * related scripts: (analysis): 'detrainment_theta_il.m', 'fracmf.m', 'detrainment_theta_il.m'; (plotting): ' entrain_detrain_thil_plot.m'
  
(II) Parcel trajectories:
  1. Distribution of vertical velocity of parcels prior to crossing cloud base
  * uses parcel trajectories provided by the model to monitor the spatial and temporal evolution of parcels
  * measures the numbers of parcels reaching different heights
  * creates histograms showing (a) the distribution of vertical velocity of parcels prior to crossing cloud base for the last time, (b) percentages of parcels that are able to reach different heights
  * related scripts: 'prcl_track_wlcl_1.m' and'prcl_track_wlcl_plot.m'

(III) Subcloud-layer analysis:
  1. Vertical velocity (w)-budget profile
  * computes each term in the vertical momentum euqation
  * creates profiles that show (a) left-hand side term and (b) sum of right-hand side terms, of the equation  
  * related scripts: 'wbudget_1_profiles.m' and 'wbudget_1_profiles.m'
  
  2. Contour plot of each term in the w-momentum equation
  * computes y- and time-average of each term in the vertical momentum equation
  * creates contour of each term and makes comparions of a specific term between different cases (e.g. 125 m vs 250 m, and 500 m vs 250 m)
  * related scripts: 'wbudget_1.m' and 'wbudget_1_plot.m'
  
  3. Profile of total turbulent kinetic energy (TKE)
  * explicitly computes perturbation of u, v, and w from their respective y-averages,
  * then calculates the resolved component of TKE
  * get total TKE by combining the resolved component calculated and the subgrid component provided by the model
  * related scripts: 'myuw_exnstar_tke.m' and 'myuw_exnstar_tke_plot.m'
  
  4. Contours of total (subgrid + resolved) turbulent forcing on vertical velocity field
  * calculates perturbation velocity (e.g. u', v', and w') with respect to the y-average, and turbulent fluxes u'w', u'w', and w'w'
  * then computes resolved component of turbulent forcing on the w-field according to Reynold's averaged  momentum equations
  * gets subgrid component of turbulent forcing on w form model output
  * obtains total (subgrid + resolved) turbulent foricng on w
  * creates contours of total turbulent foricng on w and compares them between different cases (e.g. 125 m vs 250 m and 500 m vs 250 m)
  * related scripts: 'turbflux_w_2.m' and 'turbflux_w_2_plot.m'
