# EHTW

The original atmospheric model is taken from ECHAM 5.4 (e.g. Rockner et al., 2006) and the original ocean model is taken from the dynamics of TIMCOM (Tseng & Chien,2011; Young et al., 2012). Climate model ECHAM5 has been developed from the ECMWF operational forecast model cycle 36 (EC, 1989) and a comprehensive parameterisation package developed at Hamburg (HAM). Dynamics of ECHAM is based on the ECMWF, modified to describe the newly implemented features and the changes necessary for climate experiments. After the release of the ECHAM4, the whole source code has been extensively redesigned in the major infrastructure and transferred to Fortran 95. ECHAM is now fully portable and runs on all major high performance platforms. The restart mechanism is implemented on top of netCDF and because of that absolutelly independent on the underlying architecture. The dynamical part of ECHAM is formulated in spherical harmonics. After the inter-model comparisons by Jarraud et al. (1981) and Girard and Jarraud (1982) truncated expansions in terms of spherical harmonics were adopted for the representation of dynamical fields. The transform technique developed by Eliasen et al. (1970), Orszag (1970), and Machenhauer and Rasmussen (1972) is used such that non-linear terms, including parameterizations, are evaluated at a set of almost regularly distributed grid points - the Gaussian grid. In the vertical, a flexible coordinate is used, enabling the model to use either the usual terrain following sigma coordinate (Phillips, 1957), or a hybrid coordinate for which upper-level model surfaces atten over steep terrain, becoming surfaces of constant pressure in the stratosphere (Simmons and Burridge (1981) and Simmons and Struring (1981). Moist processes are treated in a different way using a mass conserving algorithm for the transport (Lin and Rood, 1996) of the different water species and potential chemical tracers. The transport is determined on the Gaussian grid.

TIMCOM is an advanced version derived from DIECAST(e.g., Tseng et al., 2005). It solves the three-dimensional primitive equations for fluid motions on the spherical or Cartesian coordinate under hydrostatic and Boussinesq approximations. The TIMCOM uses a Z-level stretched vertical coordinate and a horizontal grid based on spherical coordinates with non-uniform latitudinal increments. Either rigid-lid or free-surface option is available. A standard nonlinear equation of state is used to calculate in-situ density in terms of potential temperature, salinity and pressure-depth. Weakly physical-based parameterization is used for the horizontal mixing. TIMCOM uses a blend of collocated and staggered grid structure with non-uniform grid increments. Second-order Leap-frog temporal scheme and fourth-order-accurate spatial discretization (or interpolation) are applied to discretize the governing equations. Efficient and accurate Error Vector Propagation (EVP) solver is used to solve the system matrix. EVP solver is also ideal for parallel computing, see Tseng et al. (2012), Yiung et al. (2012) for more model details.

The major significant model development in our current implementation is its inclusion of SIT solver (Tsuang et al., 2001). The solver has been improved for ocean cool skin simulation (Tu and Tsuang, 2005) and better formulation (Tsuang et al., 2009), and has been applied for South China Sea SST simulation (Lan et al., 2010). SIT denotes for Snow/Ice/Thermocline solver. SIT solves the vertical heat diffusion equation for temperatures in snow, ice and water column in a tri-diagonal matrix. This solver is numerically unconditional stable. Hence, melting and refreeze of a thin sea ice, as well as the warm layer and the cool skin in the upper few meters of a water column can be resolved. This can give a unique surface diurnal cycle of ocean within the climate models. Only a few climate models have explored the importance of diurnal cycle. Second, our ATM and OCN coupling frequency is every ATM time step. This is unique and to the best of our knowledge, no other climate model has used this approach because it is too expensive. The ECHAM/SIT/TIMCOM does not use any specific coupler while the OCN is directly coupled with ECHAM through SIT. So, the daily cycle of the ocean surface can easily be resolved. Finally, we don’t use any fancy interpolation scheme while our ocean grid collocates with the atmospheric grid (or its subdivision equally). So this guarantees the surface fluxes are conserved without using any specific technique. These are the superior points of the current model architecture.

We should note that SIT is not only a pseudo atmosphere-ocean coupler but also a “simplified” 1-D ocean model which can have very high resolution in the ocean top 10 m or so. The average of 0-10m SIT result is then feeding into the ocean model as the OCN boundary condition. Since it is a simplified 1-D ocean model, SIT can do more than a coupler. Our coupled simulation shall have some characteristics resolving high frequency atmospheric-oceanic dynamics.

ECHAM5/SIT：

It is a coupled AGCM (ECHAM5)(e.g. Rockner et al., 2006) with a one-column ocean model (SIT)(Tsuang et al., 2001). The original atmospheric model is taken from ECHAM 5.4. SIT denotes for Snow/Ice/Thermocline solver. SIT solves the vertical heat diffusion equation for temperatures in snow, ice and water column in a tri-diagonal matrix. This solver is numerically unconditional stable. Hence, melting and refreeze of a thin sea ice, as well as the warm layer and the cool skin in the upper few meters of a water column can be resolved. This can give a unique surface diurnal cycle of ocean within the climate models. The solver has been improved for ocean cool skin simulation (Tu and Tusang, 2005) and better formulation for skin temperature simulations of sea surface, ice-sheet surface and snow surface (Tsuang et al., 2009), and has been applied for South China Sea SST simulation (Lan et al., 2010).

Reference：

Jarraud, M., Girrard, C. and Geleyn, J.-F. (1982): Note on a possible linearization of the vorticity equation in a primitive spectral model. Research Activities in Atmospheric and Ocean Modelling, Report No. 3, Working Group on Numerical Experimentation, Geneva.

Jarraud, M., Girrard, C. and Cubasch, U. (1981): Comparison of medium-range forecasts made with the models using spectral or finite-difference techniques in the horizontal. Technical Report 23, ECMWF.

Eliasen, E., Machenhauer, B. and Rasmussen, E. (1970): On a numerical method for integration of the hydrodynamical equations with a spectral representation of the horizontal fields. Tech. Rep. 2, Institute of Theoretical Meteorology, University of Copenhagen.

Orszag, S. A. (1970): Transform method for calculation of vector coupled sums. J. Atmos. Sci., 27, 890-895.

Machenhauer, B. and Rasmussen, E. (1972): On the integration of the spectral hydrodynamical equations by a transform method. Tech. Rep. 4, Institute of Theoretical Meteorology, University of Copenhagen.

Phillips, N. A. (1957): A coordinate system having some special advantages for numerical forecasting. Journal of Meteorology, 14, 184-185.

Simmons, A. J. and Burridge, D. M. (1981): An energy and angular-momentum conserving vertical finite difference scheme and hybrid vertical coordinates. Mon. Wea. Rev., 109, 758-866.

Simmons, A. J. and Strufing, R. (1981): An energy and angular-momentum conserving finite difference scheme, hybrid coordinates and medium-range weather prediction. Technical Report 28, ECMWF, Reading, UK.

Lin, S. J. and Rood, R. B. (1996): Multidimensional flux form semi-Largrangian transport. Mon. Wea. Rev., 124, 2046-2068.

Roeckner, E., R. Brokopf, M. Esch, M. Giorgetta, S. Hagemann, L. Kornblueh, E. Manzini, U. Schlese, and U. Schulzweida (2006): Sensitivity of simulated climate to horizontal and vertical resolution in the ECHAM5 atmosphere model, J. Climate, 19, 3771-3791.

Lan, YY; Tsuang, BJ; Tu, CY; Wu, TY; Chen, YL; Hsieh, CI, 2010/04, “Observation and Simulation of Meteorology and Surface Energy Components over the South China Sea in Summers of 2004 and 2006”, TERRESTRIAL ATMOSPHERIC AND OCEANIC SCIENCES21 (2): 325-342

Tseng, Wan-Ling (2012) Role of ocean-atmosphere interaction for tropical climate variability over warm pool regions (Doctoral thesis/PhD), Christian-Albrechts-Universität Kiel, Kiel, Germany, 78 pp. (http://eprints.uni-kiel.de/19781/)

Tseng, Y. H., Dietrich, D. E., & Ferziger, J. H. 2005. Regional circulation of the Monterey Bay region- Hydrostatic versus Non-hydrostatic modeling. J. Geophys. Res., 110, doi:10.1029/2003JC002153

Tseng, Y.-H. and Chien, M.-H., 2011: Parallel Domain-decomposed Taiwan Multi-scale Community Ocean Model (PD-TIMCOM), Computers & Fluids, 45, 77-83. (http://efdl.as.ntu.edu.tw/research/timcom/index.html)

Tseng, Y. H., Shen, M.L., Jan, S., Dietrich, D.E. and Chiang, C.P.(2012), 'Validation of the Kuroshio current system in the dual-domain Pacific Ocean model framework', Progress in Oceanography doi:10.1016/j.pocean.2012.04.003.

Tsuang, B.-J.; Tu, C.-Y.; Arpe K.; 2001/03: Lake Parameterization for climate models. Max-Planck-Institute for Meteorology, Report No. 316, Hamburg. 72 pp.

Tsuang, B.-J., C.-Y. Tu, J.-L. Tsai, J.A. Dracup, K. Arpe and T. Meyers, 2009: A more accurate scheme for calculating Earth’s skin temperature. Climate Dynamics 32, 251-272, DOI 10.1007/s00382-008-0479-2.

Tu, C.-Y.; Tsuang, B.-J., 2005/11: Cool-skin simulation by a one-column ocean model, Geophys. Res. Lett., 32, L22602, doi:10.1029/2005GL024252.

Young, C.C., Tseng, Y.H., Shen, M.L., Liang, Y.C., Chien, M.H. and Chien, C.H.(2012), 'Software development of the Taiwan Multi-scale Community Ocean Model (TIMCOM)', Environmental Modelling & Software, 38, 214-219. doi: 10.1016/j.envsoft.2012.05.017
