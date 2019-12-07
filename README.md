#  fvcom

[![fvcom web site](http://www.easyicon.net/api/resizeApi.php?id=1205810&size=24)](http://fvcom.smast.umassd.edu)

###  Description

FVCOM is a prognostic, unstructured-grid, finite-volume, free-surface, 3-D primitive equation coastal ocean circulation model developed by UMASSD-WHOI joint efforts. The model consists of momentum, continuity, temperature, salinity and density equations and is closed physically and mathematically using turbulence closure submodels. The horizontal grid is comprised of unstructured triangular cells and the irregular bottom is preseented using generalized terrain-following coordinates. The General Ocean Turbulent Model (GOTM) developed by Burchard’s research group in Germany (Burchard, 2002) has been added to FVCOM to provide optional vertical turbulent closure schemes. FVCOM is solved numerically by a second-order accurate discrete flux calculation in the integral form of the governing equations over an unstructured triangular grid. This approach combines the best features of finite-element methods (grid flexibility) and finite-difference methods (numerical efficiency and code simplicity) and provides a much better numerical representation of both local and global momentum, mass, salt, heat, and tracer conservation.  The ability of FVCOM to accurately solve scalar conservation equations in addition to the topological flexibility provided by unstructured meshes and the simplicity of the coding structure has make FVCOM ideally suited for many coastal and interdisciplinary scientific applications.  
<br />
FVCOM was originally developed for the estuarine flooding/drying process in estuaries and the tidal-, buoyancy- and wind-driven circulation in the coastal region featured with complex irregular geometry and steep bottom topography. This model has been upgraded to the spherical coordinate system for basin and global applications. A non-hydrostatic version of FVCOM has been coded and is being tested.  
<br />
![fvcom](http://fvcom.smast.umassd.edu/wp-content/uploads/2013/11/fvcom1.jpg)  
<br />
The present version of FVCOM includes a number of options and components as shown in Figure above. These include (1) choice of Cartesian or spherical coordinate system, (2) a mass-conservative wet/dry point treatment for the flooding/drying process simulation, (3) the General Ocean Turbulent Model (GOTM) modules (Burchard et al., 1999; Burchard, 2002) for optional vertical turbulent mixing schemes, (4) a water quality module to simulate dissolved oxygen and other environmental indicators, (5) 4-D nudging and Reduced/Ensemble Kalman Filters (implemented in collaboration with P. Rizzoli at MIT) for data assimilation, (6) fully-nonlinear ice models (implemented by  F. Dupont), (7) a 3-D sediment transport module (based on the U.S.G.S. national sediment transport model) for estuarine and near-shore applications, and (8) a flexible biological module (FBM) for food web dynamics study. FBM includes seven groups: nutrients, autotrophy, heterotrophy, detritus, dissolved organic matter, bacteria, and other. With various pre-built functions and parameters for these groups, GBM allows users to either select a pre-built biological model (such as NPZ, NPZD, etc.) or to build their own biological model using the pre-defined pool of biological variables and parameterization functions. FVCOM was originally coded for sigma-coordinates in the vertical and now has been upgraded to a generalized terrian-following coordinate system with choices of various topographic-following coordinates. FVCOM is written with Fortran 90 with MPI parallelization, and runs efficiently on single and multi-processor machines.  
<br />
FVCOM is an open source code ocean community model that always welcomes new users. This program is only permitted for use in non-commercial academic research and education. Users are required to register in orde to receive the source codes, demo examples, and user manuals as well as some recommended postprocessing tools.  
<br />
A detailed description of FVCOM is given in ”[User Manual](http://fvcom.smast.umassd.edu/wp-content/uploads/2013/11/MITSG_12-25.pdf)” and published papers.  

###  Repos
fvcom 2.7.1  
fvcom 3.2  
fvcom 4.0  
fvcom 4.3
