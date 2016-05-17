# CESM_1.2.2_Subgrid_Triggering
This repository contains source modules necessary for running the Community Earth System Model (CESM; http://www.cesm.ucar.edu/models/cesm1.2/) version 1.2.2 with a new sub-grid deep convective trigger.  Note that these source modules are met to be run within the CESM framework for the specific release version, 1.2.2.  The scientific description and implications of including the sub-grid trigger are described in paper entitled "Representing sub-grid convective initiation in the CESM" by Tawfik, Dirmeyer, and Lawrence (2016). Here are steps for implementing the new source modules.

1)  After registering and downloading the CESM 1.2.2, create a new case by invoking the following in the command line:
> cd cesm1_2_2/scripts/
> ./create_newcase -case your_case_name -res your_resolution -compset your_component_set -mach your_machine_name

where your_case_name = whatever you want to call your case, your_resolution = any of the supported grids, your_component_set = any of the supported component sets that include both atmosphere and land models functioning for the sub-grid to work, your_machine_name = the specific machine name.  To recreate the work found in Tawfik et al. (2016) here is the case:
> ./create_newcase -case BCL_TILE_SE_CAM5_1.00 -res a%ne30np4_l%ne30np4_oi%gx1v6_r%r05_m%gx1v6_g%null_w%null -compset FAMIPC5 -mach yellowstone

2)  Once the case is created a directory with the name your_case_name will appear.  Now we can setup and build the code while including the new sub-grid source modules.
> cd your_case_name
> ./cesm_setup
> cp src.cam/* SourceMods/src.cam/
> cp src.clm/* SourceMods/src.clm/
> cp src.drv/* SourceMods/src.drv/
> cp src.share/* SourceMods/src.share/

where the cp path is coming from the files in this repository.

3)  Finally the code can be compiled and run.
> ./your_case_name.build
> ./your_case_name.submit

That's it! ** assuming your machine is a support system


Finally there are also two other output variables that can be placed in the history file, TBM and CONVFRAC which are the Convective Threshold and Sub-grid Convective Fraction from the Tawfik et al. (2016) paper.  This can be done by adding the following to your user namelist.  Open user_nl_cam and add:
> fincl = 'CONVFRAC:A', 'TBM:A'



  
