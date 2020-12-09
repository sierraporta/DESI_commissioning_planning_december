### Scripts to compute effective exposure times for DESI observations  

This repository contains scripts and files provided by D. Kirkby to compute the `DEPTH` for DESI observations

The computation is done in two steps.

**First Step*

Use the script `desi_exposures_gfa.py` (`python desi_exposures_gfa.py sv0`, or `python desi_exposures_gfa.py minisv2`) to read the GFA images and build summary table for each exposure that contains:

[`expid`, `night`, `tileid`, `exptime`, `mjdobs`, `tilera`, `tiledec`, `ngfa`, `ebv`].




**STEP 1**: The first script called "__desi_exposures_gfa.py__" is an original script by D.Kirkby used for SV0 and miniSV2. The main objective of this script is to take the outputs of the GFA condition observation for day and build a summary table that essentially contains minimum, maximum and median values, respectively:  
['airmass', 'moon_sep_deg', 'transparency', 'fwhm_asec', 'sky_mag_ab', 'fiber_fracflux']  
as well as information on each tile and exposition_ID:  

This script generates a table with these values and it is saved to disk.  

To run this script we do (for SV0):

> python desi_exposures_gfa.py sv0

This action generate the file "desi_sv0_exposures_gfa.fits".

**STEP 2**: The second script called "__depth_calculation.py__" take the output generate for the fist script and does the following (also based on D.Kirkby's calculation):      

1. Read the spectra of each tile and exposure and divide into bands b, r and z for each one.     
2. Read the DESIMODEL yields for each spectrograph from a local copy obtained using a console command:        
> svn co https://desi.lbl.gov/svn/code/desimodel/tags/0.13.0/data/throughput thru13
3. Now read the GFA measures of transparency and fiber fraction that have been obtained with the first script (“__desi_exposures_gfa.py__”)     
4. The average sky for a single exposure is estimated in photons/sec detected in each camera. Then two sky models are used (for these models you need two files that are copied locally "dark_desimodel.fits" and "dark_eso.fits"):     
- EXPSKY is the mean sky in electrons detected with 100A smoothing     
- FIDSKY is the fiducial sky "dark zenith" with 100A smoothing     
7. For each tile and each exposure and in each band (brz), then it is calculated:

                     DEPTH = EXPTIME x (TRANSP/1.0)^2 x (FRACFLUX/0.56)^2 x (FIDSKY/EXPSKY)

Using the TRANSP and FRACFLUX values from the results of the first script. Finally, these values are added to the first .fits generated with the first script.

It just runs doing:

> python depth_calculation.py sv0

This action generate the file "__desi_sv0_exposures_gfa_with_depth.fits__" and looks like as

![](Screenshot_file.png)
