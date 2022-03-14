# venetoclax-tedizolid

This is the code for running simulations for the manuscript:

-----------------------------------------------------------------------------------------------------
"An integrative systems biology approach to overcome venetoclax resistance in acute myeloid leukemia" 
-----------------------------------------------------------------------------------------------------
by Michelle Przedborski, David Sharon, Severine Cathelin, Steven Chan, and Mohammad Kohandel

The author of the code is Michelle Przedborski
----------------------------------------------

*Important details about the parameter values (e.g. order, converting to units in the manuscript, etc.) are given in the file "Make_manuscript_plots.m". 
*Simply download the codes (and the Protein.dat) file and run individual files to simulate the experimental treatment protocols.

Overview of each file:
----------------------
1. Make_manuscript_plots.m: Simulates the untreated system, ABT-199 monotherapy, Tedizolid monotherapy, and combination therapy for the nominal kinetic parameter set 
   (see Supplemental Information). This file plots the cell numbers and cell viability for each treatment condition (Figure 4) and calls makeplots_manuscript.m to 
   create the protein expression plots for each treatment condition (Figure 3). This file also contains all of the cell number and cell viability experimental data 
   (averages and standard deviations). All of the protein expression data is contained in the file "Protein.dat".
2. Make_expt_timing_plots.m: This file contains all of the cell viability measurements from the optimal drug timing experiments and creates the plots of this data, 
   including the sigmoidal fits (Figure 7).
3. Sensitivity_analysis_ABT199_Ted.m: This script runs a local sensitivity analysis around the nominal kinetic parameter set. It calls obj_ABT199_Ted_sensitivity.m
   to run the simulations, then plots the relative sensitivities for each treatment condition (Figure 5). 
4. Run_mcl1_inhibitor.m: This file simulates ABT-199 treatment in combination with Mcl-1 inhibition and plots the resulting cell viability (Figure 6).
5. Run_drug_timing_simulations.m: This file runs the combination therapy with pre-treatment simulations and plots the response, endpoint and peak cell numbers, and 
   active Caspase-3 levels (Figure 7).

Note: Files 2-5 were created from the code in File 1 and may therefore contain additional lines of code that are not necessary for creating the plots. For example, 
File 5 currently runs ABT-199/Tedizolid monotherapy simulations unnecessarily (as of 03/14/22). In the next commit of the code, these additional lines of code will be
removed.

Dependencies:
-------------
To perform numerical integration, each of the files above call the integrating functions:
1. integratingfunction_untreated_N.m: Integrating function for cell numbers in the untreated system.
2. integratingfunction_equilibrium_no_feedback.m: Integrating function for the untreated system.
3. integratingfunction_ABT199_no_feedback.m: Integrating function for ABT-199 monotherapy.
4. integratingfunction_Tedizolid_no_feedback.m: Integrating function for Tedizolid monotherapy.
5. integratingfunction_combination_ABT_Ted_no_feedback.m: Integrating function for ABT-199/Tedizolid combination therapy.
6. integratingfunction_mcl1_inhibitor.m: Integrating function for ABT-199/Mcl-1 inhibitor combination therapy.
7. EquilibriumEventFcn_no_feedback.m: Event function to stop equilibrium simulation once steady state is reached.

ODE solver:
-----------
The integration is performed using the dde15s solver (dde15s_new.m or dde15s_updated.m since they are identical). This file was written by Dr. Lawrence F. Shampine 
(Southern Methodist University, Dallas, TX) for the manuscript: Agrawal, Vikas, et al. "A dynamic mathematical model to clarify signaling circuitry underlying 
programmed cell death control in Arabidopsis disease resistance." Biotechnology progress 20.2 (2004): 426-442. The code was updated to run in MATLAB R2018b by Jacek 
Kierzenka (Mathworks) in December 2020.
