Requirements:
-------------

Install the following packages:

- https://github.com/dcorre/pyGRBaglow
- https://github.com/dcorre/pyETC
- https://github.com/dcorre/pyGRBz


Installation:
-------------

Clone the project:   
git clone https://github.com/dcorre/colibri_obs_strategy 

Install it:   
cd colibri_obs_strategy   
python setup.py develop


Information:
------------

If you want to modify the configuration of the Colibri characteristics, edit the configuration file of Colibri in the pyETC package:   
pyETC/pyETC/telescope_database/colibri.hjson  

The name of the files describing the transmission curves of the optical elements and QE of the detector you defined in that json files must be defined in pyETC/pyETC/transmissions/ and corresponding subfolders


PhotoZ estimation for an observation strategy
---------------------------------------------

You can run either the jupyter notebook or python script named Estimation_photoZ_Mock_Sample in colibri_obs_strategy/notebooks/


These files are commented, the main steps are:

1) General stuff

- Define a name for creating a new folder for your results. It will be created in pyGRBz/pyGRBz/results/   
- Define the kind of GRB simulation: 'empirical' or 'theoretical' as well as the number of GRBs    
- Create the transmissions curves, and store them in pyGRBz/pyGRBz/transmissions/colibri/   

2) Set observational strategy

- Define the filters to be used   
- So far only 2 strategies can be used, with 1, 2 or 3 channels and any combination of the filters available in each channel   
- The first set of filters is used until the GRB is detected, then it switches to the second set of filters   
- Define the time at which to start and end the observations, the exposure time for both set of filters and the dead time   

3) Simulate GRB light curves

- 2 possibilities: 
        - empirical model, simple power law both in frequencies and time   
	- using distributions from Kann et al. 2006 and 2010   
        - or self-defined distributions on the parameters    
	- Only suited to study a short period of time. Not to study efficiency at different times.   

	- theoretical model, based on the standard synchrotron model of Granot & Sari 2002   
        - using parameter distributions of D. Turpin found when fitting real GRB light curves   
	- or self-defined parameter distributins   
        - Can study the efficiency at different times with this model.   
 

- The extinction law are randomly drawn: 2/3 SMC, 1/3 MW, 1/3 LMC for moderate Av.   
- For dusty GRBs (Av>1mag), equirepartition between the 3 laws. Dusty GRBs are considered to be located at z <~4   

4) Convolve these light curves with telescope response using pyETC
------------------------------------------------------------------

- SNR considered for a detection: 3   
- Calibration uncertainty used: 0.04 and 0.06 mag for DDRAGO and CAGIRE respectively.   
- No images are produced, only the ETC is used.   


5) Run photoz MCMC
------------------

- 4 parameters to fit: z, Av, beta (spectral slope) and a scaling factor.    
- The priors range are flat and can be user defined.   
- The code is ran for each type of extinction law and without dust.   
- The MCMC is ran once for 300 steps with 30 walkers. Then the values of the best chi2 are used as initial values to start a new run for 300 steps with the 50 first steps as burnin steps (i.e. not considered in statistical analysis). This configuration can be changed, and one might want to run the MCMC only once tih more steps and a larger burnin phase in order not to get stock in a local mimima, and to visualise all plusible solutions. Especially when there is a GRB at high redshift or highly extinguished.   
- To chose the extinction law that best fits the SED, the probability of each walker at each step are summed up and we pick the extinction law with the highest probability.    


6) Outputs
----------

- Plot the zphot vs zsim, and the relative uncertainty   
- Write a short summary of the statistics in a txt file, results_summary.txt, in the results folder.   


