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
................

- Define a name for creating a new folder for your results. It will be created in pyGRBz/pyGRBz/results/   
- Define the kind of GRB simulation: 'empirical' or 'theoretical' as well as the number of GRBs    
- Create the transmissions curves, and store them in pyGRBz/pyGRBz/transmissions/colibri/   

2) Set observational strategy
.............................

- Define the filters to be used   
- So far only 2 strategies can be used, with 1, 2 or 3 channels and any combination of the filters available in each channel   
- The first set of filters is used until the GRB is detected, then it switches to the second set of filters   

3) Simulate GRB light curves
............................

