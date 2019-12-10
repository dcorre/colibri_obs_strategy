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


