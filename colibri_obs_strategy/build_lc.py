#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
from astropy.table import Table,Column,vstack
from astropy.io import ascii
import astropy.units as u
from scipy.stats import halfnorm, truncnorm,norm,uniform
import hjson
import gc
import os 
import subprocess as sp
import random
import imp


class obs_strategy():
   """

   Parameters
   ----------

   Returns
   -------

   """
   def __init__(self,wvl_step=100,wvl_min=2500,wvl_max=25000,plot=True,resdir='Default'):
        """
        Class Constructor.

        """
        # Find path to 'obs_strategy' and 'pyGRBz' python packages. Need to be installed.
        _, path, _ = imp.find_module('colibri_obs_strategy')
        _, path_pyGRBz, _ = imp.find_module('pyGRBz')

        self.wvl_step=wvl_step
        self.wvl_min=wvl_min
        self.wvl_max=wvl_max
        self.plot = plot
        self.path=path
        self.path_pyGRBz = path_pyGRBz
        self.resdir=resdir
        # Put results directly in pyGRBz results folder
        self.respath = self.path_pyGRBz+'/results/'+self.resdir+'/Obs_strategy/'
        if not os.path.exists(self.respath): os.makedirs(self.respath)

        # Write FWHM of the filters by hand. Need to automatised!
        #self.optics_fwhm=dict(g=0.26,r=0.26,i=0.37,z=0.37,y=0.37,J=0.59,H=0.73,gri=0.3,zy=0.37)

        return None

   def set_filters(self,vis1=[['gri'],['g','r','i']],vis2=[['zy'],['z','y']],nir=[['J','H'],['J','H']]):
       """ set the filters to use. Set None if not using one channel. """     

       self.filter_vis1=vis1
       self.filter_vis2=vis2
       self.filter_nir=nir

       # Compute number of channels
       nb_channels=0
       list_channels = []
       if self.filter_vis1 is not None:
           nb_channels+=1
           list_channels.append('vis1')
       if self.filter_vis2 is not None:
           nb_channels+=1
           list_channels.append('vis2')
       if self.filter_nir is not None:
           nb_channels+=1
           list_channels.append('nir')
       self.nb_channels = nb_channels
       self.list_channels = list_channels
       return None
   
   def set_strategy_obs(self,load_strategy=False,filename='Default_strategy.hjson',t_since_burst=60,t_start=0,t_end=300,t_dead=5,texp=[30,30],id_start=1, filter_set=0, telescope_name='colibri'):
       """ Set the observational strategy

       Parameters
       ----------
       load_strategy: boolean (default: False)
                      load an existing observationnal strategy (not implemented yet)

       t_since_burst: float (default: 60)
                      time since trigger (t0), in seconds

       t_start: float (default: 0)
                starting time of observations (t-t0), in seconds

       t_end: float (default: 300)
              end time for observations with respect to t0 (t-t0), in seconds

       t_dead: float (default: 5)
               dead time between observations, in seconds

       t_exp1: list float (default: 30)
               exposure time for both strategies, in seconds

       filter_set: integer (default: 0)
                 filter strategy to use

       Returns
       -------
       None. The observationnal strategy is stored in an astropy table to self.strategy_obs
       """

       if load_strategy:
            
           with open(filename,encoding='utf-8') as f:
               strategy_file=hjson.load(f)
          
           obs_name_table=[]
           nb_obs_table=[]
           obs_time_table=[]
           t_since_burst_table=[]
           channel_band=[]
           channel_mag=[]
           channel_magerr=[]
           texp_table=[]
           phot_sys=[]
           tel_name=[]
           detection=[]
           wvl_eff=[]
        
           #Setting the srategy:
           self.t_since_burst = strategy_file['T_since_trigger']   # Time since burst in seconds (t0), in seconds
           self.t_dead=strategy_file['Dead_time']        # dead time between observations
           self.nb_channels=strategy_file['Nb_channels']
           current_time=self.t_since_burst
            
           
           
           if self.nb_channels == 3: self.list_channels=['vis1','vis2','nir']
           elif self.nb_channels == 2: self.list_channels=['vis','nir']
           counter_channel=1
           for i in range(len(strategy_file['strategy'])):
               i+=1   # Python indexing starts at 0
               
               #obs_time_table.append(t_start)
                   
               for channel in self.list_channels:
                   Texp = strategy_file['strategy']['%d'%i]['Texp']
                   band = list(strategy_file['strategy']['%d'%i][channel].items())[0][0]
                   #print (band)
                   for nb_exp in range(strategy_file['strategy']['%d'%i][channel][band]):
                       tot_nb_exp=strategy_file['strategy']['%d'%i][channel][band]
                       nb_exp+=1   # Python indexing starts at 0
                       obs_name_table.append(None)
                       texp_table.append(int((Texp-(tot_nb_exp-1)*self.t_dead)/tot_nb_exp))
                       t_since_burst_table.append(current_time+int(((Texp-(tot_nb_exp-1)*self.t_dead)/tot_nb_exp+self.t_dead)*(nb_exp-1)))
                       channel_mag.append(-99.0)
                       channel_magerr.append(-99.0)
                       phot_sys.append('AB')
                       tel_name.append(telescope_name)
                       detection.append(0)
                       wvl_eff.append(-99)
                       nb_obs_table.append( i)
                       channel_band.append(band)
                   
               current_time+=strategy_file['strategy']['%d'%i]['Texp']+self.t_dead
                
               
                
       else:
           strategy_dict = defaultdict(lambda: defaultdict(dict))

           #Setting the srategy:
           self.t_since_burst = t_since_burst   # Time since burst in seconds (t0), in seconds
           self.t_start=t_start      # starting time for observations to t_since_burst (t-t0), in seconds
           self.t_end=t_end          # end time for observations to t0 (t-t0), in seconds
           self.t_dead=t_dead        # dead time between observations
           self.texp = texp[filter_set]

           # initial values of indexes to loop over the filters.
           #0 means that we start with the first band given by the user
           channel1=0
           channel2=0
           channel3=0

           starting_time_obs=[]
           i=t_start
           while i < t_end:
               starting_time_obs.append(i)
               i += texp[filter_set] + t_dead
           starting_time_obs=np.array(starting_time_obs) 


           obs_name='obs_'
           obs_number=id_start-1

           obs_name_table=[]
           nb_obs_table=[]
           obs_time_table=[]
           t_since_burst_table=[]
           channel_band=[]
           channel_mag=[]
           channel_magerr=[]
           texp_table=[]
           phot_sys=[]
           tel_name=[]
           detection=[]
           wvl_eff=[]

           for t_start in starting_time_obs:
               obs_number+=1
               obs_ref=obs_name+str(obs_number)
               if self.nb_channels == 3:
                   strategy_dict['Observation_sets'][obs_number]={'t_start': t_start, 'texp':texp[filter_set], 'filter':{self.filter_vis1[filter_set][channel1],self.filter_vis2[filter_set][channel2],self.filter_nir[filter_set][channel3]},'obs_mag':{self.filter_vis1[filter_set][channel1]: None,self.filter_vis2[filter_set][channel2]: None,self.filter_nir[filter_set][channel3]: None}, 'err_mag':{self.filter_vis1[filter_set][channel1]:None,self.filter_vis2[filter_set][channel2]:None,self.filter_nir[filter_set][channel3]:None}}


                   for i in range(self.nb_channels):
                       obs_name_table.append(None)
                       obs_time_table.append(t_start)
                       texp_table.append(texp[filter_set])
                       t_since_burst_table.append(t_since_burst+t_start)
                       channel_mag.append(-99.0)
                       channel_magerr.append(-99.0)
                       phot_sys.append('AB')
                       tel_name.append(telescope_name)
                       detection.append(0)
                       wvl_eff.append(-99)
                       nb_obs_table.append( obs_number )

                   channel_band.append(self.filter_vis1[filter_set][channel1])
                   channel_band.append(self.filter_vis2[filter_set][channel2])
                   channel_band.append(self.filter_nir[filter_set][channel3])

                   channel1 = (channel1+1) % len(self.filter_vis1[filter_set])
                   channel2 = (channel2+1) % len(self.filter_vis2[filter_set])
                   channel3 = (channel3+1) % len(self.filter_nir[filter_set])


               elif self.nb_channels == 2:
                   if self.filter_vis1 and self.filter_nir:
                       strategy_dict['Observation_sets'][obs_number]={'t_start': t_start, 'texp':texp[filter_set], 'filter':{self.filter_vis1[filter_set][channel1],self.filter_nir[filter_set][channel3]},'obs_mag':{self.filter_vis1[filter_set][channel1]: None,self.filter_nir[filter_set][channel3]: None}, 'err_mag':{self.filter_vis1[filter_set][channel1]:None,self.filter_nir[filter_set][channel3]:None}}
                   elif self.filter_vis1 and self.filter_vis2:
                       strategy_dict['Observation_sets'][obs_number]={'t_start': t_start, 'texp':texp[filter_set], 'filter':{self.filter_vis1[filter_set][channel1],self.filter_vis2[filter_set][channel2]},'obs_mag':{self.filter_vis1[filter_set][channel1]: None,self.filter_vis2[filter_set][channel2]: None}, 'err_mag':{self.filter_vis1[filter_set][channel1]:None,self.filter_vis2[filter_set][channel2]:None}}
                   elif self.filter_vis2 and self.filter_nir:
                       strategy_dict['Observation_sets'][obs_number]={'t_start': t_start, 'texp':texp[filter_set], 'filter':{self.filter_vis2[filter_set][channel2],self.filter_nir[filter_set][channel3]},'obs_mag':{self.filter_vis2[filter_set][channel2]: None,self.filter_nir[filter_set][channel3]: None}, 'err_mag':{self.filter_vis2[filter_set][channel2]:None,self.filter_nir[filter_set][channel3]:None}}

                   if self.filter_vis1:
                       channel1 = (channel1+1) % len(self.filter_vis1[filter_set])
                       channel_band.append(self.filter_vis1[filter_set][channel1])

                   if self.filter_vis2:
                       channel2 = (channel2+1) % len(self.filter_vis2[filter_set])

                       channel_band.append(self.filter_vis2[filter_set][channel2])
                   if self.filter_nir:
                       channel3 = (channel3+1) % len(self.filter_nir[filter_set])
                       channel_band.append(self.filter_nir[filter_set][channel3])


                   for i in range(self.nb_channels):
                       obs_name_table.append(None)
                       obs_time_table.append(t_start)
                       texp_table.append(texp[filter_set])
                       t_since_burst_table.append(t_since_burst+t_start)
                       channel_mag.append(-99.0)
                       channel_magerr.append(-99.0)
                       phot_sys.append('AB')
                       tel_name.append(telescope_name)
                       detection.append(0)
                       wvl_eff.append(-99)
                       nb_obs_table.append( obs_number )



                   

           #print (nb_obs_table)
           strategy_dict['time_since_burst']=t_since_burst
           strategy_dict['nobs']=len(starting_time_obs)
           #print (OrderedDict(sorted(obs_dict.items(), key=lambda t: t[0])))
           #print (strategy_dict['Observation_sets'][1]['err_mag'])
          
            
       #create astropy table
       strategy_obs = Table([obs_name_table,nb_obs_table,t_since_burst_table,texp_table,channel_band,tel_name,channel_mag,channel_magerr,detection,phot_sys],names=['Name','Nb_obs','time_since_burst','Exp_time','band','telescope','mag','mag_err','detection','phot_sys'])
       #strategy_obs['Obs_time'].unit = u.second
       strategy_obs['Exp_time'].unit = u.second
       strategy_obs['time_since_burst'].unit = u.second
       #strategy_table['mag'].unit = 'AB mag'
       #strategy_table['mag_err'].unit = 'AB mag'
       strategy_obs['mag'].format='3.3f'
       strategy_obs['mag_err'].format='3.3f'

       #print (strategy_obs)
       self.strategy_obs=strategy_obs
       ascii.write(self.strategy_obs,self.respath+'Observational_strategy.txt',overwrite=True)
       ascii.write(self.strategy_obs,self.respath+'Observational_strategy.tex',format='latex',overwrite=True)

       #with open(respath+'/Observational_strategy.hjson', 'w') as fp:
       #     hjson.dump(strategy_dict, fp)

       return strategy_obs

   def set_grb_params(self,params=None,random_flag=False,load_params=False,params_file='grb_params_emp.tex',load_min=1,load_max=10,load_distrib=True,distrib_file='Turpin_GS02.txt',num_samples=1,model='gs02',plot_distribution=False,seed=1234567890,name_prefix_start=0):
       """ set the parameters caracterising the grb afterglow

       Parameters
       ----------
       model: string (default: 'gs02')
              name of the model to use. 'gs02': Granot&Sari2002 / 'kann': empirical data from Kann2010 

       params: astropy table (default: None)
               contains the parameters to model the GRB either with Granot&Sari or with empirical model.
               either list or ['type',mean,std,a,b] a and b being the limits of allowed range. 
               'type' can be: 'constant','norm','halfnorm','truncnorm','linspace','uniform'
               See function gen_distribution for more details

       random_flag: 	boolean (default: False)
               	True: draw parameters from distributions. False: Parameters are read from files or given

       seed: 	integer (default: 1234567890)
             	seed for random number generator

       num_samples: integer (default: 1)
                    number of GRB to generate

       load_params: boolean (default: False)
                    True: load parameters from file

       params_file: string (default: 'grb_params_emp.tex')
                    Path and name of the file containing the GRB parameters to be loaded

       load_min: integer (default: 1)
                 first line to consider in parameter file. For instance, if 1 then it starts to read from the first line

       load_max: integer (default: 10)
                 last line to consider in parameter file. For instance, if 10 then it stops to read at the line 10. -1 means until the end

       load_distrib: boolean (default: True)
                     True: use parameters distribution to generate random values. If one parameter is not in this file then the values in "params" are used

       distrib_file: string (default: 'Turpin_GS02.txt')
                     name of the file containing the parameters distribution

       plot_distribution: boolean (default: False)
                          True: Plot the parameters distributions used for the current run

       name_prefix_start: integer (default: 1)
                          prefix for GRB names (do not change it)

       Returns
       --------
       None. The grb parameters are stored in an astropy table in self.grb_params

       """
       #params: either list or ['type',mean,std,a,b] a and b being the limits of allowed range

       # GRB parameters
       #----------------
       self.model=model
       self.num_samples=num_samples
       self.seed=seed
       self.model=model

       if load_params==True:
           params=ascii.read(params_file)
           params=params[load_min-1:load_max]
           print (params)

           name=params['name']

           Av_samples = params['Av_host'] 
           z_samples=params['z'] 
           extlaws_samples=params['ExtLaw']
           #data_old=ascii.read('distrib/Granot_Sari/samples300_gs02.tex')
           #Av_samples=data_old['Av_host'][load_min-1:load_max]
           #z_samples=data_old['z'][load_min-1:load_max]
           #extlaws_samples=data_old['ExtLaw'][load_min-1:load_max]

           if model=='kann':
               Rc=params['Rc']
               alpha=params['alpha']
               beta=params['beta']
               wvl0=params['wvl0']
               t0=params['t0']
               F0_z=params['F0']

           elif model=='gs02':
               n0=params['n0']
               eps_b=params['eps_b']
               eps_e=params['eps_e']
               E_iso=params['E_iso']
               p=params['p']
               Y=params['Y']
               ism_type=params['ism_type']
               Mdot_loss=params['Mdot_loss']
               Vw=params['Vw']
               eta=params['eta']

       elif load_params==False:
            
           if load_distrib==True: 
               distrib=ascii.read(self.path+'/distrib/'+distrib_file,delimiter='\t',fill_values=[('-',np.nan)])
               from scipy.stats import rv_histogram
               #require scipy >=19.0
               #print (distrib)
           if model=='kann':
               if random_flag == False:
                   Av_samples = params['Av_host']
                   extlaws_samples=random.choices(params['ExtLaw'],k=num_samples)
                   z_samples=params['z']
          
                   # mag Rc
                   Rc=params['Rc']

                   # time at which Rc is extracted
                   t0=params['t0']

                   # wavelength at which Rc is extracted
                   wvl0=params['wvl0']

                   # temporal index
                   alpha=params['alpha']
    
                   # spectral index
                   beta=params['beta']
 

               elif random_flag == True:
                   if self.seed!=-1: np.random.seed(self.seed)
                   #same_z_AV_samples=ascii.read(self.path+'/distrib/Granot_Sari/samples300_gs02.tex')
                   
                   # floats between 0 and 1
                   rand_list=np.random.rand(num_samples)

                   #z
                   if load_distrib and 'z' in distrib.colnames:
                       hist = np.histogram(distrib['z'], bins=20, density=True,range=(np.min(distrib['z']),np.max(distrib['z'])))
                       hist_dist = rv_histogram(hist)
                       z_samples=hist_dist.ppf(rand_list)
                   else:
                       z_samples = self.gen_distribution(params['z'])
                   #z_samples=same_z_AV_samples['z']

                   #Av
                   if load_distrib and 'Av' in distrib.colnames:
                       hist = np.histogram(distrib['Av'], bins=20, density=True,range=(np.min(distrib['Av']),np.max(distrib['Av'])))
                       hist_dist = rv_histogram(hist)
                       Av_samples=hist_dist.ppf(rand_list)
                   else:
                       Av_samples = self.gen_distribution(params['Av_host'])
                   #Av_samples=same_z_AV_samples['Av_host']

                   #Extinction law
                   if load_distrib and 'ExtLaw' in distrib.colnames:
                       extlaws_samples=random.choices(distrib['ExtLaw'],k=num_samples)
                   else:
                       extlaws_samples=random.choices(params['ExtLaw'],k=num_samples)

                   #Rc
                   if load_distrib and 'Rc' in distrib.colnames:
                       hist = np.histogram(distrib['Rc'], bins=13,density=True,range=(np.min(distrib['Rc']),np.max(distrib['Rc'])))
                       hist_dist = rv_histogram(hist)
                       Rc=hist_dist.ppf(rand_list)
                   else:
                       Rc=self.gen_distribution(params['Rc'])

                   #beta
                   if load_distrib and 'beta' in distrib.colnames:
                       hist = np.histogram(distrib['beta'], bins=13,density=True,range=(np.min(distrib['beta']),np.max(distrib['beta'])))
                       hist_dist = rv_histogram(hist)
                       beta=hist_dist.ppf(rand_list)
                   else:
                       beta=self.gen_distribution(params['beta'])

                   #alpha
                   if load_distrib and 'alpha' in distrib.colnames:
                       hist = np.histogram(distrib['alpha'], bins=13,density=True,range=(np.min(distrib['alpha']),np.max(distrib['alpha'])))
                       hist_dist = rv_histogram(hist)
                       alpha=hist_dist.ppf(rand_list)
                   else:
                       alpha=self.gen_distribution(params['alpha'])

                   #wvl0
                   if load_distrib and 'wvl0' in distrib.colnames:
                       hist = np.histogram(distrib['wvl0'], bins=13,density=True,range=(np.min(distrib['wvl0']),np.max(distrib['wvl0'])))
                       hist_dist = rv_histogram(hist)
                       wvl0=hist_dist.ppf(rand_list)
                   else:
                       wvl0=self.gen_distribution(params['wvl0'])

                   #t0
                   if load_distrib and 't0' in distrib.colnames:
                       hist = np.histogram(distrib['t0'], bins=13,density=True,range=(np.min(distrib['t0']),np.max(distrib['t0'])))
                       hist_dist = rv_histogram(hist)
                       t0=hist_dist.ppf(rand_list)
                   else:
                       t0=self.gen_distribution(params['t0'])

               #wavelength of reference in angstroms
               #wvl0=6400
               #time of reference in seconds 
               #t0=86.4*np.ones(num_samples)


               #correct for redshifting but not for extinction and IGM
               #it will be done later, here we want to calculate the Flux normalisation
               from astropy.cosmology import Planck15 as cosmo
               # Calculate the luminosity distance
               #cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
               z0=1
               F0_z=[]
               for mag, z,be,al in list(zip(Rc,z_samples,beta,alpha)):
                   F0_z0=3080*10**(-0.4*mag)   #Jy
                   F_z1_nu_t=F0_z0*(1+z)/(1+z0)*(cosmo.luminosity_distance(z0).value/cosmo.luminosity_distance(z).value)**2*((1+z0)/(1+z))**(-al)*((1+z)/(1+z0))**(-be)
                   #print(z,(1+z)/(1+z0)*(cosmo.luminosity_distance(z0).value/cosmo.luminosity_distance(z).value)**2*((1+z0)/(1+z))**(-al)*((1+z)/(1+z0))**(-be))
                   F0_z.append(F_z1_nu_t)
                   #print (mag,F01,z,F0_z[-1],F0_z[-1]/(1+z))

               Rc=np.array(["{0:.2e}".format(x) for x in Rc],dtype=float)
               alpha=np.array(["{0:.2f}".format(x) for x in alpha],dtype=float)
               beta=np.array(["{0:.2f}".format(x) for x in beta],dtype=float)
               F0_z=np.array(["{0:.2e}".format(x) for x in F0_z],dtype=float)
               self.plot_distrib(F0_z,50,plot_distribution,'F0 (Jy)',log=True)
               self.plot_distrib(Rc,20,plot_distribution,'Rc mag')
               self.plot_distrib(alpha,20,plot_distribution,'alpha')
               self.plot_distrib(beta,20,plot_distribution,'beta')

           elif model=='gs02':  #Granot and Sari
               if random_flag == False:
                   Av_samples = params['Av_host'] 
                   z_samples=params['z'] 
                   extlaws_samples=random.choices(params['ExtLaw'],k=num_samples)
                    
                   #Density of the environment
                   n0= params['n0']
                   # Fraction of energy to magneic field
                   eps_b=params['eps_b']

                   # Fraction of energy to electrons
                   eps_e=params['eps_e']
                   # E_iso
                   E_iso=params['E_iso']

                   # electron energy distribution power index  
                   p=params['p']
                   # Inverse compton 
                   Y=params['Y']
                   #ism_type
                   ism_type=params['ism_type']
                   #Mdot_loss
                   Mdot_loss=params['Mdot_loss']
                   #Vw
                   Vw=params['Vw']
                   # eta
                   eta=params['eta']
                   #ism_type
                   ism_type=params['ism_type']

               elif random_flag == True:
                   # Fixed the seed, or not...
                   if self.seed!=-1: np.random.seed(self.seed)

                   # floats between 0 and 1
                   rand_list=np.random.rand(num_samples)
                   
                   #z
                   if load_distrib and 'z' in distrib.colnames:
                       hist = np.histogram(distrib['z'], bins=20,density=True,range=(np.min(distrib['z']),np.max(distrib['z'])))
                       hist_dist = rv_histogram(hist)
                       z_samples=hist_dist.ppf(rand_list)
                   else:
                       z_samples = self.gen_distribution(params['z'])

                   #z
                   if load_distrib and 'Av' in distrib.colnames:
                       hist = np.histogram(distrib['Av'], bins=20,density=True,range=(np.min(distrib['Av']),np.max(distrib['Av'])))
                       hist_dist = rv_histogram(hist)
                       Av_samples=hist_dist.ppf(rand_list)
                   else:
                       Av_samples = self.gen_distribution(params['Av_host'])
                   
                   #Extinction law
                   if load_distrib and 'ExtLaw' in distrib.colnames:
                       extlaws_samples=random.choices(distrib['ExtLaw'],k=num_samples)
                   else:
                       extlaws_samples=random.choices(params['ExtLaw'],k=num_samples)

                   #n0
                   if load_distrib == True and 'n0' in distrib.colnames:
                       hist = np.histogram(np.log10(distrib['n0']), bins=20,density=True,range=(np.log10(np.min(distrib['n0'])),np.log10(np.max(distrib['n0']))))
                       hist_dist = rv_histogram(hist)
                       n0=hist_dist.ppf(rand_list)
                       #print (n0)
                       n0=10**n0
                       #print (n0)
                   else:
                       n0=self.gen_distribution(params['n0'])

                   #eps_b
                   if load_distrib and 'eps_b' in distrib.colnames:
                       #plt.figure()
                       #plt.hist(np.log10(distrib['eps_b']))
                       #plt.show()
                       hist = np.histogram(np.log10(distrib['eps_b']), bins=20,density=True,range=(np.log10(np.min(distrib['eps_b'])),np.log10(np.max(distrib['eps_b']))))
                       hist_dist = rv_histogram(hist)
                       eps_b=hist_dist.ppf(rand_list)
                       eps_b=10**eps_b
                       #plt.figure()
                       #plt.hist(eps_b)
                       #plt.show()
                   else:
                       eps_b=self.gen_distribution(params['eps_b'])

                   #eps_e
                   if load_distrib and 'eps_e' in distrib.colnames:
                       hist = np.histogram(np.log10(distrib['eps_e']), bins=20,density=True,range=(np.log10(np.min(distrib['eps_e'])),np.log10(np.max(distrib['eps_e']))))
                       hist_dist = rv_histogram(hist)
                       eps_e=hist_dist.ppf(rand_list)
                       eps_e=10**eps_e
                   else:
                       eps_e=self.gen_distribution(params['eps_e'])

                   #p
                   if load_distrib and 'p' in distrib.colnames:
                       hist = np.histogram(distrib['p'], bins=20,density=True,range=(np.min(distrib['p']),np.max(distrib['p'])))
                       hist_dist = rv_histogram(hist)
                       p=hist_dist.ppf(rand_list)
                   else:
                       p=self.gen_distribution(params['p'])

                   #eta
                   if load_distrib and 'eta' in distrib.colnames:
                       hist = np.histogram(distrib['eta'], bins=20,density=True,range=(np.min(distrib['eta']),np.max(distrib['eta'])))
                       hist_dist = rv_histogram(hist)
                       eta=hist_dist.ppf(rand_list)
                   else:
                       eta=self.gen_distribution(params['eta'])
          
                   #E_iso
                   if load_distrib and 'E_iso' in distrib.colnames:
                       hist = np.histogram(np.log10(distrib['E_iso']), bins=20,density=True,range=(np.log10(np.min(distrib['E_iso'])),np.log10(np.max(distrib['E_iso']))))
                       hist_dist = rv_histogram(hist)
                       E_iso=10**hist_dist.ppf(rand_list)
                   else:
                       E_iso=10**self.gen_distribution(params['E_iso'])

                   #E_K
                   if load_distrib and 'E_K' in distrib.colnames:
                       hist = np.histogram(np.log10(distrib['E_K']), bins=20,density=True,range=(np.log10(np.min(distrib['E_K'])),np.log10(np.max(distrib['E_K']))))
                       hist_dist = rv_histogram(hist)
                       E_K=10**hist_dist.ppf(rand_list)
                   else:
                       pass
                       # Need to implement it properly
                       #E_K=10**self.gen_distribution(params['E_K'])

                   #Y
                   if load_distrib and 'Y' in distrib.colnames:
                       hist = np.histogram(distrib['Y'], bins=20,density=True,range=(np.min(distrib['Y']),np.max(distrib['Y'])))
                       hist_dist = rv_histogram(hist)
                       Y=hist_dist.ppf(rand_list)
                   else:
                       Y=self.gen_distribution(params['Y'])

                   #Mdot_loss
                   if load_distrib and 'Mdot_loss' in distrib.colnames:
                       hist = np.histogram(distrib['Mdot_loss'], bins=20,density=True,range=(np.min(distrib['Mdot_loss']),np.max(distrib['Mdot_loss'])))
                       hist_dist = rv_histogram(hist)
                       Mdot_loss=hist_dist.ppf(rand_list)
                   else:
                       Mdot_loss=self.gen_distribution(params['Mdot_loss'])
 
                   #Vw
                   if load_distrib and 'Vw' in distrib.colnames:
                       hist = np.histogram(distrib['Vw'], bins=20,density=True,range=(np.min(distrib['Vw']),np.max(distrib['Vw'])))
                       hist_dist = rv_histogram(hist)
                       Vw=hist_dist.ppf(rand_list)
                   else:
                       Vw=self.gen_distribution(params['Mdot_loss'])

                   #ism_type
                   if load_distrib and 'ism_type' in distrib.colnames:
                       hist = np.histogram(distrib['ism_type'], bins=20,density=True,range=(np.min(distrib['ism_type']),np.max(distrib['ism_type'])))
                       hist_dist = rv_histogram(hist)
                       ism_type=hist_dist.ppf(rand_list)
                   else:
                       ism_type=self.gen_distribution(params['ism_type'])
               #print (eps_b)
               n0=np.array(["{0:.2e}".format(x) for x in n0],dtype=float)
               eps_e=np.array(["{0:.2e}".format(x) for x in eps_e],dtype=float)
               eps_b=np.array([float("{0:.2e}".format(x)) for x in eps_b],dtype=float)
               E_iso=np.array(["{0:.2e}".format(x) for x in E_iso],dtype=float)
               p=np.array(["{0:.2f}".format(x) for x in p],dtype=float)
               Y=np.array(["{0:.2f}".format(x) for x in Y],dtype=float)
               Mdot_loss=np.array(["{0:.2f}".format(x) for x in Mdot_loss],dtype=float)
               Vw=np.array(["{0:.2f}".format(x) for x in Vw],dtype=float)
               eta=np.array(["{0:.2f}".format(x) for x in eta],dtype=float)
               ism_type=np.array(["{0:d}".format(int(x)) for x in ism_type],dtype=int)
            
 
           #print (np.log10(n0))
           Av_samples=np.array(["{0:.2f}".format(x) for x in Av_samples],dtype=float)
           z_samples=np.array(["{0:.2f}".format(x) for x in z_samples],dtype=float)



           #Add a reference name for each grb
           grb_ref='GRB_'
           name=[]
           for i in range(num_samples):
               i+=name_prefix_start#+1 # starts at 0
               name.append(grb_ref+str(i))
       if model=='kann':
           grb_params=Table([name,z_samples,Av_samples,extlaws_samples,F0_z,wvl0,t0,alpha,beta,Rc],names=['name','z','Av_host','ExtLaw','F0','wvl0','t0','alpha','beta','Rc'])
           self.plot_distrib(Av_samples,20,plot_distribution,'Av_host')
           self.plot_distrib(z_samples,20,plot_distribution,'zsim')
           self.plot_distrib(F0_z,50,plot_distribution,'F0 (Jy)',log=True)
           self.plot_distrib(Rc,20,plot_distribution,'Rc mag')
           self.plot_distrib(alpha,20,plot_distribution,'alpha')
           self.plot_distrib(beta,20,plot_distribution,'beta')
       elif model=='gs02':
           self.plot_distrib(Av_samples,20,plot_distribution,'Av_host')
           self.plot_distrib(z_samples,20,plot_distribution,'zsim')
           self.plot_distrib(np.log10(n0),20,plot_distribution,'n0',log=False)
           self.plot_distrib(np.log10(eps_b),20,plot_distribution,'eps_b',log=False)
           self.plot_distrib(np.log10(eps_e),20,plot_distribution,'eps_e',log=False)
           self.plot_distrib(np.log10(E_iso),20,plot_distribution,'E_iso')
           #self.plot_distrib(np.log10(E_K),20,plot_distribution,'E_K_turpin')
           self.plot_distrib(np.log10(eta),20,plot_distribution,'eta')
           self.plot_distrib(np.log10(E_iso*(1-eta)/eta),20,plot_distribution,'E_K')
           self.plot_distrib(p,20,plot_distribution,'p')

           # Rather take E_K from Turpin because E_iso and eta degenerate
           #E_iso=E_K
           #eta=0.5*np.ones(num_samples)
           grb_params=Table([name,z_samples,Av_samples,extlaws_samples,E_iso,eta,eps_b,eps_e,p,Y,ism_type,n0,Mdot_loss,Vw],names=['name','z','Av_host','ExtLaw','E_iso','eta','eps_b','eps_e','p','Y','ism_type','n0','Mdot_loss','Vw'])

       #print (grb_params)
       self.grb_params=grb_params
       ascii.write(self.grb_params,self.respath+'Params.tex',format='latex',overwrite=True)
       ascii.write(self.grb_params,self.respath+'Params.txt',overwrite=True)
       if load_params == False: ascii.write(self.grb_params,'test.tex',format='latex',overwrite=True)

       return None

   def gen_distribution(self,params):
       """ Computes a distribution given the mean and variance in imputs

       Parameters
       ----------
       params: dict 
               ex: dict(z=['type',mean,std,a,b],beta=['type',mean,std,a,b])

               type: 'constant','norm','halfnorm','truncnorm','uniform','linspace'
               constant: mean: the constant value/ std,a,b useless
               norm: mean: mean / std: variance / a,b useless
               halfnorm: mean: mean / std: variance / a,b useless
               truncnorm: mean: mean / std: variance / a,b: lower and upper bounds for allowed interval
               uniform: mean: lower limit, std: upper bounds  of interval
               linspace: mean: mean: lower limit, std: upper bounds  of interval 
       
               the number of samples is take from self.num_samples

       Returns
       -------
       x: array
           array containing the sample values
       """
       #np.random.seed(123)
           
       num_samples = int(self.num_samples)
       distype = str(params[0])
       mean = float(params[1])
       if distype != 'constant':
           std = float(params[2])
           if distype=='truncnorm':
               a=float(params[3])
               b=float(params[4])

       #if self.seed!=-1: np.random.seed(self.seed)

       # Computes the distribution
       if distype == 'constant':
           x = mean*np.ones(num_samples)
       elif distype == 'norm':
           #x = mean + sigma * np.random.randn(num_samples)
           x = norm.rvs(loc=mean,scale=std,size=num_samples)#,random_state=None)
       elif distype == 'halfnorm':
           x = halfnorm.rvs(loc=mean,scale=std,size=num_samples)#,random_state=None)
       elif distype == 'truncnorm':
           a, b = (a - mean) / std, (b - mean) / std
           x = truncnorm.rvs(a,b,loc=mean,scale=std,size=num_samples)
       elif distype=='uniform':
           x = uniform.rvs(loc=mean,scale=std,size=num_samples)
       elif distype=='linspace':
           x =np.linspace(mean,std,num_samples)
       return x


   def plot_distrib(self,distrib,num_bins,disp,label,log=False):
       """ Plot the given distribution in 'nu_old=ascii.read('../../obs_strategy/obs_strategy/distrib/Granot_Sari/samples300_gs02.tex')
m_bins' bins
    
       Parameters
       ----------
       distrib: array
             conatins the values of each sample of the distribution

       num_bins: integer
              number of bins

       disp: boolean
          True: displays the plot / False: just save it 
 
       label: string
            Av, z, beta, NHI

       Returns
       -------
       Nothing except the plot saved in file
       """
       plt.figure()
       plotpath= self.respath+'distrib/'
       if not os.path.exists(plotpath): os.makedirs(plotpath)
       if log == True: 
           distrib=np.log10(distrib)
           label = '$log_{10}$ ' + label
       n, bins, patches = plt.hist(distrib,bins=int(num_bins),density=False, facecolor='c')#,alpha=0.5)
       plt.ylabel('N',size=18)
       if label == 'Av_host': newlabel = '$A_{V host}$'
       elif label == 'zsim': newlabel = '$z_{sim}$'
       elif label == 'eps_b': newlabel = '$log_{10}$ $(\epsilon_{B})$'
       elif label == 'eps_e': newlabel = '$log_{10}$ $(\epsilon_{e})$'
       elif label == 'E_iso': newlabel = '$log_{10}$ $(E_{iso})$'
       elif label == 'E_K': newlabel = '$log_{10}$ $(E)$'
       elif label == 'eta': newlabel = '$log_{10}$ $(\eta)$'
       elif label == 'n0': newlabel = '$log_{10}$ $(n_0)$'
       else: newlabel = label
       plt.xlabel(r'%s' % newlabel,size=18)
       plt.tick_params(labelsize=15)
       plt.tight_layout()
       plt.savefig(plotpath+'distrib_'+label+'.png')
       if disp: plt.show()

       return None

   #def set_lightcurves(self,SNR4detection=3,err_mag_sys=[0.04,0.06],nb_strategy=1,configFileVIS=os.getenv('GFTSIM_DIR')+'/obs_strategy/obs_strategy/data/input_vis',configFileNIR=os.getenv('GFTSIM_DIR')+'/obs_strategy/obs_strategy/data/input_nir',resdir=os.getenv('GFTSIM_DIR')+'/obs_strategy/obs_strategy/results'):
   def set_lightcurves(self,SNR4detection=3,err_mag_sys=[0.04,0.06],configFile='Test_mock.hjson',telescope='colibri',etc_plot=False,disp=False,verbose=False,resdir='/data/lc/', fname='',RA_J2000='0h54m50.6s',DEC_J2000='+14d05m04.5s',galactic_dust_corrected=1, t_since_burst=60,t_start=0,t_end=300,t_dead=5,texp=[30,30]):
       """ Computes the observed magnitudes and substitutes it inside the observationnal strategy table

       Parameters
       ----------
       SNR4detection: float (default: 3)
                      SNR for which we consider to have a  detection
 
       err_mag_sys: list (default: [0.04,0.06])
                    systematic error to be added in quadrature with the statistical one.
                    First element is for the visible and second for NIR channel. (values from GROND)

       Returns
       -------
       None. The updated self.observations astropy table is written in .dat and .tex files  
       """
       from pyETC.pyETC import etc

       #Load two etc, one with visible detectors and the other one with cagire
       #a_vis=etc(configFile=self.path+'obs_strategy/data/input_vis')
       #a_nir=etc(configFile=self.path+'obs_strategy/data/input_nir')
       #print (sys.getsizeof(a_vis))
       obs_dict=defaultdict(lambda: defaultdict(dict))


       nb_sample=0
       p1 = len(self.grb_params)

       etc_GFT=etc(configFile=configFile,name_telescope=telescope)
       etc_GFT.information['etc_plot']=etc_plot
       etc_GFT.information['disp']=disp
       etc_GFT.information['verbose']=verbose
       etc_GFT.information['object_type']='grb_sim'

       print ("GRB Time Band")
       for p in range(p1):

           #Initailise detection to  False
           GRB_detection_status='not'
           filter_set=0
           # Need to make it in a more practical way
           
           #Reinitialize filters for first strategy to increase detection
           self.set_strategy_obs(load_strategy=False, t_since_burst=t_since_burst, t_start=t_start, t_end=t_end, t_dead=t_dead, texp=texp,filter_set=filter_set, telescope_name=telescope)
           #print (self.strategy_obs)
           #print ('GRB params:')
           #print (self.grb_params[p])
           obs_table=self.strategy_obs.copy()
           #obs_table['Name']=self.grb_params['name'][p]
           for i in range(len(obs_table)):
               obs_table['Name'][i] = self.grb_params['name'][p]

           n=len(obs_table)
           #print (n)
           #print (obs_table)

           #GRB properties
           etc_GFT.information['grb_redshift']=self.grb_params['z'][p]
           etc_GFT.information['Av_Host']=self.grb_params['Av_host'][p]
           etc_GFT.information['host_extinction_law']=self.grb_params['ExtLaw'][p]
           if self.model == 'gs02':
               etc_GFT.information['grb_model']='gs02'
               etc_GFT.information['n0']=self.grb_params['n0'][p]
               etc_GFT.information['eps_b']=self.grb_params['eps_b'][p]
               etc_GFT.information['eps_e']=self.grb_params['eps_e'][p]
               etc_GFT.information['E_iso']=self.grb_params['E_iso'][p]
               etc_GFT.information['p']=self.grb_params['p'][p]
               etc_GFT.information['eta']=self.grb_params['eta'][p]
               etc_GFT.information['Y']=self.grb_params['Y'][p]
               etc_GFT.information['ism_type']=self.grb_params['ism_type'][p]
               etc_GFT.information['Mdot_loss']=self.grb_params['Mdot_loss'][p]
               etc_GFT.information['Vw']=self.grb_params['Vw'][p]
           elif self.model == 'kann':
               etc_GFT.information['grb_model']='SPL'
               etc_GFT.information['F0']=self.grb_params['F0'][p]
               etc_GFT.information['wvl0']=self.grb_params['wvl0'][p]
               etc_GFT.information['t0']=self.grb_params['t0'][p]
               etc_GFT.information['alpha']=self.grb_params['alpha'][p]
               etc_GFT.information['beta']=self.grb_params['beta'][p]
                   
           for j in range(n):
               #Change of strategy if GRB detected 
               if GRB_detection_status == 'detected':
                   print ('GRB detected')
                   filter_set=1
                   #print (obs_table['Name','Obs_time','time_since_burst','band'])
                   #print (obs_table['time_since_burst'][j+1],obs_table['Obs_time'][j+1])
                   counter=0
                   #print (obs_table['Nb_obs'][j-1])
                   #if not the last observation chang the strategy
                   if j < n-1:
                       while obs_table['Nb_obs'][j-1] == obs_table['Nb_obs'][j+counter]:
                           #print (obs_table['Nb_obs'][j+counter])
                           counter+=1
                       #print (counter)
                       #print (obs_table['time_since_burst'][j+counter],obs_table['Obs_time'][j+counter])
                       self.set_strategy_obs(load_strategy=False,t_since_burst=obs_table['time_since_burst'][j-1],t_start=obs_table['time_since_burst'][j+counter]-obs_table['time_since_burst'][0],t_end=t_end,t_dead=t_dead,texp=texp,filter_set=filter_set,telescope_name=telescope)
                       #print ('test')
                       #print (self.strategy_obs['Name','Obs_time','time_since_burst','band'])
                       obs_table[j+counter:]=self.strategy_obs.copy()
                       for i in range(len(obs_table[j+counter:])):
                           obs_table['Name'][j+counter+i] = self.grb_params['name'][p]
                       #obs_table['Name']=self.grb_params['name'][p]
                       GRB_detection_status='already detected'
                       #print (obs_table['Name','Obs_time','time_since_burst','band'])
               nb_sample+=1

               etc_GFT.information['etc_type']='snr'
               etc_GFT.information['Nexp']=1

               if obs_table['band'][j] in self.filter_vis1[filter_set]:
                   etc_GFT.information['channel']='DDRAGO-B'
                   err_sys=err_mag_sys[0]
               elif obs_table['band'][j] in self.filter_vis2[filter_set]:
                   etc_GFT.information['channel']='DDRAGO-R'
                   err_sys=err_mag_sys[0]
               elif obs_table['band'][j] in self.filter_nir[filter_set]:
                   etc_GFT.information['channel']='CAGIRE'
                   err_sys=err_mag_sys[1]

               etc_GFT.information['exptime']=obs_table['Exp_time'][j]
               etc_GFT.information['filter_band']=obs_table['band'][j]
               #etc_GFT.information['Fwhm_psf_opt']=self.optics_fwhm[obs_table['band'][j]]

               etc_GFT.information['t_sinceBurst']=obs_table['time_since_burst'][j]/86400


               #print (obs_table['time_since_burst'][j])
               etc_GFT.sim()
               mag=etc_GFT.information['mag']
               SNR=etc_GFT.information['SNR']
               photometry_system=etc_GFT.information['photometry_system']
               #del etc
               gc.collect()


               print (self.grb_params['name'][p],obs_table['time_since_burst'][j],obs_table['band'][j])
               if SNR >= SNR4detection:
                   err_mag = np.sqrt((float(2.5*np.log10(1.+1./SNR)))**2 + (err_sys)**2)
                   #Draw the magitude from a gaussian distribution
                   mag=np.random.normal(mag, err_mag, 1)
                   #mag=mag
                   detection=1
                   if GRB_detection_status=='not': GRB_detection_status='detected'
                   #print (float(2.5*np.log10(1.+1./SNR)),gft_dict['err_min'][band],err_mag)
               elif SNR < SNR4detection :
                   err_mag = np.sqrt((float(2.5*np.log10(1+1./SNR4detection)))**2 + (err_sys)**2)
                   #Compute limiting magnitude   
                   etc_GFT.information['etc_type']='mag'
                   etc_GFT.information['SNR']=SNR4detection
                   etc_GFT.information['Nexp']=1
                   etc_GFT.sim()
                   mag = etc_GFT.information['mag']
                   detection=0
               
               obs_table['mag'][j]=mag
               obs_table['mag_err'][j]=err_mag
               obs_table['detection'][j]=detection
               obs_table['phot_sys'][j]=photometry_system
           #if self.num_samples==1 or p==0: self.observations=obs_table.copy()
           #else: self.observations=vstack([self.observations,obs_table.copy()])

           # Write Light curves in files
           self.format_outputs(obs_table,resdir=self.path_pyGRBz+resdir, fname='', RA_J2000=RA_J2000, DEC_J2000=DEC_J2000, galactic_dust_corrected=1)
           obs_table=None
           plt.close("all")
           gc.collect()
       #print ('observations table:')
       #print (self.observations)
       #self.observations=obs_table

       #ascii.write(self.observations,respath+'/Observations.dat',overwrite=True)
       #ascii.write(self.observations,respath+'/Observations.tex',format='latex',overwrite=True)

       return None

   def sources_extraction(self,image,sextractor_pars):
       """ Extract sources from the generated image using sextractor """

       cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP, gain, pixelScale,seeing,back_type,back_value,back_size,backphoto_type,backphoto_thick,back_filterthresh,checkimage_type,checkimage_name,filter_sex,filter_name= sextractor_pars
       #print ('sextractor')
       #sp.run('sex %s -c gft.sex -CATALOG_NAME %s.cat' % (image,cat_name),shell=True)
       sp.run('sex %s -c /home/dcorre/code/etc_is/gft-sim/obs_strategy/obs_strategy/gft.sex -CATALOG_NAME %s.cat -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME /home/dcorre/code/etc_is/gft-sim/obs_strategy/obs_strategy/gft.param -DETECT_TYPE CCD -DETECT_MINAREA %d -DETECT_THRESH %d -ANALYSIS_THRESH %d -PHOT_APERTURES %d -SATUR_LEVEL %d -MAG_ZEROPOINT %f -GAIN %f -PIXEL_SCALE %f -SEEING_FWHM %f -BACK_TYPE %s -BACK_VALUE %f -BACK_SIZE %d -BACKPHOTO_TYPE %s -BACKPHOTO_THICK %d -BACK_FILTTHRESH %f -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s -FILTER %s -FILTER_NAME %s' % (image,cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP, gain, pixelScale, seeing, back_type, back_value, back_size, backphoto_type, backphoto_thick, back_filterthresh,checkimage_type,checkimage_name,filter_sex, filter_name),shell=True)

   
   def set_lightcurves_images(self,SNR4detection=3,err_mag_sys=[0.04,0.06],nb_strategy=1,configFile='Test_mock.hjson',telescope='gft',etc_plot=False,disp=False,verbose=False,ImageSize=[500,500]):
       """ Computes the observed magnitudes and substitutes it inside the observationnal strategy table

       Parameters
       ----------
       SNR4detection: float (default: 3)
                      SNR for which we consider to have a  detection
 
       err_mag_sys: list (default: [0.04,0.06])
                    systematic error to be added in quadrature with the statistical one.
                    First element is for the visible and second for NIR channel. (values from GROND)

       nb_strategy: integer (default: 1)
                    number of strategy to use: 1 or 2

        
       Returns
       -------
       None. The updated self.observations astropy table is written in .dat and .tex files  
       """
       
       from pyETC.pyETC import etc
       from ImSimpy.ImSimpy import ImageSimulator

       respath_images= self.respath+'/Images'
       respath_is= self.respath+'/IS_lc'
       if not os.path.exists(respath_images): os.makedirs(respath_images)
       if not os.path.exists(respath_is): os.makedirs(respath_is)


       # Parameters for source extraction with Sextractor
       detect_minarea = 1
       detect_thresh = 1.1
       analysis_thresh = 1.
       phot_aperture = 10
       back_type = 'AUTO'
       back_value = 0.0
       back_size = 70
       backphoto_type = 'LOCAL'
       backphoto_thick = 24
       backphoto_filterthresh = 0.0
       checkimage_type = 'APERTURES'
       checkimage_name = 'all_sources'
       satur_level = 65535
       filter_sex = 'Y'
       filter_name = '/home/corre/codes/astromatic/sextractor-2.19.5/config/gauss_1.5_3x3.conv'


       filename='observations'
       #fname_obs = "%s/%s" % (resdir,filename)

       #Load two etc, one with visible detectors and the other one with cagire
       #a_vis=etc(configFile=self.path+'obs_strategy/data/input_vis')
       #a_nir=etc(configFile=self.path+'obs_strategy/data/input_nir')
       #print (sys.getsizeof(a_vis))
       obs_dict=defaultdict(lambda: defaultdict(dict))

       # Place GRB at center of the image
       grb_coord_type='pixels'
       grb_coords=[int(ImageSize[0]/2),int(ImageSize[1]/2)]

       nb_sample=0
       p1 = len(self.grb_params)

       GFT_IS=ImageSimulator(configFile=configFile,name_telescope=telescope)
       GFT_IS.readConfigs()

       GFT_IS.config['verbose']=False
       GFT_IS.config['etc_plot']=False
       GFT_IS.config['disp']=False
       GFT_IS.config['etc_type']='snr'
       GFT_IS.config['object_type']='grb_sim'

       GFT_IS.config['galactic_extinction_law']='none'

       GFT_IS.config['IGM_extinction_model']='meiksin'
       # Resize image to be 500x500 pixels to gain computational speed
       GFT_IS.config['ImageResized']=ImageSize
       
       print ("GRB Time Band")
       for p in range(p1):

           #Initailise detection to  False
           GRB_detection_status='not'

           # Need to make it in a more practical way
           """
           #Reinitialize filters for first strategy in case of 2 startegies
           # Need to make it dynamically
           if nb_strategy==2:
               self.set_filters(vis1=['gri'],vis2=['zy'],nir=['J','H']) 
               self.set_strategy_obs(load_strategy=False,t_since_burst=60+6*3600,t_start=0,t_end=300,t_switch=0,t_dead=5,texp_1=30,texp_2=30,channels=3)
           """
           #print (self.strategy_obs)
           #print ('GRB params:')
           #print (self.grb_params[p])
           obs_table=self.strategy_obs.copy()
           obs_table['Name']=self.grb_params['name'][p]
           n=len(obs_table)

           #GRB properties
           GFT_IS.config['grb_redshift']=self.grb_params['z'][p]
           GFT_IS.config['Av_Host']=self.grb_params['Av_host'][p]
           GFT_IS.config['host_extinction_law']=self.grb_params['ExtLaw'][p]
           if self.model == 'gs02':
               GFT_IS.config['grb_model']='gs02'
               GFT_IS.config['n0']=self.grb_params['n0'][p]
               GFT_IS.config['eps_b']=self.grb_params['eps_b'][p]
               GFT_IS.config['eps_e']=self.grb_params['eps_e'][p]
               GFT_IS.config['E_iso']=self.grb_params['E_iso'][p]
               GFT_IS.config['p']=self.grb_params['p'][p]
               GFT_IS.config['eta']=self.grb_params['eta'][p]
               GFT_IS.config['Y']=self.grb_params['Y'][p]
               GFT_IS.config['ism_type']=self.grb_params['ism_type'][p]
               GFT_IS.config['Mdot_loss']=self.grb_params['Mdot_loss'][p]
               GFT_IS.config['Vw']=self.grb_params['Vw'][p]
           elif self.model == 'kann':
               GFT_IS.config['grb_model']='SPL'
               GFT_IS.config['F0']=self.grb_params['F0'][p]
               GFT_IS.config['wvl0']=self.grb_params['wvl0'][p]
               GFT_IS.config['t0']=self.grb_params['t0'][p]
               GFT_IS.config['alpha']=self.grb_params['alpha'][p]
               GFT_IS.config['beta']=self.grb_params['beta'][p]
           #print (GFT_IS.config)
           for j in range(n):

               """
               #Change of strategy if GRB detected 
               if GRB_detection_status == 'detected':
                   print ('GRB detected')
                   #print (obs_table['Name','Obs_time','time_since_burst','band'])
                   self.set_filters(vis1=['g','r','i'],vis2=['z','y'],nir=['J','H'])
                   #print (obs_table['time_since_burst'][j+1],obs_table['Obs_time'][j+1])
                   counter=0
                   #print (obs_table['Nb_obs'][j-1])
                   while obs_table['Nb_obs'][j-1] == obs_table['Nb_obs'][j+counter]:
                       #print (obs_table['Nb_obs'][j+counter])
                       counter+=1
                   #print (counter)
                   #print (obs_table['time_since_burst'][j+counter],obs_table['Obs_time'][j+counter])
                   self.set_strategy_obs(load_strategy=False,t_since_burst=obs_table['time_since_burst'][j-1],t_start=obs_table['Obs_time'][j+counter],t_end=300,t_switch=0,t_dead=5,texp_1=30,texp_2=30,channels=3)
                   #print ('test')
                   #print (self.strategy_obs['Name','Obs_time','time_since_burst','band'])
                   obs_table[j+counter:]=self.strategy_obs.copy()
                   obs_table['Name']=self.grb_params['name'][p]
                   GRB_detection_status='already detected'
                   #print (obs_table['Name','Obs_time','time_since_burst','band'])
               """
               nb_sample+=1

               """
               etc_GFT=etc(configFile=configFile,name_telescope=telescope)
               etc_GFT.information['etc_type']='snr'
               etc_GFT.information['object_type']='grb_sim'
               etc_GFT.information['Nexp']=1
               etc_GFT.information['etc_plot']=etc_plot
               etc_GFT.information['disp']=disp
               etc_GFT.information['verbose']=verbose

               #GRB properties
               etc_GFT.information['grb_redshift']=self.grb_params['z'][p]
               etc_GFT.information['Av_Host']=self.grb_params['Av_host'][p]
               if self.model == 'gs02':
                   etc_GFT.information['grb_model']='gs02'
                   etc_GFT.information['n0']=self.grb_params['n0'][p]
                   etc_GFT.information['eps_b']=self.grb_params['eps_b'][p]
                   etc_GFT.information['eps_e']=self.grb_params['eps_e'][p]
                   etc_GFT.information['E_iso']=self.grb_params['E_iso'][p]
                   etc_GFT.information['p']=self.grb_params['p'][p]
                   etc_GFT.information['eta']=self.grb_params['eta'][p]
                   etc_GFT.information['Y']=self.grb_params['Y'][p]
                   etc_GFT.information['ism_type']=self.grb_params['ism_type'][p]
                   etc_GFT.information['Mdot_loss']=self.grb_params['Mdot_loss'][p]
                   etc_GFT.information['Vw']=self.grb_params['Vw'][p]
               elif self.model == 'kann':
                   etc_GFT.information['grb_model']='SPL'
                   etc_GFT.information['F0']=self.grb_params['F0'][p]
                   etc_GFT.information['wvl0']=self.grb_params['wvl0'][p]
                   etc_GFT.information['t0']=self.grb_params['t0'][p]
                   etc_GFT.information['alpha']=self.grb_params['alpha'][p]
                   etc_GFT.information['beta']=self.grb_params['beta'][p]
               """

               GFT_IS.config['etc_type']='snr'
               GFT_IS.config['Nexp']=1
               print (p,obs_table['time_since_burst'][j],obs_table['band'][j])
               if obs_table['band'][j] in self.filter_vis1:
                   GFT_IS.config['channel']='DDRAGO-B'
                   GFT_IS.config['GainMapFile']='%s/Gain_vis.fits' % self.resdir
                   GFT_IS.config['VignettingFile']='%s/Vignetting_vis.fits' % self.resdir
                   GFT_IS.config['OffsetFile']='%s/Offset_vis.fits' % self.resdir
                   GFT_IS.config['DeadPixFile']='%s/DeadPixs_vis.fits' % self.resdir
                   GFT_IS.config['HotPixFile']='%s/HotPixs_vis.fits' % self.resdir
                   err_sys=err_mag_sys[0]
               elif obs_table['band'][j] in self.filter_vis2:
                   GFT_IS.config['channel']='DDRAGO-R'
                   GFT_IS.config['GainMapFile']='%s/Gain_vis.fits' % self.resdir
                   GFT_IS.config['VignettingFile']='%s/Vignetting_vis.fits' % self.resdir
                   GFT_IS.config['OffsetFile']='%s/Offset_vis.fits' % self.resdir
                   GFT_IS.config['DeadPixFile']='%s/DeadPixs_vis.fits' % self.resdir
                   GFT_IS.config['HotPixFile']='%s/HotPixs_vis.fits' % self.resdir
                   err_sys=err_mag_sys[0]
               elif obs_table['band'][j] in self.filter_nir:
                   GFT_IS.config['channel']='CAGIRE'
                   GFT_IS.config['GainMapFile']='%s/Gain_nir.fits' % self.resdir
                   GFT_IS.config['VignettingFile']='%s/Vignetting_nir.fits' % self.resdir
                   GFT_IS.config['OffsetFile']='%s/Offset_nir.fits' % self.resdir
                   GFT_IS.config['DeadPixFile']='%s/DeadPixs_nir.fits' % self.resdir
                   GFT_IS.config['HotPixFile']='%s/HotPixs_nir.fits' % self.resdir
                   err_sys=err_mag_sys[1]

               GFT_IS.config['exptime']=obs_table['Exp_time'][j]
               GFT_IS.config['filter_band']=obs_table['band'][j]
               GFT_IS.config['Fwhm_psf_opt']=self.optics_fwhm[obs_table['band'][j]]
               GFT_IS.config['t_sinceBurst']=obs_table['time_since_burst'][j]/86400
               GFT_IS.config['SourcesList']['file']="%s/Sources.txt" % (self.resdir)#,obs_table['band'][j])
               GFT_IS.config['output']='%s/%s/%s_%ss.fits' % (self.resdir,obs_table['Name'][0],obs_table['band'][j],obs_table['time_since_burst'][j])
               GFT_IS.config['grb_coord_type']=grb_coord_type
               GFT_IS.config['grb_coords']=grb_coords
               GFT_IS.config['PSF']['total']['method']='load'
               GFT_IS.config['PSF']['total']['file']='total_PSF/%s/PSF_total_%s.fits' % (self.resdir,obs_table['band'][j])
               
               #print (GFT_IS.config)
               GFT_IS.simulate('data')
               #print (GFT_IS.information)
               #mag=GFT_IS.information['mag']
               #SNR=GFT_IS.information['SNR']
               photometry_system=GFT_IS.information['photometry_system']
               #del etc
               gc.collect()

               #print (GFT_IS.information['grb_mag'])
               # Run sextractors
               gain=GFT_IS.information['gain']
               ZP=GFT_IS.information['zeropoint'] -2.5*np.log10(gain)
               pixelScale=GFT_IS.information['pixelScale_X']
               well_cap=GFT_IS.information['FWC']
               EXPTIME=GFT_IS.information['exptime']
               seeing=GFT_IS.information['seeing']

               ZP_sextractor = ZP + 2.5*np.log10(EXPTIME)

               checkimage_path = GFT_IS.path + '/images/'+'%s/%s/%s_%ss_%s.fits' % (self.resdir,obs_table['Name'][0],obs_table['band'][j],obs_table['time_since_burst'][j],checkimage_name)
               detected_cat_name='%s/%s_%s_%ss' % (respath_images,obs_table['Name'][0],obs_table['band'][j],obs_table['time_since_burst'][j])

               sextractor_pars = detected_cat_name, detect_minarea, detect_thresh, analysis_thresh, phot_aperture, satur_level, ZP_sextractor, gain, pixelScale, seeing, back_type, back_value, back_size, backphoto_type, backphoto_thick, backphoto_filterthresh, checkimage_type, checkimage_path,filter_sex, filter_name
               #print (GFT_IS.path + '/images/'+GFT_IS.config['output'])
               self.sources_extraction(GFT_IS.path + '/images/'+GFT_IS.config['output'],sextractor_pars)
               #self.sources_extraction('g_60s.fits',sextractor_pars)


               # Get magnitude from sextractor output
               sex_cat=ascii.read('%s.cat' % detected_cat_name)
               # Get the source at GRB position
               offset=0
               
               mask = (sex_cat['XPEAK_IMAGE'] >= 150-offset) & (sex_cat['XPEAK_IMAGE'] <= 150+offset) & (sex_cat['YPEAK_IMAGE'] >= 150-offset) & (sex_cat['YPEAK_IMAGE'] <= 150+offset)
               
               if mask.any():
                   #print (sex_cat[mask])
                   print ('GRB_mag: %.2f' % GFT_IS.information['grb_mag'])
                   print ('ISO: %.2f +/- %.2f' % (sex_cat[mask]['MAG_ISO'],sex_cat[mask]['MAGERR_ISO']))
                   print ('ISOCORR: %.2f +/- %.2f' % (sex_cat[mask]['MAG_ISOCOR'],sex_cat[mask]['MAGERR_ISOCOR']))
                   print ('APER: %.2f +/- %.2f' % (sex_cat[mask]['MAG_APER'],sex_cat[mask]['MAGERR_APER']))

                   print ('AUTO: %.2f +/- %.2f' % (sex_cat[mask]['MAG_AUTO'],sex_cat[mask]['MAGERR_AUTO']))
                   print ('SNR: %.2f' % (sex_cat[mask]['FLUX_AUTO']/sex_cat[mask]['FLUXERR_AUTO']))

                   mag, mag_err = sex_cat[mask]['MAG_ISOCOR'], sex_cat[mask]['MAGERR_ISOCOR']
                   # Add quadratically systematic uncertainty 
                   err_mag = np.sqrt((mag_err)**2 + (err_sys)**2)
                   detection=1
                   if GRB_detection_status=='not': GRB_detection_status='detected'
                   #print (float(2.5*np.log10(1.+1./SNR)),gft_dict['err_min'][band],err_mag)
               else:
                   err_mag = np.sqrt((float(2.5*np.log10(1+1./SNR4detection)))**2 + (err_sys)**2)
                   #Compute limiting magnitude 
                   mag=24  
                   #etc_GFT.information['etc_type']='mag'
                   #etc_GFT.information['SNR']=SNR4detection
                   #etc_GFT.information['Nexp']=1
                   #etc_GFT.sim()
                   #mag = etc_GFT.information['mag']
                   detection=0
               
               obs_table['mag'][j]=mag
               obs_table['mag_err'][j]=err_mag
               obs_table['detection'][j]=detection
               obs_table['phot_sys'][j]=photometry_system

           if self.num_samples==1 or p==0: self.observations=obs_table.copy()
           else: self.observations=vstack([self.observations,obs_table.copy()])
           obs_table=None
           #del obs
           #del 
           plt.close("all")
           gc.collect()
       print ('observations table:')
       print (self.observations)
       #self.observations=obs_table

       #ascii.write(self.observations,respath_etc+'/Observations.dat',overwrite=True)
       #ascii.write(self.observations,respath_etc+'/Observations.tex',format='latex',overwrite=True)

       return None

#   def format_outputs(self,resdir=os.getenv('GFTSIM_DIR')+'/obs_strategy/obs_strategy/results',fname='_pZ_format',RA_J2000='0h54m50.6s',DEC_J2000='+14d05m04.5s',galactic_dust_corrected=False):
#       """ formatting the observations to be used by photoz code """
#
#       if not os.path.exists(resdir): os.makedirs(resdir)
#
#       for obs in self.observations.group_by('Name').groups:
#           text="""#name:%s
##type:lc
##RA_J2000:%s
##DEC_J2000:%s
##MW_corrected:%s
##time_unit:s 
##z:%s
##Av_host:%s
#time_since_burst band mag mag_err zp phot_sys detection telescope
#""" % (obs['Name'].data[0],str(RA_J2000),str(DEC_J2000),galactic_dust_corrected,self.grb_params['z'][self.grb_params['name'] == obs['Name'].data[0]].data[0],self.grb_params['Av_host'][self.grb_params['name'] == obs['Name'].data[0]].data[0])
#
#           f = open('%s/%s%s.txt'  % (resdir,obs['Name'].data[0],fname), 'w')
#           f.write(text)
#           for line in obs:
#               f.write('%.2f %s %.4f %.4f %s %s %d %s\n' % (line['time_since_burst'],line['band'],line['mag'],line['mag_err'],'-',line['phot_sys'],line['detection'],line['telescope']))
#           f.close()
#  
#       return None

   def format_outputs(self,observations,resdir='/results',fname='_pZ_format',RA_J2000='0h54m50.6s',DEC_J2000='+14d05m04.5s',galactic_dust_corrected=1):
       """ formatting the observations to be used by photoz code """
       if not os.path.exists(resdir): os.makedirs(resdir)

       for obs in observations.group_by('Name').groups:
           text="""#name:%s
#type:lc
#RA_J2000:%s
#DEC_J2000:%s
#MW_corrected:%s
#time_unit:s 
#z:%s
#Av_host:%s
time_since_burst band mag mag_err zp phot_sys detection telescope
""" % (obs['Name'].data[0],str(RA_J2000),str(DEC_J2000),galactic_dust_corrected,self.grb_params['z'][self.grb_params['name'] == obs['Name'].data[0]].data[0],self.grb_params['Av_host'][self.grb_params['name'] == obs['Name'].data[0]].data[0])

           f = open('%s/%s%s.txt'  % (resdir,obs['Name'].data[0],fname), 'w')
           f.write(text)
           for line in obs:
               f.write('%.2f %s %.4f %.4f %s %s %d %s\n' % (line['time_since_burst'],line['band'],line['mag'],line['mag_err'],'-',line['phot_sys'],line['detection'],line['telescope']))
           f.close()
  
       return None

if __name__ == '__main__':
   
   a=obs_strategy()
   a.set_filters(vis1=['g','r','i'],vis2=['z','y'],nir=['J','H'])
   
   a.set_strategy_obs(load_strategy=False,t_since_burst=60,t_start=0,t_end=300,t_switch=0,t_dead=5,texp_1=30,texp_2=30,channels=3)
   params=dict(z=['truncnorm',5,2.7,0,10],
            Av_host=['truncnorm',0.2,0.3,0,2],
            E_iso=['constant',56],
            eta=['constant',0.3],
            eps_b=['constant',1e-4],
            eps_e=['constant',0.1],
            p=['constant',2.3],
            Y=['constant',0],
            ism_type=['constant',0],
            n0=['constant',1],
            Mdot_loss=['constant',1e-5],
            Vw=['constant',1000]
            )
   a.set_grb_params(params=params,
                 num_samples=300,random_flag=False,model='gs02',
                 #load_params=True,params_file='distrib/Granot_Sari/samples300_gs02.tex',load_min=1,load_max=100,
                 load_params=True,params_file='test_wp2.tex',load_min=201,load_max=300,
                 load_distrib=False,distrib_file='Granot_Sari/Turpin_GS02_v2_Ek.txt',
                 plot_distribution=True,
                 seed=-1)
   """
   a.set_strategy_obs(load_strategy=False,t_since_burst=60,t_start=0,t_end=300,t_switch=0,t_dead=5,texp_1=30,texp_2=30,channels=3)
   params=dict(z=['truncnorm',5,3,0.01,10],
            Av_host=['truncnorm',0.2,0.3,0.01,2],
            wvl0=['constant',6400],
            t0=['constant',86.4],
            Rc=['constant',12],
            #alpha=['truncnorm',0.48,0.69,0.01,3],
            alpha=['truncnorm',1.2,0.6,0.01,3],
            beta=['truncnorm',0.66,0.25,0.01,2]
            )
   a.set_grb_params(params=params,
                 num_samples=300,random_flag=False,model='kann',
                 load_params=True,params_file='distrib/Kann/samples300_kann.tex',load_min=1,load_max=300,
                 load_distrib=False,distrib_file='Kann/Rc_86s_z=1.txt',
                 plot_distribution=True,
                 seed=-1)

   """
   a.set_lightcurves(SNR4detection=3,err_mag_sys=[0.04,0.06],nb_strategy=1)
   a.format_outputs(resdir='/data/lc/',fname='test_wp')
