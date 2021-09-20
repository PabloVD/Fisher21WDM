#------------------------
# Copy the code 21cmFAST, which performs the simulations, with the relevant values of the parameters to compute derivatives and the Fisher matrix
# It needs a copy of 21cmFAST in the main directory
# Author: Pablo Villanueva Domingo
# Started July 2020
#------------------------

import os
import fileinput
import sys
import numpy as np
import time
import random

time_ini = time.time()

##################################### PARAMETERS ##########################################

# Path of the simulations
pathsims = "/projects/QUIJOTE/WDMSimulationsFisher/"

# Fiducial values for the parameters
mturn_fid = 8.
lx_fid = 40.
ngamma_fid = 4.
oneovermwdm_fid = 0.

# Step
mturn_step = 0.1
lx_step = 0.1
ngamma_step = 0.1
oneovermwdm_step = 0.1

# Number of random seeds and initial seed
n_seeds = 30
ini_seed = 0
seeds = range(ini_seed, n_seeds)

# Parameter grid
"""mturn_vals = np.linspace(7,10,num=num_sim).tolist()
lx_vals = np.linspace(38.,42.,num=num_sim).tolist()
ngamma_vals = np.linspace(3,5,num=num_sim).tolist()
mwdm_vals = np.linspace(1.,10.,num=num_sim).tolist()"""

# Parameter files paths
heat_file = '/Parameter_files/HEAT_PARAMS.H'
anal_file = '/Parameter_files/ANAL_PARAMS.H'
init_file = '/Parameter_files/INIT_PARAMS.H'
cosmo_file = '/Parameter_files/COSMOLOGY.H'

######################################################################################


##################################### FUNCTIONS #####################################

# Replace a string in a given file
def ReplaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = replaceExp
        sys.stdout.write(line)

# Copy the code with the proper parameters
def Routine(mturn,lx,ngamma,mwdm,seed,i,cutoff=False):

    print("Seed", seed, "Simulation", i)
    print("mturn_{:.3f}_lumx_{:.3f}_ngamma_{:.3f}_mwdm_{:.3f}_seed_{:d}".format(mturn,lx,ngamma,mwdm,seed))

    foldername = pathsims+"Simulation_seed_"+str(seed)+"_num_"+str(i)
    if os.path.exists(foldername):
        print(foldername+" already exists!")
    os.system("cp -R 21cmFAST-master "+foldername)

    ReplaceAll(foldername + anal_file, "#define M_TURNOVER", "#define M_TURNOVER (double) (pow(10.,{:.3f})) // Halo mass threshold for efficient star formation (in Msun). \n".format(mturn))
    ReplaceAll(foldername + anal_file, "#define N_GAMMA_UV", "#define N_GAMMA_UV (float) (pow(10.,{:.3f})) // number of ionizing photons per stellar baryon \n".format(ngamma))
    ReplaceAll(foldername + heat_file, "#define L_X", "#define L_X (double) ({:.3f}) \n".format(lx))
    if cutoff:
        ReplaceAll(foldername + cosmo_file, "#define P_CUTOFF", "#define P_CUTOFF (int) (1) // supress the power spectrum? 0= CDM; 1=WDM \n")
        ReplaceAll(foldername + cosmo_file, "#define M_WDM", "#define M_WDM (float) ({:.3f}) // mass of WDM particle in keV.  this is ignored if P_CUTOFF is set to zero \n".format(mwdm))
    else:
        ReplaceAll(foldername + cosmo_file, "#define P_CUTOFF", "#define P_CUTOFF (int) (0) // supress the power spectrum? 0= CDM; 1=WDM \n")
    ReplaceAll(foldername + init_file, "#define RANDOM_SEED", "#define RANDOM_SEED (long) ({:d}) // seed for the random number generator \n".format(seed))

    # Write file with the parameters
    np.savetxt(foldername+"/params_mturn_lumx_ngamma_mwdm_seed.txt",np.transpose([[mturn],[lx],[ngamma],[mwdm],[seed]]))



######################################################################################


##################################### MAIN ##########################################

if __name__ == "__main__":

    if not os.path.exists(pathsims):
        os.system("mkdir "+pathsims)

    mturn, lx, ngamma, mwdm = mturn_fid, lx_fid, ngamma_fid, 100.

    # For each random seed, copy the code with the fiducial values and a step above/below for each parameter, for computing the derivatives for the Fisher matrix
    for seed in seeds:

        mturn, lx, ngamma, mwdm = mturn_fid, lx_fid, ngamma_fid, 100.

        Routine(mturn,lx,ngamma,mwdm,seed,0)

        Routine(mturn+mturn_step,lx,ngamma,mwdm,seed,1)
        Routine(mturn-mturn_step,lx,ngamma,mwdm,seed,2)

        Routine(mturn,lx+lx_step,ngamma,mwdm,seed,3)
        Routine(mturn,lx-lx_step,ngamma,mwdm,seed,4)

        Routine(mturn,lx,ngamma+ngamma_step,mwdm,seed,5)
        Routine(mturn,lx,ngamma-ngamma_step,mwdm,seed,6)

        Routine(mturn,lx,ngamma,1./( oneovermwdm_fid + oneovermwdm_step ),seed,7,cutoff=True)
        Routine(mturn,lx,ngamma,1./( oneovermwdm_fid + 2.*oneovermwdm_step ),seed,8,cutoff=True)


    print("Minutes elapsed:",(time.time()-time_ini)/60.)

######################################################################################
