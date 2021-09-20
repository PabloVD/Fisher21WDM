#------------------------
# Writes and runs a sbatch script to perform the simulations
# It needs the folders already created with CreateFisherSimulations
# Author: Pablo Villanueva Domingo
# Started July 2020
#------------------------


import numpy as np
import sys,os,glob
from CreateFisherSimulations import *

seeds = range(ini_seed, n_seeds)
sims = range(0,9)

for seed in seeds:
    for i in sims:

        text="""#!/bin/bash
#SBATCH -J %d_%d
#SBATCH --time=2:00:00
######SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8G       # memory per cpu-core (4G per cpu-core is default)

cd Programs
time ./drive_logZscroll_Ts

cd ../Boxes
shopt -s extglob
rm -r -- !(updated_smoothed_deltax*|delta_T_v3_z*)  ## remove boxes files except delta and deltaT_b

"""%(seed,i)

        f = open(pathsims+"Simulation_seed_"+str(seed)+"_num_"+str(i)+"/script.sh",'w');  f.write(text);  f.close()

        os.chdir(pathsims+"Simulation_seed_"+str(seed)+"_num_"+str(i))
        if not(os.path.exists('Boxes/updated_smoothed_deltax_z006.00_200_300Mpc')):
            print("Seed", seed, "Simulation", i)
            os.chdir(pathsims+"Simulation_seed_"+str(seed)+"_num_"+str(i)+"/Programs")
            os.system('make clean; make')
            os.chdir(pathsims+"Simulation_seed_"+str(seed)+"_num_"+str(i))
            os.system('sbatch script.sh')
        else:
            print("Seed", seed, "Simulation", i," already exists")
