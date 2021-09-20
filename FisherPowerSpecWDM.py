#----------------------------------
# Code to compute the Fisher matrix from 21cmFAST simulations with WDM
# Author: Pablo Villanueva Domingo
# Started July 2020
# For details on the Fisher matrix computation, see e.g. https://arxiv.org/PS_cache/arxiv/pdf/0906/0906.4123v1.pdf
#----------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import glob, time


time_ini = time.time()

#--- PARAMETERS ---#

# Path of the 21cmFAST simulations
pathsims = "/projects/QUIJOTE/WDMSimulationsFisher/"
# Step for computing the derivatives
step = 0.1
# Names of the parameters
paramnames = [r"$log_{10}(M_{\rm turn})$", r"$log_{10}(L_{\rm X})$", r"$log_{10}(N_{\gamma})$", r"$1/m_{\rm wdm}$"]
# Simulation numbers relative to each parameter
params = [1,3,5,7]
n_params = len(params)
# If 1, choose a given set of redshifts, otherwise chooses several z
choose_z = 1
# 1 for plotting derivatives
plot_deriv = 1
# Number of random seeds
n_seeds = 30
list_seeds = list(range(n_seeds))
# Number of random seeds for computing the covariance matrix
n_seeds_cov = 300

# Select redshifts
if choose_z:
    redshifts = ["010.16"]#,"015.78","020.18"]
else:
    filefolder = pathsims+"Simulation_seed_0_num_1/Output_files/Deldel_T_power_spec/"
    folder = glob.glob(filefolder+"*")
    redshifts = []
    for fil in folder:
        redshifts.append(fil[len(filefolder)+4:len(filefolder)+10])
    redshifts.sort()
    #redshifts = redshifts[8:]
    redshifts = redshifts[8:8+12]

print("Redshifts considered:", redshifts)

n_redshifts = len(redshifts)

# 1, 2, 3 sigma CL legend
cols = ["cyan", "blue", "purple"]
lines = ["-", ":"]
customlegend = []
for i, col in enumerate(cols):
    customlegend.append(Line2D([0], [0], color=col, lw=3., linestyle="-", label=str(i+1)+r"$\sigma$"))

#--- FUNCTIONS ---#

# Load the 21 cm power spectrum from a simulation number numsim with a given random seed
def load_ps(seed, numsim):

    powsz = []
    filefolder = pathsims+"Simulation_seed_"+str(seed)+"_num_"+str(numsim)+"/Output_files/Deldel_T_power_spec/"
    folder = glob.glob(filefolder+"*")

    folder.sort()
    if len(folder)==0:
        print("Folder "+filefolder+" empty")
        return 0
    else:
        for z in redshifts:
            file = glob.glob(filefolder+"/ps_z"+z+"*")
            tabpow = np.loadtxt(file[0],unpack=True)
            powk = tabpow[1]
            #powk = np.log10(tabpow[1])
            powsz.append(powk)

        powsz = np.array(powsz)
        powsz = powsz.reshape(-1)
        return powsz

# Compute the derivative with respect to a given parameter given by numsim
def compute_deriv(numsim, step, forward_derivative = False):

    ps_plus_ar, ps_minus_ar = [], []
    for seed in list_seeds:
        ps_plus = load_ps(seed, numsim)
        ps_minus = load_ps(seed, numsim+1)
        ps_plus_ar.append(ps_plus); ps_minus_ar.append(ps_minus)
    ps_plus_ar, ps_minus_ar = np.array(ps_plus_ar), np.array(ps_minus_ar)
    ps_plus_mean, ps_minus_mean = ps_plus_ar.mean(axis=0), ps_minus_ar.mean(axis=0)
    deriv = (ps_plus_mean - ps_minus_mean)/2./step
    if forward_derivative:
        ps_fids = []
        for seed in list_seeds:
            ps_fids.append( load_ps(seed, 0) )
        ps_fid_mean = np.array(ps_fids).mean(axis=0)
        deriv = (4.*ps_minus_mean - ps_plus_mean -3.*ps_fid_mean)/2./step

    return deriv

# Compute the covariance matrix
def covariance_matrix():

    ps_ar = []
    for seed in range(n_seeds_cov):
        ps = load_ps(seed, 0)   # number of simulation 0 has parameters at the fiducial values
        ps_ar.append(ps)
    ps_ar = np.array(ps_ar)
    ps_mean = ps_ar.mean(axis=0)
    ps_centered = ps_ar - ps_mean
    ps_centered = np.transpose(ps_centered)
    covariance = np.cov(ps_centered)

    # Compute the correlation matrix
    corr = np.zeros((len(ps_mean),len(ps_mean)))
    for a in range(len(ps_mean)):
        for b in range(len(ps_mean)):
            corr[a,b] = covariance[a,b]/np.sqrt(covariance[a,a]*covariance[b,b])

    return covariance, corr

# Compute the Fisher matrix given a covariance matrix
def fisher_matrix(cov_matrix, cov_diag=False):

    n_modes = cov_matrix.shape[0]
    derivates = []

    # Compute derivatives for each parameter
    for param in params:
        # Compute forward derivative for the wdm mass (simulations 7 and 8)
        if param==7:
            for_der = True
        else:
            for_der = False
        deriv = compute_deriv(param, step, forward_derivative = for_der)
        derivates.append(deriv)

    # Plot derivatives
    if plot_deriv:
        fig, ax = plt.subplots(figsize = (4,4))
        for param in range(n_params):
            ax.plot(derivates[param],label=paramnames[param])
        ax.legend()
        fig.savefig("Plots/derivatives_nredshifts_{:d}_nseedsder_{:d}.pdf".format(n_redshifts, n_seeds), bbox_inches='tight')
        plt.close(fig)

    # Compute the Fisher matrix
    fisher = np.zeros((n_params,n_params))
    invcov = np.linalg.inv(cov_matrix)

    for i in range(n_params):
        for j in range(n_params):
            if cov_diag:
                for a in range(n_modes):
                    #fisher[i,j] += derivates[i][a]*derivates[j][a]/cov_matrix[a,a]
                    fisher[i,j] += derivates[i][a]*derivates[j][a]*invcov[a,a]
            else:
                for a in range(n_modes):
                    for b in range(n_modes):
                        fisher[i,j] += derivates[i][a]*derivates[j][b]*invcov[a,b]
            #fisher[i,j] = np.sum(derivates[i]*derivates[j]/cov_matrix)

    return fisher

# Marginalize over all parameters except 2
def marginalize(cov_params, parx, pary):
    margin_cov = np.zeros((2,2))
    margin_cov[0,0] = cov_params[parx,parx]
    margin_cov[1,1] = cov_params[pary,pary]
    margin_cov[0,1] = cov_params[parx,pary]
    margin_cov[1,0] = margin_cov[0,1]
    return margin_cov

# Compute the chi2 for two parameters
def chisquare(x, y, cov):
    sig2_x, sig2_y, sig_xy = cov[0,0], cov[1,1], cov[0,1]
    rho = sig_xy/np.sqrt(sig2_x*sig2_y)
    chi2 = (x**2./sig2_x + y**2./sig2_y -2.*rho*x*y/np.sqrt(sig2_x*sig2_y))/(1.-rho**2.)
    return chi2


#--- MAIN ---#

cov_matrix, corr_matrix = covariance_matrix()

"""
print("Determinant of the covariance matrix", np.linalg.det(cov_matrix))
v1,w1 = np.linalg.eig(cov_matrix)
print('Max eigenvalue    Cov = %.3e'%np.max(v1))
print('Min eigenvalue    Cov = %.3e'%np.min(v1))
print('Condition number  Cov = %.3e'%(np.max(v1)/np.min(v1)))
print(' ')"""

fig_cov, ax = plt.subplots(figsize = (4,4))
ax.imshow(corr_matrix, cmap="coolwarm")
fig_cov.savefig("Plots/correlation_matrix_nredshifts_{:d}_nseedscov_{:d}.pdf".format(n_redshifts, n_seeds_cov))
plt.close(fig_cov)

fig = plt.figure(figsize=(n_params*2,n_params*2),constrained_layout=True)
spec = fig.add_gridspec(ncols=n_params-1, nrows=n_params-1, wspace=0., hspace=0.)

# Compute the Fisher matrix either assuming a diagonal covariance matrix or a complete one
for cov_diag in [0,1]:

    fisher = fisher_matrix(cov_matrix, cov_diag=cov_diag)

    # The inverse of the Fisher matrix gives the covariance matrix of the parameters
    cov_params = np.linalg.inv(fisher)

    for par1 in range(n_params):
        for par2 in range(par1+1,n_params):

            ax = fig.add_subplot(spec[par2-1, par1])

            marg_cov = marginalize(cov_params, par1, par2)

            sig2_x, sig2_y = marg_cov[0,0], marg_cov[1,1]

            xvec = np.linspace(-4.*np.sqrt(sig2_x), 4.*np.sqrt(sig2_x), num=100)
            yvec = np.linspace(-4.*np.sqrt(sig2_y), 4.*np.sqrt(sig2_y), num=100)
            xx, yy = np.meshgrid(xvec, yvec)

            chi2 = chisquare(xx, yy, marg_cov)

            #ax.contourf(xvec, yvec, chi2, levels=[0., 2.3, 6.17, 11.8], colors = cols)
            ax.contour(xvec, yvec, chi2, levels=[2.3, 6.17, 11.8], colors = cols, linestyles = lines[cov_diag])
            #ax.set_aspect('equal', 'datalim')
            #ax.legend(handles=customlegend)
            ax.set_xlabel(paramnames[par1])
            ax.set_ylabel(paramnames[par2])


fig.savefig("Plots/chi2_fisher_nredshifts_{:d}_nseedsderiv_{:d}_nseedscov_{:d}".format(n_redshifts, n_seeds, n_seeds_cov)+".pdf", bbox_inches='tight')


print("Minutes elapsed:",(time.time()-time_ini)/60.)
