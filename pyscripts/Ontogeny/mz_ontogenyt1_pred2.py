import os
import tempfile
from cmdstanpy import CmdStanModel
import pickle as pkl
# Specify Stan program file
stan_file = os.path.join('/Users/apoorvasingh/Desktop/MZB-new-analysis/stan_models/Ontogeny/DDM Influx/T2 Precursor/DDdivision_independent.stan')
# stan_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/stan_models/new_models/Ontogeny/Mz_ontogenyT1.divisionTDM.stan')

# instantiate the model object
model = CmdStanModel(stan_file=stan_file)

# inspect model object
print(model)

# inspect compiled model
print(model.exe_info())

#specify data file
data_file = os.path.join('/Users/apoorvasingh/Desktop/MZB-new-analysis/data/MZcounts_ontogeny.json')
# init_dict = {'beta': 3.5, 'y0_Log': 14, 'kappa_0': 0.5, 'rho': 0.01, 'delta': 0.09, 'sigma': [1,1,1,1]}

#model fit with fixed parameters

#fit = model.sample(data=data_file, inits= init_dict,output_dir='/home/apoorva/Desktop/work/MZ_analysis/output',fixed_param=True ,show_console=True)
fit = model.sample(data=data_file, output_dir='/Users/apoorvasingh/Desktop/MZB-new-analysis/new_output',show_console=True, iter_sampling= 2000, iter_warmup= 100, chains=5, threads_per_chain=12, max_treedepth=20,adapt_delta = 0.99)

# Print summary of the fit
print(fit.summary())

# Diagnose the fit
diagnostic_output = fit.diagnose()
print(diagnostic_output)

#save output to pickle file
with open('/Users/apoorvasingh/Desktop/MZB-new-analysis/data/New output (pkl)/Ontogeny Output/DDM/T2 Precursor/DDM_T2_Division_independent.pkl', 'wb') as f:
    pkl.dump(fit, f)
    
# Use diagnose method to check for any potential problems with the fit
print(fit.diagnose())