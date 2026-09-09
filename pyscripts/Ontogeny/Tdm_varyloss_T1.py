import os
from cmdstanpy import CmdStanModel
import pickle as pkl
# Specify Stan program file
stan_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/stan_models/new_models/Ontogeny/TDM Loss/Mz_ontogenyFO.loss.stan')

# instantiate the model object
model = CmdStanModel(stan_file=stan_file)

# inspect model object
print(model)

# inspect compiled model
print(model.exe_info())

#specify data file

data_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/data/New_data/MZcounts_ontogeny.json')
# init_dict = {'beta': 3.5, 'y0_Log': 14, 'kappa_0': 0.5, 'rho': 0.01, 'delta': 0.09, 'sigma': [1,1,1,1]}

#model fit with fixed parameters

#fit = model.sample(data=data_file, inits= init_dict,output_dir='/home/apoorva/Desktop/work/MZ_analysis/output',fixed_param=True ,show_console=True)
fit = model.sample(data=data_file, output_dir='/home/apoorva/Desktop/work/MZ_analysis/output',show_progress=True,
    show_console=False, iter_sampling= 1000, iter_warmup= 100, chains=5, parallel_chains=5, threads_per_chain=12,  max_treedepth=20,
    adapt_delta=0.99)

# Print summary of the fit
print(fit.summary())

# Diagnose the fit
diagnostic_output = fit.diagnose()
print(diagnostic_output)

#save output to csv file

# fit.save_csvfiles('/home/apoorva/Desktop/work/MZ_analysis/data/fit_T1.csv')


with open('/home/apoorva/Desktop/work/MZ_analysis/data/New output (pkl)/Ontogeny Output/TDM_loss_FO.pkl', 'wb') as f:
    pkl.dump(fit, f)
    
# print(fit.stan_variables())

# Use diagnose method to check for any potential problems with the fit
print(fit.diagnose())