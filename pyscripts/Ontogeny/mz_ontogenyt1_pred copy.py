import os
from cmdstanpy import CmdStanModel
import pickle as pkl
# Specify Stan program file
stan_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/stan_models/new_models/Ontogeny/DDM Influx/T2 Precursor/logistic_growth_varyinginflux_check.stan')
# stan_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/stan_models/new_models/Ontogeny/Mz_ontogenyT1.divisionTDM.stan')

# instantiate the model object
model = CmdStanModel(stan_file=stan_file)

# inspect model object
print(model)

# inspect compiled model
print(model.exe_info())

#specify data file

data_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/data/MZcounts_ontogeny.json')
init_dict = {'phi': 0.01, 'M_stable': 2485012, 'K': 2285012, 'sigma': 0.5}

#model fit with fixed parameters

#fit = model.sample(data=data_file, inits= init_dict,output_dir='/home/apoorva/Desktop/work/MZ_analysis/output',fixed_param=True ,show_console=True)
fit = model.sample(data=data_file, output_dir='/home/apoorva/Desktop/work/MZ_analysis/new_output',show_console=True, inits= init_dict, iter_sampling= 1000, iter_warmup= 100, chains=4, threads_per_chain=5)

# Print summary of the fit
print(fit.summary())

# Diagnose the fit
diagnostic_output = fit.diagnose()
print(diagnostic_output)

#save output to csv file

# fit.save_csvfiles('/home/apoorva/Desktop/work/MZ_analysis/data/output_logistic_growth1')


with open('/home/apoorva/Desktop/work/MZ_analysis/data/New output (pkl)/Ontogeny Output/DDM/T2 Precursor/DDM_T2_LogisticGrowth_varyinginflux_check.pkl', 'wb') as f:
    pkl.dump(fit, f)
    
# print(fit.stan_variables())

# Use diagnose method to check for any potential problems with the fit
print(fit.diagnose())