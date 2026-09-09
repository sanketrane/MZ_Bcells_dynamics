import os
from cmdstanpy import CmdStanModel
import pickle as pkl
# Specify Stan program file
stan_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/stan_models/new_models/Ontogeny/Spline Fits /T2_emp_func.stan')

# instantiate the model object
model = CmdStanModel(stan_file=stan_file)

# inspect model object
print(model)

# inspect compiled model
print(model.exe_info())

#specify data file

data_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/data/New_data/data_T2fits.json')
# init_dict = {'theta_0': 13.5, 'nu': 0.5, 'sig': 0.5}

#model fit with fixed parameters

#fit = model.sample(data=data_file, inits= init_dict,output_dir='/home/apoorva/Desktop/work/MZ_analysis/output',fixed_param=True ,show_console=True)
fit = model.sample(data=data_file, output_dir='/home/apoorva/Desktop/work/MZ_analysis/new_output',show_console=True, iter_sampling= 2000, iter_warmup= 100, chains=4)

# Print summary of the fit
print(fit.summary())

# Diagnose the fit
diagnostic_output = fit.diagnose()
print(diagnostic_output)

#save output to csv file

# fit.save_csvfiles('/home/apoorva/Desktop/work/MZ_analysis/data/fit_T1.csv')


with open('/home/apoorva/Desktop/work/MZ_analysis/data/emperical_fit_T2.pkl', 'wb') as f:
    pkl.dump(fit, f)
    
# print(fit.stan_variables())

# Use diagnose method to check for any potential problems with the fit
print(fit.diagnose())