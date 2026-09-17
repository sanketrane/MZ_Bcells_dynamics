import os
import subprocess
from cmdstanpy import CmdStanModel
import pickle as pkl


def needs_recompile(exe_path: str) -> bool:
    """Recompile when executable is missing, not executable, or wrong platform binary."""
    if not os.path.exists(exe_path):
        return True
    if not os.access(exe_path, os.X_OK):
        return True
    try:
        file_info = subprocess.run(
            ["file", exe_path],
            capture_output=True,
            text=True,
            check=False,
        ).stdout
        # macOS needs Mach-O binaries; an ELF binary indicates a stale Linux build.
        return "ELF" in file_info
    except Exception:
        # If detection fails, avoid forcing rebuild and let CmdStanPy handle normal flow.
        return False

# Specify Stan program file
stan_file = os.path.join('/Users/apoorvasingh/Desktop/MZB-new-analysis/stan_models/Ontogeny/DDM Influx/T2 Precursor/DDdivision_independent.stan')
# instantiate the model object
model = CmdStanModel(stan_file=stan_file)

if needs_recompile(model.exe_file):
    print(f"Recompiling Stan model for this machine: {model.exe_file}")
    model.compile(force=True)

# inspect model object
print(model)

# inspect compiled model
print(model.exe_info())

#specify data file

data_file = os.path.join('/Users/apoorvasingh/Desktop/MZB-new-analysis/data/New_data/MZcounts_ontogeny.json')
# init_dict = {'beta': 3.5, 'y0_Log': 14, 'kappa_0': 0.5, 'rho': 0.01, 'delta': 0.09, 'sigma': [1,1,1,1]}

#model fit with fixed parameters

#fit = model.sample(data=data_file, inits= init_dict,output_dir='/Users/apoorvasingh/Desktop/MZB-new-analysis/output',fixed_param=True ,show_console=True)
fit = model.sample(data=data_file, output_dir='/Users/apoorvasingh/Desktop/MZB-new-analysis/new_output',show_console=True, iter_sampling= 2000, iter_warmup= 100, chains=4, threads_per_chain=12)

# Print summary of the fit
print(fit.summary())

# Diagnose the fit
diagnostic_output = fit.diagnose()
print(diagnostic_output)

with open('/Users/apoorvasingh/Desktop/MZB-new-analysis/data/New output (pkl)/Ontogeny Output/DDM/T2 Precursor/DDM_T2_Division_independent.pkl', 'wb') as f:
    pkl.dump(fit, f)
    
# print(fit.stan_variables())

# Use diagnose method to check for any potential problems with the fit
print(fit.diagnose())