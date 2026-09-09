import os
import argparse
from cmdstanpy import CmdStanModel
import pickle as pkl

def main(stan_file, data_file, output_dir, model_name):
    # instantiate the model object
    model = CmdStanModel(stan_file=stan_file)

    # inspect model object
    print(model)

    # inspect compiled model
    print(model.exe_info())

    # model fit with fixed parameters
    fit = model.sample(data=data_file, output_dir=output_dir, show_progress=True,
                       show_console=False, iter_sampling=1000, iter_warmup=200, chains=5,
                       threads_per_chain=10)

    # Print summary of the fit
    print(fit.summary())

    # Diagnose the fit
    diagnostic_output = fit.diagnose()
    print(diagnostic_output)

    # Save output to pickle file
    output_file = os.path.join('/home/apoorva/Desktop/work/MZ_analysis/data/New output (pkl)/Steady State/SHM', f'{model_name}.pkl')
    with open(output_file, 'wb') as f:
        pkl.dump(fit, f)

    # Use diagnose method to check for any potential problems with the fit
    print(fit.diagnose())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Stan models with CmdStanPy.')
    parser.add_argument('--stan_file', required=True, help='Path to the Stan model file.')
    parser.add_argument('--data_file', type=str, required=True, help='Path to the data file.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the output.')
    parser.add_argument('--model_name', type=str, required=True, help='Name of the model.')


    args = parser.parse_args()

    main(args.stan_file, args.data_file, args.output_dir, args.model_name)