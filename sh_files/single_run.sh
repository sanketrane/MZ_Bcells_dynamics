
while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

./stan_models/ki67_SHM_${modelname} sample num_warmup=300 num_samples=500 data file=datafiles/MZ_data.Rdump output file=save_csv/shm_${modelname}.csv 
