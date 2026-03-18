#!/bin/bash
cd /rds/general/user/bs218/home || {
    echo "Error: Cannot access /rds/general/user/bs218/home"
    exit 1
}

# Create logs directory if it doesn't exist
mkdir -p logs

# Set log file with timestamp
LOGFILE="logs/parameter_sweep_SU_N_non_uniform_shifts_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to log file
exec 1> >(tee -a "$LOGFILE")
exec 2>&1

echo "Starting parameter sweep at $(date)"
echo "Log file: $LOGFILE"
echo "-----------------------------------"

job_count=0
username="bs218"  # Set your username here
JOBS_PER_PARAM_SET=1  # Default value, will be updated dynamically
MAX_JOBS=48         # Maximum jobs allowed
MIN_JOBS=40          # Minimum job threshold

# Function to count the number of comma-separated entries in a string
count_entries() {
    echo "$1" | awk -F, '{print NF}'
}

# Function to get an entry from a comma-separated list by position (1-based index)
get_entry() {
    echo "$1" | cut -d, -f"$2"
}

# Function to submit a job with given parameters
submit_job() {
    # Required parameters
    u=$1
    mu_list=$2
    config_folder=$3
    beta_list=$4      # Comma-separated list of beta values

    # Optional parameters with defaults
    mcmc_mode=$5
    [ -z "$mcmc_mode" ] && mcmc_mode=0

    density_mode=$6
    [ -z "$density_mode" ] && density_mode=false

    lnz_mode=$7
    [ -z "$lnz_mode" ] && lnz_mode=false

    seed_offset=$8
    [ -z "$seed_offset" ] && seed_offset=2222092514329988

    # Change alpha_shift to alpha_shift_vec (now a vector)
    alpha_shift_vec=$9
    [ -z "$alpha_shift_vec" ] && alpha_shift_vec="1.90,2.10,2.00"

    # shift to access params beyond 9
    shift 9

    # Add these parameters
    start_order=$1
    [ -z "$start_order" ] && start_order=1

    end_order=$2
    [ -z "$end_order" ] && end_order=7

    cores_cpp=$3
    [ -z "$cores_cpp" ] && cores_cpp=8

    # Add real_type parameter
    real_type=$4
    [ -z "$real_type" ] && real_type="double"

    # Add RG_MCMC parameter
    rg_mcmc=$5
    [ -z "$rg_mcmc" ] && rg_mcmc=0

    # Add N_MIN_ORDER and N_MAX_ORDER parameters
    n_min_order=$6
    [ -z "$n_min_order" ] && n_min_order=3

    n_max_order=$7
    [ -z "$n_max_order" ] && n_max_order=4

    # Fixed parameters from simpler params.txt
    n1_list="32,24,18,10,7,4,2"
    n2_list="32,24,18,10,7,4,2"
    cheby_degree=36
    greens_func_exp_cut=32.0
    mcmc_cut="1e-14"
    eps="1e-14"
    memory_cpp_in_gb=32
    seed_interval=100000

    # Calculate JOBS_PER_PARAM_SET based on start_order and end_order
    JOBS_PER_PARAM_SET=$(expr $end_order - $start_order + 1)
    # Update the global variable for use in the main loop
    export JOBS_PER_PARAM_SET

    u_value="$u"

    # mu_h=$(echo "$mu + 0.0001" | bc -l | awk '{printf "%.6f", $0}')
    # mu_l=$(echo "$mu - 0.0001" | bc -l | awk '{printf "%.6f", $0}')
    # Process mu list to calculate mu_h and mu_l
    mu_h_list=""
    mu_l_list=""
    mu_count=$(count_entries "$mu_list")
    mu_idx=1

    while [ $mu_idx -le $mu_count ]; do
        mu_val=$(get_entry "$mu_list" $mu_idx)
        # Calculate high and low values for this beta
        mu_h_val=$(echo "$mu_val + 0.0002" | bc -l)
        mu_l_val=$(echo "$mu_val - 0.0002" | bc -l)

        # Append to our comma-separated lists
        if [ -z "$mu_h_list" ]; then
            mu_h_list="$mu_h_val"
            mu_l_list="$mu_l_val"
        else
            mu_h_list="$mu_h_list,$mu_h_val"
            mu_l_list="$mu_l_list,$mu_l_val"
        fi

        mu_idx=$(expr $mu_idx + 1)
    done

    # Create a mode signature for folder name
    mode_sig=""
    [ "$density_mode" = "true" ] && mode_sig="${mode_sig}D"
    [ "$lnz_mode" = "true" ] && mode_sig="${mode_sig}L"
    [ -z "$mode_sig" ] && mode_sig="none"

    # Add a timestamp to ensure uniqueness
    timestamp=$(date +%Y%m%d%H%M%S)

    # Get the first beta value for folder name
    first_beta=$(get_entry "$beta_list" 1)

    # Create folder name based on parameters
    folder_name="U${u}_mu${mu_list}_mcmc${mcmc_mode}_beta${first_beta}_modes${mode_sig}_rg${rg_mcmc}_${timestamp}"

    echo "Setting MU=$mu_list, U=$u_value, MCMC_mode=$mcmc_mode, BETA list=$beta_list"
    echo "Setting modes: DENSITY=$density_mode, LNZ=$lnz_mode"
    echo "Setting SEED_OFFSET=$seed_offset"
    echo "Setting ALPHA_SHIFT_VEC=$alpha_shift_vec"
    echo "Setting START_ORDER=$start_order, END_ORDER=$end_order, CORES_CPP=$cores_cpp"
    echo "Setting REAL_TYPE=$real_type"
    echo "Setting RG_MCMC=$rg_mcmc"
    echo "Setting N_MIN_ORDER=$n_min_order, N_MAX_ORDER=$n_max_order"
    echo "This configuration will create $JOBS_PER_PARAM_SET jobs"

    # Ensure the config folder exists
    if [ ! -d "$config_folder" ]; then
        echo "Error: Config folder $config_folder does not exist"
        return 1
    fi

    # Enter the config folder
    cd "$config_folder" || { echo "Failed to enter $config_folder"; return 1; }

    # Create folder for this parameter combination
    mkdir -p "$folder_name"

    # Copy params.txt and submit.pbs to the new folder
    cp params.txt "$folder_name/"
    cp submit.pbs "$folder_name/"

    # Update params.txt with new values in the new folder - only parameters that exist
    sed -i "s/^N1=.*$/N1=$n1_list/g" "$folder_name/params.txt"
    sed -i "s/^N2=.*$/N2=$n2_list/g" "$folder_name/params.txt"
    sed -i "s/^CHEBY_DEGREE=.*$/CHEBY_DEGREE=$cheby_degree/g" "$folder_name/params.txt"
    sed -i "s/^GREENS_FUNC_EXP_CUT=.*$/GREENS_FUNC_EXP_CUT=$greens_func_exp_cut/g" "$folder_name/params.txt"
    sed -i "s/^MCMC_CUT=.*$/MCMC_CUT=$mcmc_cut/g" "$folder_name/params.txt"
    sed -i "s/^EPS=.*$/EPS=$eps/g" "$folder_name/params.txt"
    sed -i "s/^MU_h=.*$/MU_h=$mu_h_list/g" "$folder_name/params.txt"
    sed -i "s/^MU_l=.*$/MU_l=$mu_l_list/g" "$folder_name/params.txt"
    sed -i "s/^MU=.*$/MU=$mu_list/g" "$folder_name/params.txt"
    sed -i "s/^REAL_TYPE=.*$/REAL_TYPE=$real_type/g" "$folder_name/params.txt"
    sed -i "s/^MCMC_mode=.*$/MCMC_mode=$mcmc_mode/g" "$folder_name/params.txt"
    sed -i "s/^CORES_CPP=.*$/CORES_CPP=$cores_cpp/g" "$folder_name/params.txt"
    sed -i "s/^U=.*$/U=$u_value/g" "$folder_name/params.txt"
    sed -i "s/^BETA=.*$/BETA=$beta_list/g" "$folder_name/params.txt"
    sed -i "s/^ALPHA_SHIFT_VEC=.*$/ALPHA_SHIFT_VEC=$alpha_shift_vec/g" "$folder_name/params.txt"
    sed -i "s/^DENSITY_MODE=.*$/DENSITY_MODE=$density_mode/g" "$folder_name/params.txt"
    sed -i "s/^LNZ_MODE=.*$/LNZ_MODE=$lnz_mode/g" "$folder_name/params.txt"
    sed -i "s/^RG_MCMC=.*$/RG_MCMC=$rg_mcmc/g" "$folder_name/params.txt"
    sed -i "s/^MEMORY_CPP_IN_GB=.*$/MEMORY_CPP_IN_GB=$memory_cpp_in_gb/g" "$folder_name/params.txt"
    sed -i "s/^SEED_OFFSET=.*$/SEED_OFFSET=$seed_offset/g" "$folder_name/params.txt"
    sed -i "s/^SEED_INTERVAL=.*$/SEED_INTERVAL=$seed_interval/g" "$folder_name/params.txt"
    sed -i "s/^START_ORDER=.*$/START_ORDER=$start_order/g" "$folder_name/params.txt"
    sed -i "s/^END_ORDER=.*$/END_ORDER=$end_order/g" "$folder_name/params.txt"
    sed -i "s/^N_MIN_ORDER=.*$/N_MIN_ORDER=$n_min_order/g" "$folder_name/params.txt"
    sed -i "s/^N_MAX_ORDER=.*$/N_MAX_ORDER=$n_max_order/g" "$folder_name/params.txt"

    # Submit the job from within the folder
    cd "$folder_name" || { echo "Failed to enter $folder_name"; return 1; }
    qsub submit.pbs
    cd ../..  # Return to the original directory (parent of config folder)

    # Increment job counter for the parameter set
    job_count=$(expr $job_count + 1)

    echo "Parameter set #$job_count submitted for U=$u, mu=$mu_list, MCMC_mode=$mcmc_mode in $config_folder/$folder_name"
    echo "Using BETA values: $beta_list"
    echo "This launches $JOBS_PER_PARAM_SET jobs (orders $start_order to $end_order)"

    # Wait between job submissions to avoid flooding the scheduler
    sleep 180
}

# Function to get parameters based on mu value
# Returns: U alpha_shift_vec
get_parameters_for_mu() {
    mu=$1
    # Format: U alpha_shift_vec (as comma-separated vector)
    case $mu in
        "0.0,0.5,1.0") echo "2.3 1.2,1.5,1.8" ;;
    esac
}

# Set up parameters as individual variables
job_param_count=14

# Base parameters
BETA_LIST="3.0,2.2,1.3"
SEED_PREFIX=784

# Seed range constants
MIN_SEED=22192833664439
MAX_SEED=55192833664439
RANGE=$((MAX_SEED - MIN_SEED + 1))

generate_random_seed_suffix() {
    local random_bytes
    local seed

    random_bytes=$(od -An -N8 -tu8 < /dev/urandom | tr -d ' ')
    seed=$(echo "$MIN_SEED + $random_bytes % $RANGE" | bc)

    echo "$seed"
}

# Array of all mu values
MU_VALUES=("0.0,0.5,1.0")

param_num=1

# Generate parameters
echo "Generating parameters ..."
for mu in "${MU_VALUES[@]}"; do
    # Get parameters for this mu
    params=$(get_parameters_for_mu $mu)
    U=$(echo $params | cut -d' ' -f1)
    alpha_shift_vec=$(echo $params | cut -d' ' -f2)

    # First set: density mode
    density_mode="false"
    lnz_mode="true"

    for ii in 0; do
        SEED_SUFFIX=$(generate_random_seed_suffix)
        # Added RG_MCMC (0), N_MIN_ORDER (3), N_MAX_ORDER (4) at the end
        eval "job_param_${param_num}=\"$U $mu config_arrayjobs_SU_N_non_uniform_shifts $BETA_LIST -1 $density_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_shift_vec 6 6 64 double 0 6 7\""
        param_num=$((param_num + 1))
        SEED_PREFIX=$((SEED_PREFIX + 1))
        sleep 0.1
    done
    for ii in 0; do
        SEED_SUFFIX=$(generate_random_seed_suffix)
        # Added RG_MCMC (0), N_MIN_ORDER (3), N_MAX_ORDER (4) at the end
        eval "job_param_${param_num}=\"$U $mu config_arrayjobs_SU_N_non_uniform_shifts $BETA_LIST -1 $density_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_shift_vec 6 6 64 double 0 5 6\""
        param_num=$((param_num + 1))
        SEED_PREFIX=$((SEED_PREFIX + 1))
        sleep 0.1
    done
    for ii in 0; do
        SEED_SUFFIX=$(generate_random_seed_suffix)
        # Added RG_MCMC (0), N_MIN_ORDER (3), N_MAX_ORDER (4) at the end
        eval "job_param_${param_num}=\"$U $mu config_arrayjobs_SU_N_non_uniform_shifts $BETA_LIST -1 $density_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_shift_vec 6 6 64 double 0 4 5\""
        param_num=$((param_num + 1))
        SEED_PREFIX=$((SEED_PREFIX + 1))
        sleep 0.1
    done

    for mode in 0; do
        SEED_SUFFIX=$(generate_random_seed_suffix)
        # Added RG_MCMC (0), N_MIN_ORDER (3), N_MAX_ORDER (4) at the end
        eval "job_param_${param_num}=\"$U $mu config_multicores_SU_N_non_uniform_shifts $BETA_LIST $mode $density_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_shift_vec 1 1 8 double 1 6 7\""
        param_num=$((param_num + 1))
        SEED_PREFIX=$((SEED_PREFIX + 1))
        sleep 0.1
    done
    for mode in 0; do
        SEED_SUFFIX=$(generate_random_seed_suffix)
        # Added RG_MCMC (0), N_MIN_ORDER (3), N_MAX_ORDER (4) at the end
        eval "job_param_${param_num}=\"$U $mu config_multicores_SU_N_non_uniform_shifts $BETA_LIST $mode $density_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_shift_vec 2 4 64 double 1 6 7\""
        param_num=$((param_num + 1))
        SEED_PREFIX=$((SEED_PREFIX + 1))
        sleep 0.1
    done
    for mode in 1 2 3 4 5 6 7 8 9; do
        SEED_SUFFIX=$(generate_random_seed_suffix)
        # Added RG_MCMC (0), N_MIN_ORDER (3), N_MAX_ORDER (4) at the end
        eval "job_param_${param_num}=\"$U $mu config_multicores_SU_N_non_uniform_shifts $BETA_LIST $mode $density_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_shift_vec 1 1 64 double 1 6 7\""
        param_num=$((param_num + 1))
        SEED_PREFIX=$((SEED_PREFIX + 1))
        sleep 0.1
    done
done

# Function to count current jobs for the user
count_jobs() {
    qselect -u "$username" | wc -l
}

# Function to submit a job from the parameter sets
submit_next_job() {
    # Calculate which parameter set to use (job_count + 1)
    param_idx=$(expr $job_count + 1)

    # Get the corresponding parameter set using eval
    params=$(eval echo \${job_param_$param_idx})

    # Submit the job with those parameters
    # shellcheck disable=SC2086
    submit_job $params
}

# Extract the start_order and end_order from a parameter set to calculate jobs
get_jobs_in_param_set() {
    param_idx=$1
    params=$(eval echo \${job_param_$param_idx})

    # Extract 10th parameter (start_order) and 11th parameter (end_order)
    # Split the parameters by space
    param_array=($params)

    # Parameters are 0-indexed in the array
    start_order=${param_array[9]}
    end_order=${param_array[10]}

    # Calculate jobs in this parameter set
    expr $end_order - $start_order + 1
}

# For initial display, calculate the jobs from the first parameter set
FIRST_PARAM_JOBS=$(get_jobs_in_param_set 1)

# Main loop to monitor and maintain job count
echo "Starting job monitor for user $username."
echo "Target: minimum $MIN_JOBS jobs, maximum $MAX_JOBS jobs."
echo "Jobs per parameter set varies based on order range (first set: $FIRST_PARAM_JOBS jobs)"
echo "Total parameter sets to submit: $job_param_count"
echo "Mode configuration: Running density and lnZ modes for all mu values"
echo "Seed pattern: Starting with prefix $SEED_PREFIX, random suffix in range $MIN_SEED-$MAX_SEED"

while true; do
    current_jobs=$(count_jobs)
    echo "$(date): Current job count: $current_jobs (Target: $MIN_JOBS-$MAX_JOBS)"

    # Check if all parameter sets have been submitted
    if [ $job_count -ge $job_param_count ]; then
        echo "ALL PARAMETER SETS HAVE BEEN SUBMITTED ($job_count of $job_param_count)"
        echo "Exiting job monitor."
        break  # Exit the while loop once all parameter sets are submitted
    fi

    # Get the number of jobs that would be added by the next parameter set
    next_param_idx=$(expr $job_count + 1)
    next_jobs=$(get_jobs_in_param_set $next_param_idx)

    # Submit new jobs as long as adding the jobs for the next parameter set won't exceed MAX_JOBS
    # and we still have parameter sets to submit
    while [ "$(expr $current_jobs + $next_jobs)" -lt "$MAX_JOBS" ] && [ $job_count -lt $job_param_count ]; do
        echo "Submitting parameter set #$(expr $job_count + 1) (adding $next_jobs jobs). Current jobs: $current_jobs, Maximum: $MAX_JOBS"

        submit_next_job

        # Update current job count after submission
        current_jobs=$(count_jobs)

        # Update next_param_idx and next_jobs for the next iteration
        next_param_idx=$(expr $job_count + 1)
        if [ $next_param_idx -le $job_param_count ]; then
            next_jobs=$(get_jobs_in_param_set $next_param_idx)
        else
            next_jobs=0
        fi

        echo "Parameter set submitted. New job count: $current_jobs"
    done

    echo "Next check in 3 minutes ($(date -d '+3 minutes'))"
    sleep 180
done

echo "Job monitor completed successfully."