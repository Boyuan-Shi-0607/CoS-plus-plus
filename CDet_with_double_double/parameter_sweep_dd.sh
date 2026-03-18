#!/bin/bash
cd /rds/general/user/bs218/home || {
    echo "Error: Cannot access /rds/general/user/bs218/home"
    exit 1
}

# Create logs directory if it doesn't exist
mkdir -p logs

# Set log file with timestamp
LOGFILE="logs/parameter_sweep_dd_$(date +%Y%m%d_%H%M%S).log"

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
    mu=$2
    config_folder=$3  # Either "config_bipar_arrayjobs" or "config_bipar_multicores"
    beta_list=$4      # Comma-separated list of beta values

    # Optional parameters with defaults
    mcmc_mode=$5
    [ -z "$mcmc_mode" ] && mcmc_mode=0

    density_mode=$6
    [ -z "$density_mode" ] && density_mode=false

    energy_mode=$7
    [ -z "$energy_mode" ] && energy_mode=false

    double_occupancy_mode=$8
    [ -z "$double_occupancy_mode" ] && double_occupancy_mode=true

    compressibility_mode=$9
    [ -z "$compressibility_mode" ] && compressibility_mode=false

    # shift to access params beyond 9
    shift 9
    staggered_magnetization_mode=$1
    [ -z "$staggered_magnetization_mode" ] && staggered_magnetization_mode=false

    lnz_mode=$2
    [ -z "$lnz_mode" ] && lnz_mode=false

    seed_offset=$3
    [ -z "$seed_offset" ] && seed_offset=214409251432

    # Add these four parameters
    alpha_up=$4
    [ -z "$alpha_up" ] && alpha_up=1.97

    alpha_down=$5
    [ -z "$alpha_down" ] && alpha_down=1.97

    delta_up=$6
    [ -z "$delta_up" ] && delta_up=-0.1

    delta_down=$7
    [ -z "$delta_down" ] && delta_down=0.1

    # Add these three new parameters
    start_order=$8
    [ -z "$start_order" ] && start_order=1

    end_order=$9
    [ -z "$end_order" ] && end_order=1

    shift 9
    cores_cpp=$1
    [ -z "$cores_cpp" ] && cores_cpp=100

    # Add real_type parameter
    real_type=$2
    [ -z "$real_type" ] && real_type="long double"

    # Translate "long" to "long double"
    if [ "$real_type" = "long" ]; then
        real_type="long double"
    fi

    cheby_degree=35
    greens_func_exp_cut=61.0
    mcmc_cut="1e-28"
    eps="1e-28"
    tol_fpm="1e-26"

    # Calculate JOBS_PER_PARAM_SET based on start_order and end_order
    JOBS_PER_PARAM_SET=$(expr $end_order - $start_order + 1)
    # Update the global variable for use in the main loop
    export JOBS_PER_PARAM_SET

    times_cpp_order9="08"
    times_cpp_order8="08"

    u_value="$u"

    mu_h=$(echo "$mu + 0.001" | bc -l | awk '{printf "%.6f", $0}')
    mu_l=$(echo "$mu - 0.001" | bc -l | awk '{printf "%.6f", $0}')
    u_h="$(echo "$u+0.0001" | bc)"
    u_l="$(echo "$u-0.0001" | bc)"
    H_h="0.001"
    H_l="-0.001"

    # Process beta list to calculate beta_h and beta_l
    beta_h_list=""
    beta_l_list=""
    beta_count=$(count_entries "$beta_list")
    beta_idx=1

    while [ $beta_idx -le $beta_count ]; do
        beta_val=$(get_entry "$beta_list" $beta_idx)
        # Calculate high and low values for this beta
        beta_h_val=$(echo "$beta_val + 0.0002" | bc -l)
        beta_l_val=$(echo "$beta_val - 0.0002" | bc -l)

        # Append to our comma-separated lists
        if [ -z "$beta_h_list" ]; then
            beta_h_list="$beta_h_val"
            beta_l_list="$beta_l_val"
        else
            beta_h_list="$beta_h_list,$beta_h_val"
            beta_l_list="$beta_l_list,$beta_l_val"
        fi

        beta_idx=$(expr $beta_idx + 1)
    done

    # Create a mode signature for folder name
    mode_sig=""
    [ "$density_mode" = "true" ] && mode_sig="${mode_sig}D"
    [ "$energy_mode" = "true" ] && mode_sig="${mode_sig}E"
    [ "$double_occupancy_mode" = "true" ] && mode_sig="${mode_sig}O"
    [ "$compressibility_mode" = "true" ] && mode_sig="${mode_sig}C"
    [ "$staggered_magnetization_mode" = "true" ] && mode_sig="${mode_sig}S"
    [ "$lnz_mode" = "true" ] && mode_sig="${mode_sig}L"
    [ -z "$mode_sig" ] && mode_sig="none"

    # Add a timestamp to ensure uniqueness
    timestamp=$(date +%Y%m%d%H%M%S)

    # Get the first beta value for folder name
    first_beta=$(get_entry "$beta_list" 1)

    # Create folder name based on parameters
    folder_name="U${u}_mu${mu}_mcmc${mcmc_mode}_beta${first_beta}_modes${mode_sig}_dd_${timestamp}"

    echo "Setting MU=$mu, U=$u_value, MCMC_mode=$mcmc_mode, BETA list=$beta_list"
    echo "Setting modes: DENSITY=$density_mode, ENERGY=$energy_mode, DOUBLE_OCCUPANCY=$double_occupancy_mode"
    echo "Setting modes: COMPRESSIBILITY=$compressibility_mode, STAGGERED_MAGNETIZATION=$staggered_magnetization_mode, LNZ=$lnz_mode"
    echo "Setting SEED_OFFSET=$seed_offset"
    echo "Setting ALPHA_UP=$alpha_up, ALPHA_DOWN=$alpha_down, DELTA_UP=$delta_up, DELTA_DOWN=$delta_down"
    echo "Setting START_ORDER=$start_order, END_ORDER=$end_order, CORES_CPP=$cores_cpp"
    echo "Setting REAL_TYPE=$real_type (CHEBY_DEGREE=$cheby_degree, GREENS_FUNC_EXP_CUT=$greens_func_exp_cut, MCMC_CUT=$mcmc_cut, EPS=$eps, TOL_FPM=$tol_fpm)"
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

    # Update params.txt with new values in the new folder - improved with anchored patterns
    sed -i "s/^MU=.*$/MU=$mu/g" "$folder_name/params.txt"
    sed -i "s/^MU_h=.*$/MU_h=$mu_h/g" "$folder_name/params.txt"
    sed -i "s/^MU_l=.*$/MU_l=$mu_l/g" "$folder_name/params.txt"
    sed -i "s/^U=.*$/U=$u_value/g" "$folder_name/params.txt"
    sed -i "s/^U_h=.*$/U_h=$u_h/g" "$folder_name/params.txt"
    sed -i "s/^U_l=.*$/U_l=$u_l/g" "$folder_name/params.txt"
    sed -i "s/^BETA=.*$/BETA=$beta_list/g" "$folder_name/params.txt"
    sed -i "s/^BETA_h=.*$/BETA_h=$beta_h_list/g" "$folder_name/params.txt"
    sed -i "s/^BETA_l=.*$/BETA_l=$beta_l_list/g" "$folder_name/params.txt"
    sed -i "s/^MCMC_mode=.*$/MCMC_mode=$mcmc_mode/g" "$folder_name/params.txt"
    sed -i "s/^DENSITY_MODE=.*$/DENSITY_MODE=$density_mode/g" "$folder_name/params.txt"
    sed -i "s/^ENERGY_MODE=.*$/ENERGY_MODE=$energy_mode/g" "$folder_name/params.txt"
    sed -i "s/^DOUBLE_OCCUPANCY_MODE=.*$/DOUBLE_OCCUPANCY_MODE=$double_occupancy_mode/g" "$folder_name/params.txt"
    sed -i "s/^COMPRESSIBILITY_MODE=.*$/COMPRESSIBILITY_MODE=$compressibility_mode/g" "$folder_name/params.txt"
    sed -i "s/^STAGGERED_MAGNETIZATION_MODE=.*$/STAGGERED_MAGNETIZATION_MODE=$staggered_magnetization_mode/g" "$folder_name/params.txt"
    sed -i "s/^LNZ_MODE=.*$/LNZ_MODE=$lnz_mode/g" "$folder_name/params.txt"
    sed -i "s/^SEED_OFFSET=.*$/SEED_OFFSET=$seed_offset/g" "$folder_name/params.txt"
    sed -i "s/^ALPHA_UP=.*$/ALPHA_UP=$alpha_up/g" "$folder_name/params.txt"
    sed -i "s/^ALPHA_DOWN=.*$/ALPHA_DOWN=$alpha_down/g" "$folder_name/params.txt"
    sed -i "s/^DELTA_UP=.*$/DELTA_UP=$delta_up/g" "$folder_name/params.txt"
    sed -i "s/^DELTA_DOWN=.*$/DELTA_DOWN=$delta_down/g" "$folder_name/params.txt"
    sed -i "s/^H_h=.*$/H_h=$H_h/g" "$folder_name/params.txt"
    sed -i "s/^H_l=.*$/H_l=$H_l/g" "$folder_name/params.txt"
    sed -i "s/^TIMES_CPP_IN_H_ORDER_9=.*$/TIMES_CPP_IN_H_ORDER_9=$times_cpp_order9/g" "$folder_name/params.txt"
    sed -i "s/^TIMES_CPP_IN_H_ORDER_8=.*$/TIMES_CPP_IN_H_ORDER_8=$times_cpp_order8/g" "$folder_name/params.txt"
    sed -i "s/^START_ORDER=.*$/START_ORDER=$start_order/g" "$folder_name/params.txt"
    sed -i "s/^END_ORDER=.*$/END_ORDER=$end_order/g" "$folder_name/params.txt"
    sed -i "s/^CORES_CPP=.*$/CORES_CPP=$cores_cpp/g" "$folder_name/params.txt"
    sed -i "s/^REAL_TYPE=.*$/REAL_TYPE=$real_type/g" "$folder_name/params.txt"
    sed -i "s/^CHEBY_DEGREE=.*$/CHEBY_DEGREE=$cheby_degree/g" "$folder_name/params.txt"
    sed -i "s/^GREENS_FUNC_EXP_CUT=.*$/GREENS_FUNC_EXP_CUT=$greens_func_exp_cut/g" "$folder_name/params.txt"
    sed -i "s/^MCMC_CUT=.*$/MCMC_CUT=$mcmc_cut/g" "$folder_name/params.txt"
    sed -i "s/^EPS=.*$/EPS=$eps/g" "$folder_name/params.txt"
    sed -i "s/^TOL_FPM=.*$/TOL_FPM=$tol_fpm/g" "$folder_name/params.txt"

    # Submit the job from within the folder
    cd "$folder_name" || { echo "Failed to enter $folder_name"; return 1; }
    qsub submit.pbs
    cd ../..  # Return to the original directory (parent of config folder)

    # Increment job counter for the parameter set
    job_count=$(expr $job_count + 1)

    echo "Parameter set #$job_count submitted for U=$u, μ=$mu, MCMC_mode=$mcmc_mode in $config_folder/$folder_name"
    echo "Using BETA values: $beta_list"
    echo "This launches $JOBS_PER_PARAM_SET jobs (orders $start_order to $end_order)"

    # Wait between job submissions to avoid flooding the scheduler
    sleep 180
}

# Function to get parameters based on mu and beta_list combination
# Returns: U alpha_up alpha_down delta_up delta_down
get_parameters_for_mu_beta() {
    mu=$1
    beta_list=$2
    # Format: U alpha_up alpha_down delta_up delta_down
    case "${mu}_${beta_list}" in
        "1.8_10.0,4.0") echo "6.0 2.58 2.58 -0.195 0.195" ;;
        "2.2_3.5,2.0") echo "5.0 2.5 2.5 0.0 0.0" ;;
        # Add more mu_beta combinations as needed
    esac
}

# Arrays of paired mu and beta_list values (treated as single parameter sets)
MU_BETA_PAIRS=(
    "1.8:10.0,4.0"
)

# Set up parameters as individual variables
job_param_count=1
SEED_PREFIX=366

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

param_num=1

# Process each mu_beta pair as a single parameter set
echo "Generating parameters ..."
for pair in "${MU_BETA_PAIRS[@]}"; do
    # Split the pair into mu and beta_list
    mu=$(echo $pair | cut -d':' -f1)
    beta_list=$(echo $pair | cut -d':' -f2)
    
    # Get parameters for this mu and beta_list combination
    params=$(get_parameters_for_mu_beta $mu "$beta_list")
    
    # Check if parameters were found
    if [ -z "$params" ]; then
        echo "Warning: No parameters defined for mu=$mu, beta_list=$beta_list - skipping"
        continue
    fi
    
    U=$(echo $params | cut -d' ' -f1)
    alpha_up=$(echo $params | cut -d' ' -f2)
    alpha_down=$(echo $params | cut -d' ' -f3)
    delta_up=$(echo $params | cut -d' ' -f4)
    delta_down=$(echo $params | cut -d' ' -f5)
    
    echo "Processing parameter set: mu=$mu, beta_list=$beta_list with U=$U"

    # Set mode flags
    density_mode="false"
    energy_mode="false"
    double_occupancy_mode="false"
    compressibility_mode="false"
    staggered_magnetization_mode="false"
    lnz_mode="true"

    # for ii in 0 1 2 3; do
        # SEED_SUFFIX=$(generate_random_seed_suffix)
        # eval "job_param_${param_num}=\"$U $mu config_bipar_multicores_dd $beta_list 0 $density_mode $energy_mode $double_occupancy_mode $compressibility_mode $staggered_magnetization_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_up $alpha_down $delta_up $delta_down 6 7 128 double\""
        # param_num=$((param_num + 1))
        # SEED_PREFIX=$((SEED_PREFIX + 1))
        # sleep 0.1
    # done
    # for ii in 0 1 2; do
        # SEED_SUFFIX=$(generate_random_seed_suffix)
        # eval "job_param_${param_num}=\"$U $mu config_bipar_multicores_dd $beta_list 0 $density_mode $energy_mode $double_occupancy_mode $compressibility_mode $staggered_magnetization_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_up $alpha_down $delta_up $delta_down 2 5 64 double\""
        # param_num=$((param_num + 1))
        # SEED_PREFIX=$((SEED_PREFIX + 1))
        # sleep 0.1
    # done
    # for ii in 0; do
        # SEED_SUFFIX=$(generate_random_seed_suffix)
        # eval "job_param_${param_num}=\"$U $mu config_bipar_multicores_dd $beta_list 0 $density_mode $energy_mode $double_occupancy_mode $compressibility_mode $staggered_magnetization_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_up $alpha_down $delta_up $delta_down 1 1 4 double\""
        # param_num=$((param_num + 1))
        # SEED_PREFIX=$((SEED_PREFIX + 1))
        # sleep 0.1
    # done
    # for mode in 1 2 3 4 5 6 7; do
        # SEED_SUFFIX=$(generate_random_seed_suffix)
        # eval "job_param_${param_num}=\"$U $mu config_bipar_multicores_dd $beta_list $mode $density_mode $energy_mode $double_occupancy_mode $compressibility_mode $staggered_magnetization_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_up $alpha_down $delta_up $delta_down 1 2 64 double\""
        # param_num=$((param_num + 1))
        # SEED_PREFIX=$((SEED_PREFIX + 1))
        # sleep 0.1
    # done
    for ii in 0; do
        SEED_SUFFIX=$(generate_random_seed_suffix)
        eval "job_param_${param_num}=\"$U $mu config_bipar_arrayjobs_short8_dd $beta_list 0 $density_mode $energy_mode $double_occupancy_mode $compressibility_mode $staggered_magnetization_mode $lnz_mode ${SEED_PREFIX}${SEED_SUFFIX} $alpha_up $alpha_down $delta_up $delta_down 9 9 2000 double\""
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

    # Extract 17th parameter (start_order) and 18th parameter (end_order)
    # Split the parameters by space
    param_array=($params)

    # Parameters are 0-indexed in the array
    start_order=${param_array[16]}
    end_order=${param_array[17]}

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
echo "Mode configuration: Running all four modes (lnZ, density, staggered_magnetization, double_occupancy) for all mu values"
echo "Seed pattern: Starting with prefix 500, random suffix in range $MIN_SEED-$MAX_SEED"

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
