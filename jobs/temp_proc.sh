#!/bin/bash
#SBATCH --job-name=proc_211_lkg_tnl_U           # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=heikki.tolppanen@oulu.fi    # Where to send mail	
#SBATCH --nodes=1                               # Run all processes on a single node	
#SBATCH --ntasks=20                             # Run a single task		
#SBATCH --cpus-per-task=1                       # Number of CPU cores per task
#SBATCH --mem=40gb                              # Job memory request
#SBATCH --time=24:00:00                         # Time limit hrs:min:sec
#SBATCH --output=jobs/logs/%x_%j_log.txt        # Standard output and error log

file_to_run="LRUFUP_Vault/Data/0PAPER/LeakageTunneling/InSubspace/OffResonance211/AnhHeatmap/code.jl"

echo "$file_to_run"; hostname; date; echo

start_time=$(date +%s)  # Record start time
julia -p 20 -e "
using Distributed
@everywhere begin
    using Pkg
    Pkg.activate(\".\")
end
include(\"$file_to_run\")
"
end_time=$(date +%s)    # Record end time
elapsed_time=$((end_time - start_time))  # Calculate elapsed time

echo
./jobs/_elapsed_time.sh $elapsed_time
date
    