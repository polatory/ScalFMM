source("common.R")
source("gen_times_taskdep.R")
source("gen_efficiencies_taskdep.R")
source("gen_speedup_plots.R")
source("gen_parallel_efficiency_plots.R")
source("gen_normalized_time_plots.R")
source("gen_comm.R")
source("gen_gantt.R")

###
# Generate display of bars with the time spent in Task, Runtime and Idle.
###
gen_times_taskdep("loutre.db")
gen_efficiencies("loutre.db")
gen_comm("loutre.db")
gen_speedup("loutre.db")
gen_pareff("loutre.db")
gen_normalized_time("loutre.db")
gen_gantt("canard.db")
