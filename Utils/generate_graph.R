source("common.R")
source("gen_times_taskdep.R")
source("gen_efficiencies_taskdep.R")
source("gen_speedup_plots.R")
###
# Generate display of bars with the time spent in Task, Runtime and Idle.
###
gen_times_taskdep("loutre.db")
gen_efficiencies("loutre.db")
gen_speedup("loutre.db")
