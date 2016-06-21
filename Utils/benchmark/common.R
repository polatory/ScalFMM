get_data_subset <- function(f, n, h, p, b)
{
    d <- read.csv(file=f,comment.char = "#", sep=",", quote = "\"", head=TRUE,
                  dec=".", colClasses=
	    c("factor",	# model
		"factor",	# algorithm
		"integer",	# nnode
		"integer",	# nthread
		"integer",	# npart
		"integer",	# height
		"integer",	# bsize
		"numeric",	# global_time
		"numeric",	# runtime_time
		"numeric",	# task_time
		"numeric",	# idle_time
		"numeric",	# scheduling_time
		"numeric",	# communuication_time
		"numeric",	# communuication_vol
		"numeric"	# rmen
		))

    d$npart     <- ordered(d$npart)
    d$height    <- ordered(d$height)
    d$bsize     <- ordered(d$bsize)

    if (n)
        d <- subset(d, npart == n)
    if (h)
        d <- subset(d, height == h)
    if (b)
        d <- subset(d, bsize == b | algo == "simple-mpi")

    return (d)
}

get_breaks_runtime <- function()
{
    return (c('implicit', 'explicit', 'implicit limited', 'simple-mpi', 'implicit-two-by'))
}

get_labels_runtime <- function()
{
    return (c('Implicit', 'Explicit', 'Implicit Limited', 'Simple MPI', 'Implicit two by'))
}

get_colors_runtime <- function()
{
    return (c('implicit'  = "#266d83",
              'explicit'   = "#e20025",
			  'implicit limited' = "#bd02b6",
			  'simple-mpi' = "#9aff4f",
			  'implicit-two-by' = "#50de11"))
}

get_bsize_reference <- function()
{
	return (2000)
}
get_one_node_reference_algorithm <- function()
{
	return ("starpu")
}
get_output_directory <- function()
{
	return ("output")
}
# Timings
get_breaks_timings <- function()
{
    return (c('runtime_time', 'communication_time', 'scheduling_time', 'pipeline_time', 'task_time', 'parallel_time'))
}

get_labels_timings <- function()
{
    return (c('Runtime', 'Communication', 'Scheduling', 'Pipeline', 'Task', 'Parallel'))
}

get_colors_timings <- function()
{
    return (c('task_time'       = "#619dff",
              'runtime_time'    = "#01ba38",
              'pipeline_time'   = "#f9766e",
              'scheduling_time' = "#fdbc43",
		      'communication_time' = "#9D2FEA",
              'parallel_time'   = "#000000"))
}
npart_labeller <- function(value)
{
    return (paste("N =", format(as.numeric(value), scientific=TRUE)))
}

height_labeller <- function(value)
{
    return (paste("h =", value))
}
group_size_labeller <- function(value)
{
    return (paste("bs =", value))
}
nnode_labeller <- function(value)
{
    return (paste("node =", value))
}
