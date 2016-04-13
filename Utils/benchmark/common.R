get_data_subset <- function(f, n, h, p)
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
		"numeric"	# rmen
		))

    d$npart     <- ordered(d$npart)
    d$height    <- ordered(d$height)
    d$bsize     <- ordered(d$bsize)

    if (n)
        d <- subset(d, npart == n)
    if (h)
        d <- subset(d, height == h)

    return (d)
}
get_db_with_time <- function(f, s, n, h, m)
{
    d <- read.csv(file=f,comment.char = "#", sep=",", quote = "\"", head=TRUE,
                  dec=".", colClasses=
	    c("integer",	# nthreads
		"factor",	# name
		"factor",	# scheme
                "factor",	# compiler
		"factor",	# runtime
		"factor",	# model
		"integer",	# npart
		"integer",	# height
		"integer",	# bsize
                "integer",      # niter
		"numeric",	# time
                "integer"       # rmem
		))

    d$npart     <- ordered(d$npart)
    d$height    <- ordered(d$height)
    d$bsize     <- ordered(d$bsize)

    d <- subset(d, scheme %in% s)
    if (n)
        d <- subset(d, npart == n)
    if (h)
        d <- subset(d, height == h)
    d <- subset(d, model == m)

    d <- within(d, "compiler_runtime" <- paste(compiler, runtime, sep='-'))

    names <- c(
        'GCC-OMP4',
        'Klang',
        'Klang-C',
        'Klang-CP',
        'StarPU',
        'StarPU-C',
        'StarPU-CP'
    )
    d <- subset(d, name %in% names)

    return (d)
}

# OMP Compiler/runtime breaks, colors...
get_breaks_runtime <- function()
{
    return (c('implicit', 'explicit'))
}

get_labels_runtime <- function()
{
    return (c('Implicit', 'Explicit'))
}

get_colors_runtime <- function()
{
    return (c('implicit'  = "#266d83",
              'explicit'   = "#e20025"))
}

# Scheme breaks, colors ...
get_breaks_scheme <- function()
{
    return (c('tb-omp4#task#dep', 'tb-omp4#task#dep-P','tb-omp4#task#dep-C',
              'tb-omp4#task#dep-CP'))
}

get_shapes_scheme <- function()
{
    return(c('tb-omp4#task#dep'    = 0,
             'tb-omp4#task#dep-P'  = 2,
             'tb-omp4#task#dep-C'  = 10,
             'tb-omp4#task#dep-CP' = 8))
}

get_ltypes_scheme <- function()
{
    return (c('tb-omp4#task#dep'    = "solid",
              'tb-omp4#task#dep-P'  = "dotdash",
              'tb-omp4#task#dep-C'  = "solid",
              'tb-omp4#task#dep-CP' = "dashed"))
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
