library(plyr)
library(ggplot2)
library(scales)

gen_group_size_plot_comm <- function(d, model_wanted)
{
	d <- subset(d, model == model_wanted)
	d$global_time <- d$global_time/1000
	d <- d[d$communication_vol >= 0,]
	if(nrow(d) > 0)
	{
	g <- ggplot(data=d,aes_string(x="bsize", y="communication_vol", color="algo", group="algo"))
    g <- g + geom_line()
	g <- g + facet_wrap(npart ~ height ~ nnode, scales="free",
						labeller = labeller(npart = as_labeller(npart_labeller),
											height = as_labeller(height_labeller),
											nnode = as_labeller(nnode_labeller),
											.default=label_both,
											.multi_line=FALSE))

    # Set our own colors, linetypes and point shapes.
	g <- g + scale_color_manual(name="Algorithm",
								breaks=get_breaks_runtime(),
								labels=get_labels_runtime(),
								values=get_colors_runtime())

    # Set X/Y labels.
	g <- g + xlab("bsize")
	g <- g + ylab("Communication volume (MB)")

    # Save generated plot.
	output <- paste(get_output_directory(), "/", model_wanted, "-bsize-volume.pdf", sep="")
	ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
	}
}
gen_group_size_plot_speed <- function(d, model_wanted)
{
	d <- subset(d, model == model_wanted)
	d$global_time <- d$global_time/1000
	d <- d[d$global_time > 0,]
	if(nrow(d) > 0)
	{
	g <- ggplot(data=d,aes_string(x="bsize", y="global_time", color="algo", group="algo"))
    g <- g + geom_line()
	g <- g + facet_wrap(npart ~ height ~ nnode, scales="free",
						labeller = labeller(npart = as_labeller(npart_labeller),
											height = as_labeller(height_labeller),
											nnode = as_labeller(nnode_labeller),
											.default=label_both,
											.multi_line=FALSE))

    # Set our own colors, linetypes and point shapes.
	g <- g + scale_color_manual(name="Algorithm",
								breaks=get_breaks_runtime(),
								labels=get_labels_runtime(),
								values=get_colors_runtime())

    # Set X/Y labels.
	g <- g + xlab("bsize")
	g <- g + ylab("Time (s)")
	#g <- g + scale_x_continuous(breaks=c(8, 16, 32, 64, 128, 512, 1024, 2048, 4096))

    # Save generated plot.
	output <- paste(get_output_directory(), "/", model_wanted, "-bsize-speed.pdf", sep="")
	ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
	}
}

gen_group_size <- function(dbfile)
{
    data <- get_data_subset(dbfile, 0L, 0L, "False", 0L)

	data <- subset(data, algo != get_one_node_reference_algorithm() & algo != "simple-mpi")
	all_nnode <- unique(subset(data, bsize != get_bsize_reference())$nnode)
	data <- subset(data, nnode %in% all_nnode)

	all_model <- unique(data$model)
	for (i in 1:length(all_model))
	{
		gen_group_size_plot_speed(data, all_model[i])
		gen_group_size_plot_comm(data, all_model[i])
	}
}
