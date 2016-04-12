library(plyr)
library(ggplot2)
library(scales)

calc_speedup <- function(data, ref_algo)
{
    # XXX: probably suboptimal
    data_ref <- subset(data, algo == ref_algo)
    for (i in 1:nrow(data)) {
        tmp_ref <- subset(data_ref, npart == data$npart[i] & height == data$height[i] & nnode == data$nnode[i])
		#tmp_ref <- subset(tmp_ref, nthreads == data$nthreads[i])
        data$speedup[i] <- tmp_ref$global_time / data$global_time[i]
    }
    return (data)
}

gen_speedup_taskdep_plot <- function(d, model_wanted)
{
	d <- subset(d, model == model_wanted)
	d <- calc_speedup(d, "explicit")
	g <- ggplot(data=d,aes_string(x="nnode", y="speedup", color="algo"))
    g <- g + geom_line()
	g <- g + facet_wrap(npart ~ height, scales="free",
						labeller = labeller(npart = as_labeller(npart_labeller),
											height = as_labeller(height_labeller),
											.default=label_both,
											.multi_line=FALSE))

    # Set our own colors, linetypes and point shapes.
	g <- g + scale_color_manual(name="Algorithm",
								breaks=get_breaks_runtime(),
								labels=get_labels_runtime(),
								values=get_colors_runtime())

    # Set X/Y labels.
	g <- g + xlab("Number of nodes")
	g <- g + ylab("Speedup")

    # Set y-axis range
    #g <- g + ylim(ylimits)

    # Save generated plot.
	output <- paste(model_wanted, "-speedup.pdf", sep="")
	ggsave(output, g, width=29.7, height=21, units=c("cm"), device=cairo_pdf)
}

gen_speedup <- function(dbfile)
{
    data <- get_data_subset(dbfile, 0L, 0L, "False")
    #output <- paste(output_dir, node, "-1M-7-cube-speedup.pdf", sep="")

	all_model <- unique(data$model)
	for (i in 1:length(all_model))
	{
		gen_speedup_taskdep_plot(data, all_model[i])
	}
}
