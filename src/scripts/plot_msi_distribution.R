# Title     : TODO
# Objective : TODO
# Created by: jinlf
# Created on: 9/21/17
library(ggplot2)
args <- commandArgs(T)
matrix_file = args[1]
area_file = args[2]
out_pdf = args[3]
msi_data = read.csv(matrix_file, sep = '\t', row.names = 1)
area_data = read.csv(area_file, sep = '\t', row.names = 1)
index_vec = seq(0, 50)
mode(index_vec) = 'character'
colnames(msi_data) = index_vec
msi_data = as.data.frame(t(as.matrix(msi_data)))

plot_list_data = list()
loci_names = names(msi_data)
for (loci in loci_names)
{
    start_index = as.numeric(area_data[loci,][1])
    end_index = as.numeric(area_data[loci,][2])
    plot_list_data[[loci]] = data.frame(repeatCount = as.vector(seq(start_index, end_index)), rate = msi_data[, loci][start_index : end_index])
}

pdf(out_pdf)

for (loci in loci_names)
{
#     loci = "A14"
    tmp_loci_data = plot_list_data[[loci]]
    Sys.sleep(1)
    plot_title_name = paste(c(loci, "repeatCount", "distribution_plot"), collapse = "_")
    p=ggplot(tmp_loci_data, aes(x = tmp_loci_data[, 1], y = tmp_loci_data[, 2])) +
        geom_line(colour = "#78a355") +
        geom_point(colour = "red", size = 1.15) +
        geom_text(aes(label = tmp_loci_data[, 2]), hjust = 1.0, vjust = - 1.0, size = 1.1) +
        xlab("repeat Count") +
        ylab("support rate") +
        labs(title = plot_title_name) +
        theme(plot.title = element_text(size = 7, colour = "blue", face = "bold", hjust = 0.5), axis.title.x = element_text(size = 5.7, colour = "black", face = "bold"), axis.text.x = element_text(colour = "red", size = 4.3), axis.title.y = element_text(size = 5.7, colour = "black", face = "bold"), axis.text.y = element_text(colour = "red", size = 4.3)) +
        scale_x_continuous(breaks = seq(0, 100 , 1))
    plot(p)


    # ggsave(scatter, file = out_pdf, width = 50, height = 50, limitsize = F)
}
dev.off()