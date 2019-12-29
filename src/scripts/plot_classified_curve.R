# Title     : TODO
# Objective : TODO
# Created by: jinlf
# Created on: 1/16/18

library(ggplot2)
args <- commandArgs(T)
baseline = args[1]
mss = args[2]
msh = args[3]
plot_name = args[4]
out_pdf = args[5]
base_data = read.csv(baseline, sep = '\t')
mss_data = read.csv(mss, sep = '\t')
msh_data = read.csv(msh, sep = '\t')
# plot_name = "MONO-27_25_29"
p <- ggplot() +
    theme(plot.title = element_text(size = 10, colour = "blue", face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 7, colour = "black", face = "bold"),
    axis.text.x = element_text(colour = "red", size = 5),
    axis.title.y = element_text(size = 7, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "red", size = 5),
    legend.text = element_text(size = 5), legend.title = element_text(size = 7)) +
    scale_x_continuous(breaks = seq(0, 100 , 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    xlab("sample") +
    ylab("sum frequency  of the repeat counts") +
    labs(title = plot_name)
p <- p +
    geom_line(data = base_data, aes(y = base_data[, 3], x = base_data[, 1], color = "green"), size = 0.5) +
    geom_point(data = base_data, aes(x = base_data[, 1] , y = base_data[, 3], color = "green"), size = 0.5) +
    geom_text(aes(x = base_data$X, y = base_data$rate, label = as.vector(base_data$sample)) , hjust = - 0.2, vjust = - 0.2,
    size = 1.25)
p <- p +
    geom_line(data = mss_data, aes(y = mss_data[, 3], x = mss_data[, 1], color = "blue"), size = 0.5) +
    geom_point(data = mss_data, aes(x = mss_data[, 1] , y = mss_data[, 3], color = "blue"), size = 0.5) +
    geom_text(aes(x = mss_data$X, y = mss_data$rate, label = as.vector(mss_data$sample)), hjust = - 0.2, vjust = - 0.2,
    size = 1.25)
p <- p +
    geom_line(data = msh_data, aes(y = msh_data[, 3], x = msh_data[, 1], color = "purple"), size = 0.5) +
    geom_point(data = msh_data, aes(x = msh_data[, 1] , y = msh_data[, 3], color = "purple"), size = 0.5) +
    geom_text(aes(x = msh_data$X, y = msh_data$rate, label = as.vector(msh_data$sample)), hjust = - 0.2, vjust = - 0.2,
    size = 1.25)
p <- p + scale_color_discrete(name = "sample class", labels = c("msi-stable", "normal-sample", "msi-h"))
ggsave(p, file = out_pdf, width = 8.5, height = 5, limitsize = F)