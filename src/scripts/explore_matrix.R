# Title     : TODO
# Objective : TODO
# Created by: jinlf
# Created on: 2/7/18

getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}
get_loci_data <- function(eg_file, loci) {
    src = read.table(eg_file,head=T, row.names=1)
    src.t = as.data.frame(t(as.matrix(src)))
    loci_data = src.t[[loci]]
    samples = c()
    for(i in 1:length(loci_data)){
        tmp=rep(i-1,loci_data[i])
        samples=c(samples,tmp)
    }

    return(samples)
}

desc <- function(loci_data) {
    loci_mode = getmode(loci_data)
    print(paste("mode", ": ", loci_mode))
    print("====")
    print(summary(loci_data))
}

data_file = "/home/jinlf/test/msi-test/msidetector_test/20180207/RND-NA18548-P4_170808_003958_msi_locis_freq_matrix.txt"
loci = "chr9-135773000-135773018-A18"
loci_data = get_loci_data(data_file,loci)
desc(loci_data)


