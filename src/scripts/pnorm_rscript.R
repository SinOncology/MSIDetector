# Title     : TODO
# Objective : TODO
# Created by: jinlf
# Created on: 9/20/17

args <- commandArgs(T)
query = args[1]
mean_value = args[2]
sd_value = args[3]
mode(query)= 'numeric'
mode(mean_value)= 'numeric'
mode(sd_value)= 'numeric'

pvalue = (pnorm(q=query, mean=mean_value,sd = sd_value,log.p = F))
if(query > mean_value){
    pvalue = 1-pvalue
}
log_pvalue = -log10(pvalue)
print(pvalue)
print(log_pvalue)
