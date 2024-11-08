######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
# needed to calculate ESS values
library(coda)
library("methods")


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the parameters
params = read.table("rates.txt", header=TRUE, sep="\t", nrows = 0)
param_labels = labels(params)[[2]]

# get the names of all SISCO first (of three) run log files
log <- list.files(path="../out", pattern="Log.*.log", full.names = TRUE)

data = data.frame()

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  
  # Read in the SISCO *.logs
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$Sample)/10)), ]

  # calculate ess values
  ess <- effectiveSize(t)

    if (min(ess[2:3])<100){
    print("masco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(ess[2:6]),filename1))
  }else{
      dfname <- data.frame(filename = filename1)
      
      # get the number of the run
      tmp = strsplit(filename1, split="_")
      tmp = strsplit(tmp[[1]][2], split="S")
      
      runnumber = as.numeric(tmp[[1]][2])
      param_index = which(params$run==runnumber)
      
      # loop over the Ne0's
      for (j in seq(2,length(param_labels))){
        mean = mean(t[,param_labels[j]])
        median = mean(t[,param_labels[j]])
        hpd = HPDinterval(as.mcmc(t[,param_labels[j]]))
        value_name = strsplit(param_labels[j], split="\\.")[[1]][[1]]
        true=params[param_index, param_labels[j]]
        
        data = rbind(data, data.frame(name = value_name, true=true, mean=mean, median=median, 
                                      lower=hpd[1,"lower"], upper=hpd[1,"upper"]))
      }
    
  }
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the rate ratios
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

p_ne <- ggplot(data[data$name=="Capacity",]) +
  geom_point(aes(x=true, y=median))+
  geom_abline(color="red") +
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle("present Ne") + theme_minimal()+
  theme(legend.position="none")
plot(p_ne)

p_ne2 <- ggplot(data[data$name=="CarryingProportion",]) +
  geom_point(aes(x=true, y=median))+
  geom_abline(color="red") +
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle("present Ne") + theme_minimal()+
  theme(legend.position="none")
plot(p_ne2)


p_growth <- ggplot(data[data$name=="GrowthRate",]) +
  geom_point(aes(x=true, y=-median))+
  geom_abline(color="red") +
  geom_errorbar(aes(x=true, ymin=-lower, ymax=-upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle("growth rate") + theme_minimal()+
  theme(legend.position="none")
plot(p_growth)

p_migration <- ggplot(data[data$name=="f_migrationRatesSkyline",]) +
  geom_point(aes(x=exp(true), y=median))+
  geom_errorbar(aes(x=exp(true), ymin=lower, ymax=upper), alpha=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  ylab("estimated") + xlab("true") + ggtitle("migration rates") + theme_minimal()+
  theme(legend.position="none")
plot(p_migration)



p1 = ggarrange(p_ne, p_ne2, p_growth, p_migration, ncol = 2, labels = c("A", "B", "C"))
ggsave(plot=p1, filename = "../../../MascotSkyline-Text/Figures/exponential.pdf", height=8, width=9)

