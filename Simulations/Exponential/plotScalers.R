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
library(ggpubr)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the parameters
params = read.table("rates.txt", header=TRUE, sep="\t", nrows = 0)
param_labels = labels(params)[[2]]

# get the names of all SISCO first (of three) run log files
log <- list.files(path="../out", pattern="Expo.*.log", full.names = TRUE)

states <- 3

data = data.frame()

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  
  # Read in the SISCO *.logs
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$m1)/10)), ]

  # calculate ess values
  if (length(t$Sample>100)){
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
                                        lower=hpd[1,"lower"], upper=hpd[1,"upper"], run=runnumber))
        }
      }
  }
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the rate ratios
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tot.ne = data[data$name=="NeNull","lower"]<data[data$name=="NeNull","true"] &
  data[data$name=="NeNull","upper"]>data[data$name=="NeNull","true"]
cov.ne = sum(tot.ne[!is.na(tot.ne)])/sum(!is.na(tot.ne))

tot.growth = data[data$name=="GrowthRate","lower"]<(-1)*data[data$name=="GrowthRate","true"] &
  data[data$name=="GrowthRate","upper"]>(-1)*data[data$name=="GrowthRate","true"]
cov.growth = sum(tot.growth[!is.na(tot.growth)])/sum(!is.na(tot.growth))


tot.mig = data[data$name=="f_migrationRatesSkyline","lower"]<exp(data[data$name=="f_migrationRatesSkyline","true"]) &
  data[data$name=="f_migrationRatesSkyline","upper"]>exp(data[data$name=="f_migrationRatesSkyline","true"])
cov.mig = sum(tot.mig[!is.na(tot.mig)])/sum(!is.na(tot.mig))


p_ne <- ggplot(data[data$name=="NeNull",]) +
  geom_point(aes(x=true, y=median))+
  geom_abline(color="red") +
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true")+
  ggtitle(paste("log(present Ne) (cov =", round(cov.ne, 2), ")")) +
  
  theme_minimal()+
  theme(legend.position="none")
plot(p_ne)

p_growth <- ggplot(data[data$name=="GrowthRate",]) +
  geom_point(aes(x=true, y=-median))+
  geom_abline(color="red") +
  geom_errorbar(aes(x=true, ymin=-lower, ymax=-upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + theme_minimal()+
  ggtitle(paste("Growth rate (cov =", round(cov.growth, 2), ")")) +
  theme(legend.position="none")
plot(p_growth)

p_migration <- ggplot(data[data$name=="f_migrationRatesSkyline",]) +
  geom_point(aes(x=exp(true), y=median))+
  geom_errorbar(aes(x=exp(true), ymin=lower, ymax=upper), alpha=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  ylab("estimated") + xlab("true") + 
  ggtitle(paste("Migration Rates (cov =", round(cov.mig, 2), ")")) +
  theme_minimal()+
  theme(legend.position="none")
plot(p_migration)


p1 = ggarrange(p_ne, p_growth, p_migration, ncol = 3, labels = c("A", "B", "C"))
ggsave(plot=p1, filename = "../../../MascotSkyline-Text/Figures/exponential.pdf", height=4, width=10)

# check for outliers
tmp = data[data$name=="GrowthRate",]
tmp$diff = abs(tmp$median+tmp$true)
print(tmp[tmp$diff>0.5, "run"])

