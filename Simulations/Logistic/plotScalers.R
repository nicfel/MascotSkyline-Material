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
log <- list.files(path="./out", pattern="*.log", full.names = TRUE)

states <- 3

start.Ne <- T
start.growth <- T
start.mig <- T

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  
  # Read in the SISCO *.logs
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$m1)/10)), ]

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
      for (i in seq(2,length(param_labels))){
          if (startsWith(param_labels[[i]], "carrying")){
            name_log = gsub("carrying", "carrying", param_labels[[i]])
            median_val = median(t[,name_log])
            hpd = HPDinterval(as.mcmc(t[,name_log]))
            new.data = data.frame(true=params[param_index,param_labels[[i]],], median=median_val,lower=hpd[1,"lower"],upper=hpd[1,"upper"])
            if (start.Ne){
              Ne = new.data
              start.Ne = F
            }else{
              Ne = rbind(Ne, new.data)
            }
          }
      }
      # loop over the logits's
      for (i in seq(2,length(param_labels))){
        if (startsWith(param_labels[[i]], "logit")){
          name_log = gsub("logit", "logit", param_labels[[i]])
          median_val = median(t[,name_log])
          hpd = HPDinterval(as.mcmc(t[,name_log]))
          new.data = data.frame(true=params[param_index,param_labels[[i]],], median=median_val,lower=hpd[1,"lower"],upper=hpd[1,"upper"])
          if (start.Ne){
            logit = new.data
            start.Ne = F
          }else{
            logit = rbind(Ne, new.data)
          }
        }
      }
      
      
      # loop over the growth's
      for (i in seq(2,length(param_labels))){
        if (startsWith(param_labels[[i]], "growth")){
          name_log = gsub("Ne0", "Nenull", param_labels[[i]])
          median_val = median(t[,name_log])
          hpd = HPDinterval(as.mcmc(t[,name_log]))
          new.data = data.frame(true=-params[param_index,param_labels[[i]],], median=median_val,lower=hpd[1,"lower"],upper=hpd[1,"upper"])
          if (start.growth){
            growth = new.data
            start.growth = F
          }else{
            growth = rbind(growth, new.data)
          }
        }
      }
      
      # loop over the migrations
      c = 6
      for (a in seq(1,states)){
        for (b in seq(1,states)){
          if (a!=b){
            mig_label = paste("f_migration.state",a-1,"_to_state",b-1,sep="")
            migration = paste("f_migration.state",a-1,"_to_state",b-1,sep="")
            median_val = median(t[,migration])
            hpd = HPDinterval(as.mcmc(t[,migration]))
            new.data = data.frame(true=params[param_index,mig_label], median=median_val,lower=hpd[1,"lower"],upper=hpd[1,"upper"])
            if (start.mig){
              mig = new.data
              start.mig = F
            }else{
              mig = rbind(mig, new.data)
            }
            c = c-1
          }
        }
      }
  
  }
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the rate ratios
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

p_ne <- ggplot()+
  geom_point(data=Ne, aes(x=true, y=median))+
  geom_segment(data=Ne, aes(x = -2, y = -2, xend = 2, yend = 2), color="red") +
  geom_errorbar(data=Ne, aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle("carrying capacity") +
  theme(legend.position="none")
plot(p_ne)

p_growth <- ggplot()+
  geom_point(data=growth, aes(x=true, y=median))+
  geom_segment(data=growth, aes(x = -0.2, y = -0.2, xend = 1.5, yend = 1.5), color="red") +
  geom_errorbar(data=growth, aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle("growth rate") +
  theme(legend.position="none")
plot(p_growth)

p_growth <- ggplot()+
  geom_point(data=mig, aes(x=true, y=median))+
  geom_segment(data=Ne, aes(x = -2, y = -2, xend = 2, yend = 2), color="red") +
  geom_errorbar(data=mig, aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle("migration rates") +
  theme(legend.position="none")
plot(p_growth)

p_logit <- ggplot()+
  geom_point(data=logit, aes(x=true, y=median))+
  geom_segment(data=Ne, aes(x = -2, y = -2, xend = 2, yend = 2), color="red") +
  geom_errorbar(data=logit, aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle("migration rates") +
  theme(legend.position="none")
plot(p_logit)





