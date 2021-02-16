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
library(colorblindr)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the parameters
params = read.table("rates.txt", header=TRUE, sep="\t", nrows = 0)
param_labels = labels(params)[[2]]

# get the names of all SISCO first (of three) run log files
log <- list.files(path="./out", pattern="Sky.*.log", full.names = TRUE)

states <- 2

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
  t <- t[-seq(1,ceiling(length(t$m1)/2)), ]

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
          if (startsWith(param_labels[[i]], "Ne")){
            name_log = gsub("Ne0", "NeLog", param_labels[[i]])
            name_log = gsub("state0.", "state0", name_log)
            name_log = gsub("state1.", "state1", name_log)
            
            # get state
            if (grepl("state0", name_log)){
              state="state 1"
              time = as.numeric(strsplit(name_log, split="state0")[[1]][[2]])
            }else if (grepl("state1", name_log)){
              state= "state 2"
              time = as.numeric(strsplit(name_log, split="state1")[[1]][[2]])
            }

            median_val = median(t[,name_log])
            hpd = HPDinterval(as.mcmc(t[,name_log]))
            new.data = data.frame(true=params[param_index,param_labels[[i]],], median=median_val,lower=hpd[1,"lower"],upper=hpd[1,"upper"], runnr=runnumber,time=time,state=state)
            if (start.Ne){
              Ne = new.data
              start.Ne = F
            }else{
              Ne = rbind(Ne, new.data)
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


# compute the coverage

Ne$inint = 0
Ne[which(Ne$true>=Ne$lower & Ne$true<Ne$upper), "inint"] = 1
mig$inint = 0
mig[which(mig$true>=mig$lower & mig$true<mig$upper), "inint"] = 1
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the rate ratios
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

p_ne <- ggplot()+
  geom_point(data=Ne, aes(x=true, y=median))+
  geom_segment(data=Ne, aes(x = -5, y = -5, xend = 5, yend = 5), color="red") +
  geom_errorbar(data=Ne, aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle(paste("Ne, coverage =", mean(Ne$inint))) +
  theme(legend.position="none")
plot(p_ne)


p_growth <- ggplot()+
  geom_point(data=mig, aes(x=true, y=median))+
  geom_segment(data=mig, aes(x = 0, y = 0, xend = 5, yend = 5), color="red") +
  geom_errorbar(data=mig, aes(x=true, ymin=lower, ymax=upper), alpha=0.5)+
  ylab("estimated") + xlab("true") + ggtitle(paste("migration rates, coverage =", mean(mig$inint))) +
  scale_x_log10(limits=c(0.01,5))+
  scale_y_log10(limits=c(0.01,5))+
  theme(legend.position="none")
plot(p_growth)


# plot the skyline graphs
Ne_plot = Ne[which(Ne$runnr<13 & Ne$runnr>0),]

p_skyline <- ggplot(data=Ne_plot)+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper,fill=state, group=state), alpha=0.2)+
  geom_line(aes(x=time, y=median,color=state, group=state, linetype = "median estimate"))+
  geom_line(aes(x=time, y=true,color=state, group=state, linetype = "true value"))+
  facet_wrap(runnr~., ncol=5) +
  ylab("log Ne") + 
  xlab("time") + 
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme_minimal()  +
  scale_color_OkabeIto()+
  scale_fill_OkabeIto()

plot(p_skyline)

ggsave(plot=p_skyline,paste("../../Figures/skygrid_trends.pdf", sep=""),width=9, height=4.5)


