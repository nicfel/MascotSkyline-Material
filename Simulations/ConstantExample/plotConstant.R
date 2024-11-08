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
log <- list.files(path="./out", pattern="*.log", full.names = TRUE)

data = data.frame()

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  # Read in the SISCO *.logs
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$Sample)/10)), ]
  
  if (length(t$Sample)>10){
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
        label = gsub("\\.state",".",param_labels[j])
        mean = mean(t[,label])
        median = mean(t[,label])
        hpd = HPDinterval(as.mcmc(t[,label]))
        value_name = strsplit(param_labels[j], split="\\.")[[1]][[1]]
        true=params[param_index, param_labels[j]]
        
        if (grepl("Ne", value_name)){
          data = rbind(data, data.frame(name = value_name, true=true, mean=mean, median=median, 
                                        lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                        time = as.numeric(strsplit(param_labels[j], split="\\.")[[1]][[3]]),
                                        state = strsplit(param_labels[j], split="\\.")[[1]][[2]],
                                        run=runnumber))
        }else{
          data = rbind(data, data.frame(name = value_name, true=true, mean=mean, median=median, 
                                        lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                        time = NA,
                                        state=NA,
                                        run=runnumber))
        }
      }
    }
  }
}


set.seed(1223)

runs=unique(data$run)
selected_runs <- sample(runs, 4)

# plot the skyline graphs
p_skyline <- ggplot(data=data[data$name=="SkylineNe" & data$run %in% selected_runs & data$state!="combined",], aes(x=time, color=state, fill=state, group=state))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
  geom_line(aes(y=median,linetype = "median estimate"))+
  geom_line(aes(y=true, linetype = "true value"))+
  facet_wrap(run~., ncol=2) +
  ylab("log Ne") + 
  xlab("time") + 
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme_minimal()  +
  scale_color_OkabeIto()+
  scale_fill_OkabeIto()
labels <- setNames(paste("replicate", 1:length(selected_runs)), sort(selected_runs))

# Use the labels in facet_wrap
p_skyline <- p_skyline + 
  facet_wrap(~ run, ncol=2, labeller = labeller(run = labels))

plot(p_skyline)

ggsave(plot=p_skyline,paste("../../../MascotSkyline-Text/Figures/constant_skygrid_trends.pdf", sep=""),width=8, height=6)

