# load the libraries
library(ggplot2)
library(coda)
library("methods")
library(colorblindr)
library(ggtree)
library(treeio)
library(ggpubr)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the parameters from the rates.txt file
params = read.table("rates.txt", header=TRUE, sep="\t", nrows = 0)
# get all the headers
param_labels = labels(params)[[2]]

log <- list.files(path="../out", pattern="SIR.*1mascot.log", full.names = TRUE)

# define the number of states in the params file
states <- 2

# initialize the data frames
data = data.frame()
mig = data.frame()
incidence = data.frame()
truemig = data.frame()

# get the full trees to later compute the number of migration events
truetrees <- list.files(path="./master", pattern="*.single.tree", full.names = TRUE)

# initialize the plots
plots = list()
pl = 1

# loop over the trees files
for (tr in seq(1, length(truetrees))){
  print(truetrees[[tr]]) # print the number of the tree file
  trees = read.beast(truetrees[[tr]]) # read in the tree file
  edges = trees@phylo$edge # get the edges of the tree
  nodes = as.numeric(trees@data$node) # get the nodes of the tree
  
  tmp = strsplit(truetrees[[tr]], split="_") # get the run number
  tmp = strsplit(tmp[[1]][2], split="S") # get the run number  
  runnumber = as.numeric(tmp[[1]][2]) # get the run number

  # if (runnumber>100){ # if the run number is larger than 50, it is a higher migration run
  #   sampling="lower migration" # set the sampling to higher migration
  #   if (runnumber<106){ 
  #     plots[[runnumber-95]] = ggtree(trees, aes(color=location), size=0.25) + scale_color_OkabeIto() +theme(legend.position = "none") # plot the tree
  #   }
  # }else{
  #   sampling="higher migration"
  #   if (runnumber<6){
  #     plots[[runnumber]] = ggtree(trees, aes(color=location), size=0.25) + scale_color_OkabeIto() +theme(legend.position = "none")
  #   }
  # }
  
  nosamples = paste(params[params$run==runnumber, "samples"][[1]],"samples")
  sampling = params[params$run==runnumber, "avg_migration"][[1]]
  r0 = params[params$run==runnumber, "random_r0"][[1]]
  
  events = 0; # initialize the number of migration events
  for (i in seq(1, length(edges[,1]))){ # loop over all the edges
    parent_edge = which(edges[,2]==edges[i,1]) # get the parent edge
    
    if (length(parent_edge)==1){ # if there is a parent edge
      from = which(nodes==edges[i,1]) # get the node number of the parent
      to = which(nodes==edges[i,2]) # get the node number of the child
      if (trees@data$location[from]!=trees@data$location[to]){ # if the parent and child are in different locations
        events=events+1 # increase the number of migration events
      }
    }
  }
  truemig = rbind(truemig, 
                  data.frame(events=events, tips = length(trees@phylo$tip.label), 
                             sampling=sampling, 
                             nosamples=nosamples,
                             random_ro=r0,
                             run=runnumber, 
                             rootTime=min(as.numeric(trees@data$time)))) # add the number of migration events to the data frame
}
  
# p.trees1 = ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]], ncol = 5, labels = c("F", "G", "H", "I", "J"))
# p.trees2 = ggarrange(plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]], ncol = 5, labels = c("F", "G", "H", "I", "J"))
# 
# p.trees = ggarrange(p.trees1, p.trees2, ncol=1)
# 
# plot(p.trees)
# ggsave(plot=p.trees, file="../../../MascotSkyline-Text/Figures/SIR_trees.pdf", width = 9, height=6)


# support = read.table("nodesupport.tsv", header=T, sep="\t") # read in the node support file
# support$ratio=c() # initialize the ratio of migration events per tip
# 
# uni_runs = unique(support$run) # get the unique run numbers
# for (i in seq(1, length(uni_runs))){ # loop over the unique run numbers
#   ind = which(truemig$run==uni_runs[[i]]) # get the index of the run number in the truemig data frame
#   ind2 = which(support$run==uni_runs[[i]]) # get the index of the run number in the support data frame  
#   support[ind2, "ratio"] = truemig[ind, "events"]/truemig[ind, "tips"] # add the ratio of migration events per tip to the support data frame
# }
# 
# 
# for (i in seq(1, length(support$dsa))){
#   support[i,"support"] = paste(params[params$run==support[i,"run"], "samples"][[1]],"samples")
#   support[i,"sampling"] = paste(params[params$run==support[i,"run"], "samples"][[1]],"samples")
# }
# support[support$method=="dta", "method"]="DTA"
# support[support$method=="mascot", "method"]="MASCOT-Skyline"
# 
# 
# p.density = ggplot(support, aes(x=support, fill=method, color=method)) + 
#   geom_histogram(alpha=0.4, position = 'identity', bins=5) + 
#   theme_minimal() + 
#   # scale_y_log10() +
#   facet_grid(sampling~.)+
#   scale_fill_OkabeIto() +  scale_color_OkabeIto() +xlab("posterior support for true node location") + ylab("count")+
#   scale_color_manual(values = c("MASCOT-Skyline" = "skyblue", "DTA" = "darkorange")) +
#   scale_fill_manual(values = c("MASCOT-Skyline" = "skyblue", "DTA" = "darkorange"))
# 
# plot(p.density)
# ggsave(plot=p.density, file="../../../MascotSkyline-Text/Figures/SIR_nodesupport.pdf", width = 9, height=6)
# 
# # plot the mean node support per run on the y-axis to the ratio of migration events per tip on the x-axis
# p.support = ggplot(support, aes(x=ratio, y=support, color=method)) + 
#   scale_x_log10() +facet_grid(sampling~.)+theme_minimal() + geom_smooth() +
#   scale_color_OkabeIto() +xlab("migration events per tip") + ylab("posterior support for true node location")
#   
# 
# plot(p.support)


# = ggplot(support, aes(x=ratio, y=support)) + 
#   geom_point(alpha=0.1) +scale_x_log10() +facet_grid(method~.)+theme_minimal()



# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  # Read in the SISCO *.logs
  t1 <- read.table(filename1, header=TRUE, sep="\t")
  # t2 <- read.table(gsub("1mas", "2mas", filename1), header=TRUE, sep="\t")
  # t3 <- read.table(gsub("1mas", "3mas", filename1), header=TRUE, sep="\t")
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
  # t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
  # t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]
  
  t=t1
  # t = rbind(t1, t2, t3)
  # calculate ess values
  if (length(t$posterior)>20){
    ess <- effectiveSize(t)
  
    if (min(ess[2:3])<100){
      print(sprintf("ESS value is %f for file %s",min(ess[2:3]),filename1))
    }else{
        dfname <- data.frame(filename = filename1)
        
        # get the number of the run
        tmp = strsplit(filename1, split="_")
        tmp = strsplit(tmp[[1]][2], split="S")
        
        runnumber = as.numeric(tmp[[1]][2])
        
        require("rjson")
        
        nosamples = paste(params[params$run==runnumber, "samples"][[1]],"samples")
        sampling = params[params$run==runnumber, "avg_migration"][[1]]
        r0 = params[params$run==runnumber, "random_r0"][[1]]
        
        last_time = 0
        
        jsondat <- fromJSON(file=gsub("1mascot.log","master.json", gsub("../out/","./master/",log[[i]])))
        tot_incidence = c()
        times = jsondat$t
        for (j in seq(1, length(jsondat$I))){
          incidence = rbind(incidence, data.frame(time=times, y=jsondat$I[[j]], 
                                                  state=paste("state", j-1, sep=""),
                                                  sampling=sampling,
                                                  nosamples=nosamples,
                                                  random_ro=r0,
                                                  run=runnumber))
          tot_incidence = c(tot_incidence, sum(diff(jsondat$I[[j]])==1))
        }
        
  
        
        param_index = which(params$run==runnumber)
        param_index = param_index[length(param_index)]
        
        last_time = truemig[which(truemig$run==runnumber), "rootTime"]

        # loop over the Ne0's if the the run number is x01 x02 x03 x04 x05 (take modulo 100)
        if (runnumber%%100<6 & runnumber%%100>0){
          for (j in seq(1, length(jsondat$I))){
            for (k in seq(1, 26)){
              name = paste("SkylineNe.state", j-1, ".", k, sep="")
              mean = mean(t[,name])
              median = mean(t[,name])
              hpd = HPDinterval(as.mcmc(t[,name]))
              value_name = strsplit(name, split="\\.")[[1]][[1]]
              data = rbind(data, data.frame(time=(last_time+t$Tree.height[1]) - (k-1)/25*t$Tree.height[1], name = value_name, mean=mean, median=median, 
                                            lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                            state=paste("state", j-1, sep=""),
                                            sampling=sampling,
                                            nosamples=nosamples,
                                            random_ro=r0,
                                            run=runnumber))
              
            }
          }
        }       
        
        # loop over the Ne0's
        for (j in seq(7,length(param_labels))){
          mean = mean(t[,param_labels[j]])
          median = mean(t[,param_labels[j]])
          hpd = HPDinterval(as.mcmc(t[,param_labels[j]]))
          value_name = strsplit(param_labels[j], split="\\.")[[1]][[1]]
          true=params[param_index, param_labels[j]]
          from = paste("R0",strsplit(param_labels[j], split="[._]")[[1]][[3]], sep=".")
          to = paste("R0",strsplit(param_labels[j], split="[._]")[[1]][[5]], sep=".")
          
          from_cp = tot_incidence[[as.numeric(gsub("state","",strsplit(param_labels[j], split="[._]")[[1]][[3]])) +1]]
          to_cp = tot_incidence[[as.numeric(gsub("state","",strsplit(param_labels[j], split="[._]")[[1]][[5]])) +1]]
          
          mig = rbind(mig, data.frame(true=true, mean=mean, median=median, 
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                          sampling=sampling,
                                          nosamples=nosamples,
                                          random_ro=r0,
                                          from_state=strsplit(param_labels[j], split="[._]")[[1]][[3]],
                                          to_state = strsplit(param_labels[j], split="[._]")[[1]][[5]],
                                          from = params[param_index, from],
                                          to = params[param_index, to],
                                          from_cp = from_cp,
                                          to_cp = to_cp,
                                          run=runnumber, method="MASCOT-Skyline"))
        }
        
  
    }
  }
}

p_skyline <- ggplot(data=data[data$name=="SkylineNe",], aes(x=time, color=state, fill=state, group=state))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
  geom_line(aes(y=median,linetype = "median estimate"))+
  geom_line(data=incidence[is.element(incidence$run, unique(data$run)),], aes(x=time, y=log(y/100), linetype = "true value", color=state))+
  facet_wrap(run~., ncol=5) +
  ylab("log Ne") + 
  xlab("time") +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme_minimal()  +
  scale_color_OkabeIto()+
  scale_fill_OkabeIto() +
  coord_cartesian(ylim=c(-8,3)) +
  scale_y_continuous(sec.axis = sec_axis(~ 100 * exp(.), name="Prevalence",
                                         breaks = c(1, 10, 100, 1000, 10000) # Log-transformed tick positions
  )) # Adjusted transformation for secondary y-axis

plot(p_skyline)

ggsave(plot=p_skyline, file="../../../MascotSkyline-Text/Figures/SIR_Ne_prevalence.pdf", width = 20, height=15)

# 
# p_skyline <- ggplot(data=data[data$name=="SkylineNe" & data$run>80,], aes(x=time, color=state, fill=state, group=state))+
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
#   geom_line(aes(y=median,linetype = "median estimate"))+
#   geom_line(data=incidence[incidence$run>80,], aes(x=time, y=log(y/25), linetype = "true value"), color="black")+
#   facet_grid(run~state) +
#   ylab("log Ne") +
#   xlab("time") +
#   theme(legend.position="none",
#         strip.background = element_blank(),
#         strip.text.x = element_blank()
#   ) +
#   theme_minimal()  +
#   scale_color_OkabeIto()+
#   scale_fill_OkabeIto() +
#   coord_cartesian(ylim=c(-4,8))
# 
# plot(p_skyline)




log <- list.files(path="../out", pattern="DTASIR.*1dta.log", full.names = TRUE)

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")
  
  # Read in the SISCO *.logs
  t1 <- read.table(filename1, header=TRUE, sep="\t")
  # t2 <- read.table(gsub("1dta", "2dta", filename1), header=TRUE, sep="\t")
  # t3 <- read.table(gsub("1dta", "3dta", filename1), header=TRUE, sep="\t")
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
  # t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
  # t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]
  
  t = t1
  # t = rbind(t1, t2, t3)
  # calculate ess values
  ess <- effectiveSize(t)
  
  if (min(ess[2:3])<50){
    print(sprintf("ESS value is %f for file %s",min(ess[2:3]),filename1))
  }else{
    dfname <- data.frame(filename = filename1)
    
    # get the number of the run
    tmp = strsplit(filename1, split="_")
    tmp = strsplit(tmp[[1]][2], split="S")
    
    runnumber = as.numeric(tmp[[1]][2])

    param_index = which(params$run==runnumber)
    param_index = param_index[length(param_index)]
    
    nosamples = paste(params[params$run==runnumber, "samples"][[1]],"samples")
    sampling = params[params$run==runnumber, "avg_migration"][[1]]
    r0 = params[params$run==runnumber, "random_r0"][[1]]
    
    # loop over the Ne0's
    for (j in seq(7,length(param_labels))){
      label = gsub("f_migrationRatesSkyline", "geo.rates", param_labels[j])
      label = gsub("_to_", ".", label)
      mean = mean(t[,label]*t[,"geo.meanRate"])
      median = mean(t[,label]*t[,"geo.meanRate"])
      hpd = HPDinterval(as.mcmc(t[,label]*t[,"geo.meanRate"]))
      true=params[param_index, param_labels[j]]
      from = paste("R0",strsplit(param_labels[j], split="[._]")[[1]][[3]], sep=".")
      to = paste("R0",strsplit(param_labels[j], split="[._]")[[1]][[5]], sep=".")
      from_to_cp = mig[mig$sampling==sampling & mig$nosamples==nosamples &
                       mig$random_ro==r0 & mig$run==runnumber &
                       mig$from_state==strsplit(param_labels[j], split="[._]")[[1]][[3]] & 
                       mig$to_state==strsplit(param_labels[j], split="[._]")[[1]][[5]],
                        c("from_cp", "to_cp")]
      if (nrow(from_to_cp)==0){
        from_to_cp[1,c(1,2)]=NA
      }
        
      mig = rbind(mig, data.frame(true=true, mean=mean, median=median, 
                                  lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                  sampling = sampling,
                                  nosamples=nosamples,
                                  random_ro=r0,
                                  from_state=strsplit(param_labels[j], split="[._]")[[1]][[3]],
                                  to_state = strsplit(param_labels[j], split="[._]")[[1]][[5]],
                                  from = params[param_index, from],
                                  to = params[param_index, to],
                                  from_cp = from_to_cp[[1]],
                                  to_cp = from_to_cp[[2]],
                                  run=runnumber, method="DTA"))
    }
  }
}

mig$ratio = NaN
# compute the average ratio for all the unique combinations of sampling and nosamples
for (a in unique(mig$sampling)){
  for (b in unique(mig$nosamples)){
    for (c in unique(mig$random_ro)){
      indices = which(truemig$sampling==a & truemig$nosamples==b & truemig$random_ro==c)
      ratio = mean(truemig[indices, "events" ]/truemig[indices, "tips" ])
      mig[mig$sampling==a & mig$nosamples==b & mig$random_ro==c, "ratio"] = ratio
    }
  }
}
require("ggpubr")

mig$offset = 0
mig[which(mig$method=="DTA"), "offset"]= 1

lambda=1
# Quantiles for the exponential distribution
q_exp <- function(p, lambda) {
  -log(1-p) / lambda
}
# 2.5th and 97.5th percentiles
lower_bound <- q_exp(0.025, lambda)
upper_bound <- q_exp(0.975, lambda)
# Width of the 95% interval
width_95 <- upper_bound - lower_bound
# Median of the exponential distribution
median_exp <- log(2) / lambda
# Compute the width of the 95% interval divided by the median
cv_exp <- width_95 / median_exp

mig$samplingname = "low migration"
mig[mig$sampling==25, "samplingname"] = "high migration"

# reorder the sampling facets
mig$samplingname = factor(mig$samplingname, levels=c("low migration", "high migration"))

set.seed(123)  # Set a seed for reproducibility
mig <- mig[sample(nrow(mig)), ]

text.frame = unique(mig[, c("sampling", "nosamples", "random_ro", "ratio")])
text.frame$samplingname = "low migration"
text.frame[text.frame$sampling==25, "samplingname"] = "high migration"




mig$name = paste(mig$samplingname, mig$nosamples, sep="\n")

mig[mig$run>400 & mig$run<501, "name"] = "low migration\nrandom r0\n250 samples"
mig[mig$run>500 & mig$run<601, "name"] = "low migration\nequal sampling rate\nrandom population size"
mig[mig$run>600, "name"] = "low migration\nequal sampling rate and\npopulation size"

# reorder the name facets to have the above three facets come last
mig$name = factor(mig$name, levels=c("low migration\n250 samples", "low migration\n500 samples", "high migration\n250 samples","high migration\n500 samples","low migration\nrandom r0\n250 samples", "low migration\nequal sampling rate\nrandom population size", "low migration\nequal sampling rate and\npopulation size"))

p_migration <- ggplot(mig) +
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
  geom_point(aes(x=true, y=median, color=method))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(aes(x=true, y=median, color=method, fill=method), method="lm")+
  stat_cor(aes(x=true, y=median, color=method, fill=method, label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_wrap(name~., ncol=4)+
  coord_cartesian(xlim=c(0.05,250), ylim=c(0.05,250))+
  ylab("estimated migration rate") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  theme(legend.position=c(0.9,.25))


plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates.pdf", width = 9, height=6)




corr_mig = mig[mig$method=="DTA",]
corr_mig$method = "DTA with case correction"
corr_mig$median = corr_mig$median * corr_mig$to_cp/corr_mig$from_cp
corr_mig$upper = corr_mig$upper * corr_mig$to_cp/corr_mig$from_cp
corr_mig$lower = corr_mig$lower * corr_mig$to_cp/corr_mig$from_cp
corr_mig = rbind(corr_mig, mig)

set.seed(123)  # Set a seed for reproducibility
corr_mig <- corr_mig[sample(nrow(corr_mig)), ]


p_migration <- ggplot(corr_mig, aes(x=true, y=median, ymax=upper, ymin=lower, 
                                                            color=method,fill=method)) +
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_errorbar(alpha=0.25)+
  geom_point()+
  geom_smooth(method="lm")+
  stat_cor(aes(label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_wrap(name~., ncol=4)+
  coord_cartesian(xlim=c(0.05,250), ylim=c(0.05,250))+
  ylab("estimated migration rate") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699"))+
  theme(legend.position=c(0.9,.25))
plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_correction.pdf", width = 9, height=6)







p_migration <- ggplot(mig[mig$random_ro==1,]) +
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
  geom_point(aes(x=true, y=median, color=method))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(aes(x=true, y=median, color=method, fill=method), method="lm")+
  stat_cor(aes(x=true, y=median, color=method, fill=method, label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_grid(samplingname~nosamples)+
  coord_cartesian(xlim=c(0.05,250), ylim=c(0.05,250))+
  ylab("estimated migration rate") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))
# add the ratio of migration events per tip to the plot
  geom_text(data=text.frame[text.frame$random_ro==1,], aes(x=1, y=0.1, label=paste("migration events per tip:", round(ratio, 2))), color="black", size=3, hjust=0)
plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_randomR0.pdf", width = 6, height=3.5)

p_migration <- ggplot(mig[mig$random_ro==0,],aes(x=true, y=(upper-lower)/median, color=method, fill=method)) +
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(yintercept=cv_exp, color="red")+
  geom_smooth()+
  facet_grid(samplingname~nosamples)+
  # coord_cartesian(xlim=c(0.05,250))+
  ylab("coefficient of variation") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))
plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_cv.pdf", width = 9, height=7)


p_migration <- ggplot(mig[mig$random_ro==1,],aes(x=true, y=(upper-lower)/median, color=method, fill=method)) +
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(yintercept=cv_exp, color="red")+
  geom_smooth()+
  facet_grid(samplingname~nosamples)+
  # coord_cartesian(xlim=c(0.05,250))+
  ylab("coefficient of variation") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))
plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_randomR0_cv.pdf", width = 6, height=3.5)



mig$ratio_x = 0
for (i in seq(1, length(mig$ratio))){
  mig[i, "ratio_x"] = truemig[truemig$run==mig[i, "run"], "events" ]/truemig[truemig$run==mig[i, "run"], "tips" ]
}


p_mig_ratio <- ggplot(mig[mig$random_ro==0,], aes(x=ratio_x, y=median/true, color=method, fill=method)) +
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method="lm")+
  facet_grid(samplingname~nosamples)+
  # coord_cartesian(xlim=c(0.05,250), ylim=c(0.05,250))+
  ylab("median estimated migration rate over simulated rate") + xlab("migration events in tree per tip") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))

  # geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
plot(p_mig_ratio)
ggsave(plot=p_mig_ratio, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_ratios.pdf", width = 9, height=7)


# p_samples <- ggplot(samples, aes(x=sample, y=infected)) +
#   geom_point()+
#   # scale_x_log10()+
#   # scale_y_log10()+
#   geom_smooth(method="lm")+
#   stat_cor(aes(label = ..r.label..),  p.accuracy = 0.01, r.accuracy = 0.01)+
#   facet_grid(sampling~.)+
#   # coord_cartesian(xlim=c(0.01,10), ylim=c(0.01,10))+
#   ylab("estimated") + xlab("true") + ggtitle("migration rates") + theme_minimal() +
#   scale_color_OkabeIto()+
#   scale_fill_OkabeIto()
# plot(p_samples)








p_migration <- ggplot(mig[mig$random_ro==1,], aes(x=to/from, y=median/true, ymax=upper/true, ymin=lower/true, 
                                                  color=method,fill=method)) +
  geom_errorbar(alpha=0.25)+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(method="loess")+
  # stat_cor(aes(label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  coord_cartesian(xlim=c(0.01,100), ylim=c(0.01,100))+
  facet_grid(samplingname~nosamples)+
  ylab("median estimate over true value") + xlab("R0 in sink over R0 in source") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))
# add the ratio of migration events per tip to the plot
  # geom_text(data=text.frame[text.frame$random_ro==1,], aes(x=1, y=0.1, label=paste("migration events per tip:", round(ratio, 2))), color="black", size=3, hjust=0)
plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_randomR0_Rcorr.pdf", width = 6, height=3.5)





p_migration <- ggplot(corr_mig[corr_mig$random_ro==0,], aes(x=true, y=median, ymax=upper, ymin=lower, 
                                                  color=method,fill=method)) +
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_errorbar(alpha=0.25)+
  geom_point()+
  geom_smooth(method="lm")+
  stat_cor(aes(label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+

  facet_grid(samplingname~nosamples)+
  coord_cartesian(xlim=c(0.05,250), ylim=c(0.05,250))+
  ylab("estimated migration rate") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699"))
# add the ratio of migration events per tip to the plot
  # geom_text(data=text.frame[text.frame$random_ro==0,], aes(x=1, y=0.1, label=paste("migration events per tip:", round(ratio, 2))), color="black", size=3, hjust=0)
plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_correction.pdf", width = 9, height=7)




corr_mig[corr_mig$run>600,"samplingname"] = "same"
p_migration <- ggplot(corr_mig[corr_mig$run>500,], aes(x=true, y=median, ymax=upper, ymin=lower, 
                                                            color=method,fill=method)) +
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_errorbar(alpha=0.25)+
  geom_point()+
  geom_smooth(method="lm")+
  stat_cor(aes(label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  
  facet_grid(samplingname~nosamples)+
  coord_cartesian(xlim=c(0.05,250), ylim=c(0.05,250))+
  ylab("estimated migration rate") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699"))
# add the ratio of migration events per tip to the plot
# geom_text(data=text.frame[text.frame$random_ro==0,], aes(x=1, y=0.1, label=paste("migration events per tip:", round(ratio, 2))), color="black", size=3, hjust=0)
plot(p_migration)


