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

# define the number of states in the params file
states <- 2

# initialize the data frames
incidence = data.frame()
truemig = data.frame()
totmigdir = data.frame()

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

  events = 0; # initialize the number of migration events
  # initialize vector of 0 and 2 states
  events_dir = rep(0, 2)
  for (i in seq(1, length(edges[,1]))){ # loop over all the edges
    parent_edge = which(edges[,2]==edges[i,1]) # get the parent edge

    if (length(parent_edge)==1){ # if there is a parent edge
      from = which(nodes==edges[i,1]) # get the node number of the parent
      to = which(nodes==edges[i,2]) # get the node number of the child
      if (trees@data$location[from]!=trees@data$location[to]){ # if the parent and child are in   different locations
        events_dir[as.numeric(trees@data$location[from])+1] = events_dir[as.numeric(trees@data$location[from])+1]+1
        events=events+1 # increase the number of migration events
      }
    }
  }
  truemig = rbind(truemig, 
                  data.frame(
                    events=events,
                    events_ratio=(events_dir[1]+1)/(events_dir[2]+1),
                             tips = length(trees@phylo$tip.label), 
                             run=runnumber, 
                             rootTime=min(as.numeric(trees@data$time)))) # add the number of migration events to the data frame

}


log <- list.files(path="../out", pattern="SIR.*1mascot.log", full.names = TRUE)

data = data.frame()
mig = data.frame()
ratio=data.frame()
incidence=data.frame()

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  # Read in the SISCO *.logs
  t1 <- read.table(filename1, header=TRUE, sep="\t")
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]

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
        
        last_time = 0
        
        jsondat <- fromJSON(file=gsub("1mascot.log","master.json", gsub("../out/","./master/",log[[i]])))
        tot_incidence = c()
        times = jsondat$t
        
        for (j in seq(1, length(jsondat$I))){
          if (runnumber%%100<6 & runnumber%%100>0){
            incidence = rbind(incidence, data.frame(time=times, y=jsondat$I[[j]], 
                                                  state=paste("state", j-1, sep=""),
                                                  run=runnumber))
          }
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
                                            run=runnumber))
              
            }
          }
        }       
        
        # loop over the Ne0's
        for (j in seq(7,8)){
          val_mean = mean(t[,param_labels[j]])
          val_median = median(t[,param_labels[j]])
          hpd = HPDinterval(as.mcmc(t[,param_labels[j]]))
          value_name = strsplit(param_labels[j], split="\\.")[[1]][[1]]
          true=params[param_index, param_labels[j]]

          from_cp = tot_incidence[[as.numeric(gsub("state","",strsplit(param_labels[j], split="[._]")[[1]][[3]])) +1]]
          to_cp = tot_incidence[[as.numeric(gsub("state","",strsplit(param_labels[j], split="[._]")[[1]][[5]])) +1]]
          
          mig = rbind(mig, data.frame(true=true, mean=val_mean, median=val_median, 
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                          from_state=strsplit(param_labels[j], split="[._]")[[1]][[3]],
                                          to_state = strsplit(param_labels[j], split="[._]")[[1]][[5]],
                                          from_cp = from_cp,
                                          to_cp = to_cp,
                                          run=runnumber, method="MASCOT-Skyline"))
        }

        
        ratio_mean = mean(t[,param_labels[7]]/t[,param_labels[8]])
        ratio_median = median(t[,param_labels[7]]/t[,param_labels[8]])
        hpd = HPDinterval(as.mcmc(t[,param_labels[7]]/t[,param_labels[8]]))
        value_name = strsplit(param_labels[7], split="\\.")[[1]][[1]]
        true=params[param_index, param_labels[7]]/params[param_index, param_labels[8]]

        from_cp = tot_incidence[[as.numeric(gsub("state","",strsplit(param_labels[7], split="[._]")[[1]][[3]])) +1]]
        to_cp = tot_incidence[[as.numeric(gsub("state","",strsplit(param_labels[8], split="[._]")[[1]][[5]])) +1]]

        mean_vals = (t[,param_labels[7]]+t[,param_labels[8]])/2
        mean_hpd = HPDinterval(as.mcmc(mean_vals))
        
        before = length(ratio$true)
        
        ratio = rbind(ratio, data.frame(true=true, mean=ratio_mean, median=ratio_median, 
                                        lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                        from_state=strsplit(param_labels[7], split="[._]")[[1]][[3]],
                                        to_state = strsplit(param_labels[8], split="[._]")[[1]][[5]],
                                        from_cp = from_cp,
                                        to_cp = to_cp,
                                        mean_mean = mean(mean_vals),
                                        mean_median = median(mean_vals),
                                        mean_lower = mean_hpd[1,"lower"],
                                        mean_upper = mean_hpd[1,"upper"],
                                        true_mean = (params[param_index, param_labels[7]]+params[param_index, param_labels[8]])/2,
                                        run=runnumber, method="MASCOT-Skyline"))


        if (length(ratio$true)!=before+1){
          das
        }
    }
  }
}



log <- list.files(path="../out", pattern="DTASIR.*1dta.log", full.names = TRUE)
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
  
  if (min(ess[2:3])<100){
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
    
    label_7 = gsub("f_migrationRatesSkyline", "geo.rates", param_labels[7])
    label_7 = gsub("_to_", ".", label_7)
    label_8 = gsub("f_migrationRatesSkyline", "geo.rates", param_labels[8])
    label_8 = gsub("_to_", ".", label_8)
    

    # loop over the Ne0's
    for (j in seq(7,8)){
      label = gsub("f_migrationRatesSkyline", "geo.rates", param_labels[j])
      label = gsub("_to_", ".", label)
      
      norm_correction = (t[,label_7]+t[,label_8])/2

      val_mean = mean(t[,label]*t[,"geo.meanRate"]/norm_correction)
      val_median = median(t[,label]*t[,"geo.meanRate"]/norm_correction)
      hpd = HPDinterval(as.mcmc(t[,label]*t[,"geo.meanRate"]/norm_correction))
      true=params[param_index, param_labels[j]]
      from = paste("R0",strsplit(param_labels[j], split="[._]")[[1]][[3]], sep=".")
      to = paste("R0",strsplit(param_labels[j], split="[._]")[[1]][[5]], sep=".")
      from_to_cp = mig[mig$run==runnumber &
                       mig$from_state==strsplit(param_labels[j], split="[._]")[[1]][[3]] &
                       mig$to_state==strsplit(param_labels[j], split="[._]")[[1]][[5]],
                        c("from_cp", "to_cp")]
      if (nrow(from_to_cp)==0){
        from_to_cp[1,c(1,2)]=NA
      }
        
      mig = rbind(mig, data.frame(true=true, mean=val_mean, median=val_median, 
                                  lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                  from_state=strsplit(param_labels[j], split="[._]")[[1]][[3]],
                                  to_state = strsplit(param_labels[j], split="[._]")[[1]][[5]],
                                  from_cp = from_to_cp[[1]],
                                  to_cp = from_to_cp[[2]],
                                  run=runnumber, method="DTA"))
    }
    ratio_mean = mean(t[,label_7]/t[,label_8])
    ratio_median = median(t[,label_7]/t[,label_8])
    hpd = HPDinterval(as.mcmc(t[,label_7]/t[,label_8]))
    true=params[param_index, param_labels[7]]/params[param_index, param_labels[8]]
    from = paste("R0",strsplit(param_labels[7], split="[._]")[[1]][[3]], sep=".")
    to = paste("R0",strsplit(param_labels[7], split="[._]")[[1]][[5]], sep=".")
    from_to_cp = mig[mig$run==runnumber &
                     mig$from_state==strsplit(param_labels[7], split="[._]")[[1]][[3]] &
                     mig$to_state==strsplit(param_labels[7], split="[._]")[[1]][[5]],
                      c("from_cp", "to_cp")]

    if (nrow(from_to_cp)==0){
      from_to_cp[1,c(1,2)]=NA
    }

    mean_vals = t[,"geo.meanRate"]
    mean_hpd = HPDinterval(as.mcmc(mean_vals))

    ratio = rbind(ratio, data.frame(true=true, mean=ratio_mean, median=ratio_median, 
                                    lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                    from_state=strsplit(param_labels[7], split="[._]")[[1]][[3]],
                                    to_state = strsplit(param_labels[7], split="[._]")[[1]][[5]],
                                    from_cp = from_to_cp[[1]],
                                    to_cp = from_to_cp[[2]],
                                    mean_mean = mean(mean_vals),
                                    mean_median = median(mean_vals),
                                    mean_lower = mean_hpd[1,"lower"],
                                    mean_upper = mean_hpd[1,"upper"],
                                    true_mean = (params[param_index, param_labels[7]]+params[param_index, param_labels[8]])/2,
                                    run=runnumber, method="DTA"))


  }
}

require("ggpubr")

mig[mig$run>0 & mig$run<101, "name"] = "low migration\n250 samples"
mig[mig$run>100 & mig$run<201, "name"] = "low migration\n500 samples"
mig[mig$run>200 & mig$run<301, "name"] = "high migration\n250 samples"
mig[mig$run>300 & mig$run<401, "name"] = "high migration\n500 samples"
mig[mig$run>400 & mig$run<501, "name"] = "low migration\nrandom R0\n250 samples"
mig[mig$run>500 & mig$run<601, "name"] = "high migration\nrandom R0\n250 samples"
mig[mig$run>600 & mig$run<701, "name"] = "low migration\neven sampling rate"
mig[mig$run>700 & mig$run<801, "name"] = "high migration\neven sampling rate"
mig[mig$run>800 & mig$run<901, "name"] = "low migration\nconstant sampling\n250 samples"
mig[mig$run>900 & mig$run<1001, "name"] = "high migration\nconstant sampling\n250 samples"

# # reorder the name facets to have the above three facets come last
mig$name = factor(mig$name, levels=c("low migration\n250 samples", "low migration\n500 samples", 
                                     "high migration\n250 samples","high migration\n500 samples",
                                     "low migration\nrandom R0\n250 samples", "high migration\nrandom R0\n250 samples", 
                                     "low migration\neven sampling rate", "high migration\neven sampling rate",
                                     "low migration\nconstant sampling\n250 samples", "high migration\nconstant sampling\n250 samples"))

# make  a dataframe with the correlation coefficients and the HPD values make the y value method dependet
corr_migration = data.frame()
for (m in unique(mig$method)){
  for (n in unique(mig$name)){
    tmp = mig[mig$method==m & mig$name==n,]
    coverage = sum(tmp$lower<tmp$true & tmp$upper>tmp$true)/length(tmp$true)
    correlation_coeff = cor(log(tmp$median), log(tmp$true))
    if (m=="DTA"){
      y=1
    }else{
      y=0
    }
    corr_migration = rbind(corr_migration, 
                        data.frame(method=m, name=n, y=y,
                                   coverage=coverage, 
                                   correlation_coeff=correlation_coeff))
  }
}
corr_migration$name = factor(corr_migration$name, levels(mig$name))

p_migration <- ggplot(mig) +
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
  geom_point(aes(x=true, y=median, color=method))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(aes(x=true, y=median, color=method, fill=method), method="lm")+
  geom_text(data=corr_migration, aes(x=10, y=exp(y*1.7-10),color=method,
                                     label=paste("R=",round(correlation_coeff,2), 
                                                 ", ", "cov=", round(coverage,2), sep="")), size=4)+
  # stat_cor(aes(x=true, y=median, color=method, fill=method, label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_wrap(name~., ncol=5)+
  # coord_cartesian(xlim=c(0.05,250), ylim=c(0.05,250))+
  ylab("estimated migration rate") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  theme(legend.position = "top")

plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates.pdf", width = 12, height=9)


set.seed(123)  # Set a seed for reproducibility
ratio <- ratio[sample(nrow(ratio)), ]
ratio$name =""
ratio[ratio$run>0 & ratio$run<101, "name"] = "low migration\n250 samples"
ratio[ratio$run>100 & ratio$run<201, "name"] = "low migration\n500 samples"
ratio[ratio$run>200 & ratio$run<301, "name"] = "high migration\n250 samples"
ratio[ratio$run>300 & ratio$run<401, "name"] = "high migration\n500 samples"
ratio[ratio$run>400 & ratio$run<501, "name"] = "low migration\nrandom R0\n250 samples"
ratio[ratio$run>500 & ratio$run<601, "name"] = "high migration\nrandom R0\n250 samples"
ratio[ratio$run>600 & ratio$run<701, "name"] = "low migration\neven sampling rate"
ratio[ratio$run>700 & ratio$run<801, "name"] = "high migration\neven sampling rate"
ratio[ratio$run>800 & ratio$run<901, "name"] = "low migration\nconstant sampling\n250 samples"
ratio[ratio$run>900 & ratio$run<1001, "name"] = "high migration\nconstant sampling\n250 samples"

# # reorder the name facets to have the above three facets come last
ratio$name = factor(ratio$name, levels=c("low migration\n250 samples", "low migration\n500 samples", 
                                     "high migration\n250 samples","high migration\n500 samples",
                                     "low migration\nrandom R0\n250 samples", "high migration\nrandom R0\n250 samples", 
                                     "low migration\neven sampling rate", "high migration\neven sampling rate",
                                     "low migration\nconstant sampling\n250 samples", "high migration\nconstant sampling\n250 samples"))

# p_migration <- ggplot(ratio) +
#   geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
#   geom_point(aes(x=true, y=median, color=method))+
#   scale_x_log10()+
#   scale_y_log10()+
#   geom_abline(color="red") +
#   geom_smooth(aes(x=true, y=median, color=method, fill=method), method="lm")+
#   stat_cor(aes(x=true, y=median, color=method, fill=method, label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
#   facet_wrap(name~., ncol=4)+
#   # coord_cartesian(xlim=c(0.001,100), ylim=c(0.001,100))+
#   ylab("estimated migration rate") + xlab("simulated migration rate") + theme_minimal() +
#   scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
#   scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))
# 
# 
# plot(p_migration)
# ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationratesratios.pdf", width = 9, height=6)


p_migration <- ggplot(mig,aes(x=true, y=(upper-lower)/median, color=method, fill=method)) +
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(yintercept=5.285402, color="red")+
  geom_smooth()+
  facet_wrap(name~., ncol=5)+
  # coord_cartesian(xlim=c(0.05,250))+
  ylab("HPD width divided by median estimate") + xlab("simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  theme(legend.position = "top")
plot(p_migration)
ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_cv.pdf", width = 12, height=9)

corr_migration = data.frame()
for (m in unique(ratio$method)){
  for (n in unique(ratio$name)){
    tmp = ratio[ratio$method==m & ratio$name==n,]
    coverage = sum(tmp$mean_lower<tmp$true_mean & tmp$mean_upper>tmp$true_mean)/length(tmp$true_mean)
    correlation_coeff = cor(log(tmp$true_mean), log(tmp$mean_median), method = "pearson")
    if (m=="DTA"){
      y=1
    }else{
      y=0
    }
    corr_migration = rbind(corr_migration, 
                           data.frame(method=m, name=n, y=y,
                                      coverage=coverage, 
                                      correlation_coeff=correlation_coeff))
  }
}
corr_migration$name = factor(corr_migration$name, levels(mig$name))

p_mean <- ggplot(ratio, aes(x=true_mean, y=mean_median, color=method, fill=method)) +
  geom_errorbar(aes(x=true_mean, ymin=mean_lower, ymax=mean_upper), alpha=0.25)+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(method="lm")+
  geom_text(data=corr_migration, aes(x=20, y=exp(y*1.2-3),color=method,
                                     label=paste("R=",round(correlation_coeff,2), 
                                                 ", ", "cov=", round(coverage,2), sep="")), size=4)+
  
  # stat_cor(aes(x=true_mean, y=mean_median, color=method, fill=method, label = ..r.label..), 
  # label.x = 0, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_wrap(name~., ncol=5)+
  # coord_cartesian(xlim=c(0.001,100), ylim=c(0.001,100))+
  ylab("mean estimated migration rate") + xlab("mean simulated migration rate") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  theme(legend.position="top")
plot(p_mean)
ggsave(plot=p_mean, file="../../../MascotSkyline-Text/Figures/SIR_migrationrates_mean.pdf", width = 12, height=9)



corr_mig = ratio[ratio$method=="DTA",]
corr_mig$method = "DTA with case correction"
corr_mig$median = corr_mig$median * corr_mig$to_cp/corr_mig$from_cp
corr_mig$upper = corr_mig$upper * corr_mig$to_cp/corr_mig$from_cp
corr_mig$lower = corr_mig$lower * corr_mig$to_cp/corr_mig$from_cp
# remove all lines that contain Na's
corr_mig = corr_mig[!is.na(corr_mig$median),]


corr_mig = rbind(corr_mig, ratio)



set.seed(123)  # Set a seed for reproducibility
corr_mig <- corr_mig[sample(nrow(corr_mig)), ]

corr_migration = data.frame()
for (m in unique(corr_mig$method)){
  for (n in unique(corr_mig$name)){
    tmp = corr_mig[corr_mig$method==m & corr_mig$name==n,]
    coverage = sum(tmp$lower<tmp$true & tmp$upper>tmp$true)/length(tmp$true)
    correlation_coeff = cor(log(tmp$median), log(tmp$true))
    if (m=="DTA"){
      y=1
    }else if (m=="DTA with case correction"){
      y=0
    }else{
      y=-1
    }
    corr_migration = rbind(corr_migration, 
                           data.frame(method=m, name=n, y=y,
                                      coverage=coverage, 
                                      correlation_coeff=correlation_coeff))
  }
}
corr_migration$name = factor(corr_migration$name, levels(mig$name))


p_corrmig <- ggplot(corr_mig, aes(x=true,
                                    color=method,fill=method)) +
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
  geom_point(aes(x=true, y=median, color=method))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(aes(x=true, y=median), method="lm")+
  # stat_cor(aes(label = ..r.label.., y=median), label.x = -2, p.accuracy = 0.01, r.accuracy = 0.01)+
  geom_text(data=corr_migration, aes(x=10, y=exp(y*2-9),color=method,
                                     label=paste("R=",round(correlation_coeff,2),
                                                 ", ", "cov=", round(coverage,2), sep="")), size=4)+
  facet_wrap(name~., ncol=5)+
  # coord_cartesian(xlim=c(0.01,100), ylim=c(0.01,100))+
  ylab("estimated migration rate ratio") + xlab("simulated migration rate ratio") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510", "DTA with case correction"="#336699")) +
  theme(legend.position = "top")
plot(p_corrmig)

ggsave(plot=p_corrmig, file="../../../MascotSkyline-Text/Figures/SIR_migrationratesratio.pdf", width = 12, height=9)

# make a new data frame with the ratio of migration rates for DTA and MASCOT-Skyline and the 
# ratio of migration events.
ratio_mig = data.frame()
for (a in unique(ratio$run)){
  for (b in unique(ratio$method)){
    # get the corresponding ratio of migration events
    mig_row = ratio[which(ratio$run==a & ratio$method==b),]
    if (length(mig_row$true==1)){
      ratio_mig = rbind(ratio_mig, 
        data.frame(runnumber = a,
          method = b,
          ratio=truemig[truemig$run==a, "events_ratio"],
          median = mig_row$median,
          lower = mig_row$lower,
          upper = mig_row$upper,
          name = mig_row$name
        ))
    }
  }
}

p_events = ggplot(data=ratio_mig, aes(x=ratio, y=median, ymax=upper, ymin=lower, color=method, fill=method)) +
  geom_point()+
  geom_errorbar(alpha=0.25)+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(method="lm")+
  # stat_cor(aes(x=ratio, y=median, color=method, fill=method, label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_wrap(name~., ncol=4)+
  # coord_cartesian(xlim=c(0.001,100), ylim=c(0.001,100))+
  ylab("estimated migration rate ratio") + xlab("true event ratio") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510")) +
  theme(legend.position="top")


plot(p_events)
ggsave(plot=p_events, file="../../../MascotSkyline-Text/Figures/SIR_migrationevents.pdf", width = 12, height = 9)




data[data$run>0 & data$run<101, "name"] = "low migration\n250 samples"
data[data$run>100 & data$run<201, "name"] = "low migration\n500 samples"
data[data$run>200 & data$run<301, "name"] = "high migration\n250 samples"
data[data$run>300 & data$run<401, "name"] = "high migration\n500 samples"
data[data$run>400 & data$run<501, "name"] = "low migration\nrandom R0\n250 samples"
data[data$run>500 & data$run<601, "name"] = "high migration\nrandom R0\n250 samples"
data[data$run>600 & data$run<701, "name"] = "low migration\neven sampling rate"
data[data$run>700 & data$run<801, "name"] = "high migration\neven sampling rate"
data[data$run>800 & data$run<901, "name"] = "low migration\nconstant sampling"
data[data$run>900 & data$run<1001, "name"] = "high migration\nconstant sampling"
#take the modulo of the run to denote left to right in the plot later on
data$runplot = data$run%%100

# do the same for incidence
incidence[incidence$run>0 & incidence$run<101, "name"] = "low migration\n250 samples"
incidence[incidence$run>100 & incidence$run<201, "name"] = "low migration\n500 samples"
incidence[incidence$run>200 & incidence$run<301, "name"] = "high migration\n250 samples"
incidence[incidence$run>300 & incidence$run<401, "name"] = "high migration\n500 samples"
incidence[incidence$run>400 & incidence$run<501, "name"] = "low migration\nrandom R0\n250 samples"
incidence[incidence$run>500 & incidence$run<601, "name"] = "high migration\nrandom R0\n250 samples"
incidence[incidence$run>600 & incidence$run<701, "name"] = "low migration\neven sampling rate"
incidence[incidence$run>700 & incidence$run<801, "name"] = "high migration\neven sampling rate"
incidence[incidence$run>800 & incidence$run<901, "name"] = "low migration\nconstant sampling"
incidence[incidence$run>900 & incidence$run<1001, "name"] = "high migration\nconstant sampling"
incidence$runplot = incidence$run%%100

subset1 = data[data$run<501,]

p_skyline <- ggplot(data=subset1, aes(x=time, color=state, fill=state, group=state))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
  geom_line(aes(y=median,linetype = "median estimate"))+
  geom_line(data=incidence[is.element(incidence$run, unique(subset1$run)),], aes(x=time, y=log(y/100), linetype = "true value", color=state))+
  facet_grid(name~runplot) +
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
ggsave(plot=p_skyline, file="../../../MascotSkyline-Text/Figures/SIR_Ne_prevalence1.pdf", width = 15, height=10)

# do the same for all runs above 500
subset2 = data[data$run>500,]
p_skyline <- ggplot(data=subset2, aes(x=time, color=state, fill=state, group=state))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
  geom_line(aes(y=median,linetype = "median estimate"))+
  geom_line(data=incidence[is.element(incidence$run, unique(subset2$run)),], aes(x=time, y=log(y/100), linetype = "true value", color=state))+
  facet_grid(name~runplot) +
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
ggsave(plot=p_skyline, file="../../../MascotSkyline-Text/Figures/SIR_Ne_prevalence2.pdf", width = 15, height=10)


# compute the number of migraiton event per sample for the group of runs above
for (a in levels(mig$name)){
  b = mig[mig$name==a,"run"]
  avg_events = mean(truemig[which(truemig$run %in% b), "events"]/truemig[which(truemig$run %in% b), "tips"])
  avg_tips = mean(truemig[which(truemig$run %in% b), "tips"])
  print(sprintf("The average number of migration events per sample for %s is %f is %f", a, avg_events, avg_tips))
}




support = read.table("nodesupport.tsv", header=T, sep="\t") # read in the node support file
support$ratio=c() # initialize the ratio of migration events per tip

uni_runs = unique(support$run) # get the unique run numbers
for (i in seq(1, length(uni_runs))){ # loop over the unique run numbers
  ind = which(truemig$run==uni_runs[[i]]) # get the index of the run number in the truemig data frame
  ind2 = which(support$run==uni_runs[[i]]) # get the index of the run number in the support data frame
  support[ind2, "ratio"] = truemig[ind, "events"]/truemig[ind, "tips"] # add the ratio of migration events per tip to the support data frame
}
support[support$method=="dta", "method"]="DTA"
support[support$method=="mascot", "method"]="MASCOT-Skyline"

support[support$run>0 & support$run<101, "name"] = "low migration\n250 samples"
support[support$run>100 & support$run<201, "name"] = "low migration\n500 samples"
support[support$run>200 & support$run<301, "name"] = "high migration\n250 samples"
support[support$run>300 & support$run<401, "name"] = "high migration\n500 samples"
support[support$run>400 & support$run<501, "name"] = "low migration\nrandom R0\n250 samples"
support[support$run>500 & support$run<601, "name"] = "high migration\nrandom R0\n250 samples"
support[support$run>600 & support$run<701, "name"] = "low migration\neven sampling rate"
support[support$run>700 & support$run<801, "name"] = "high migration\neven sampling rate"
support[support$run>800 & support$run<901, "name"] = "low migration\nconstant sampling"
support[support$run>900 & support$run<1001, "name"] = "high migration\nconstant sampling"

support$name = factor(support$name, levels=c("low migration\n250 samples", "low migration\n500 samples", 
                                         "high migration\n250 samples","high migration\n500 samples",
                                         "low migration\nrandom R0\n250 samples", "high migration\nrandom R0\n250 samples", 
                                         "low migration\neven sampling rate", "high migration\neven sampling rate",
                                         "low migration\nconstant sampling\n250 samples", "high migration\nconstant sampling\n250 samples"))


p.density = ggplot(support, aes(x=support, fill=method, color=method)) +
  geom_histogram(alpha=0.4, position = 'identity', bins=10) +
  theme_minimal() +
  # scale_y_log10() +
  xlab("posterior support for true node location") + ylab("count")+
  facet_wrap(name~., ncol=5)+
  # coord_cartesian(xlim=c(0.001,100), ylim=c(0.001,100))+
  ylab("density") + xlab("support for true node state") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  theme(legend.position="top")
  
plot(p.density)
ggsave(plot=p.density, file="../../../MascotSkyline-Text/Figures/SIR_nodesupport.pdf", width = 12, height=9)



p.ratio = ggplot(support, aes(x=ratio, fill=method, color=method)) +
  geom_histogram(alpha=0.4, position = 'identity', bins=10) +
  theme_minimal() +
  scale_x_log10() +
  xlab("posterior support for true node location") + ylab("count")+
  facet_wrap(name~., ncol=5)+
  # coord_cartesian(xlim=c(0.001,100), ylim=c(0.001,100))+
  ylab("density") + xlab("number of migration events per sampled tip") + theme_minimal() +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  theme(legend.position="top")
plot(p.ratio)

