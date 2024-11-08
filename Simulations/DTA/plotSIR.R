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
log <- list.files(path="../out", pattern="SIR.*1mascot.log", full.names = TRUE)

states <- 3

data = data.frame()
mig = data.frame()
incidence = data.frame()

samples = data.frame()

# get the number of migration events
truetrees <- list.files(path="./master", pattern="*.single.tree", full.names = TRUE)
truemig = data.frame()

for (tr in seq(1, length(truetrees))){
  print(tr)
  trees = read.beast(truetrees[[tr]])
  edges = trees@phylo$edge
  nodes = as.numeric(trees@data$node)
  
  tmp = strsplit(truetrees[[tr]], split="_")
  tmp = strsplit(tmp[[1]][2], split="S")
  
  runnumber = as.numeric(tmp[[1]][2])

  if (runnumber>50){
    sampling="higher migration"
  }else{
    sampling="lower migration"
  }
  
    
  events = 0;
  for (i in seq(1, length(edges[,1]))){
    parent_edge = which(edges[,2]==edges[i,1])
    
    if (length(parent_edge)==1){
      from = which(nodes==edges[i,1])
      to = which(nodes==edges[i,2])
      if (trees@data$location[from]!=trees@data$location[to]){
        events=events+1
      }
    }
  }
  truemig = rbind(truemig, data.frame(events=events, tips = length(trees@phylo$tip.label), sampling=sampling, run=runnumber))
}
  



# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
# for (i in seq(1,1,1)){
    
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
  if (length(t$posterior>20)){
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
        
        if (runnumber>50){
          sampling="higher migration"
        }else{
          sampling="lower migration"
        }
        
        last_time = 0;
        
        jsondat <- fromJSON(file=gsub("1mascot.log","master.json", gsub("../out/","./master/",log[[i]])))
        times = jsondat$t[seq(1,length(jsondat$t), 10)]
        for (j in seq(1, length(jsondat$I))){
          ind = which(diff(jsondat$sample[[j]])==1)
          samples = rbind(samples, data.frame(run=runnumber,  
                                              state=paste("state", j-1, sep=""), 
                                              sample=length(ind),infected =sum(diff(jsondat$I[[j]])==1),
                                              sampling=sampling
                                              ))
          if (length(ind)>0){
            last_time = max(last_time, jsondat$t[ind[length(ind)] +1])
          }
          incidence = rbind(incidence, data.frame(time=times, y=jsondat$I[[j]][seq(1,length(jsondat$t), 10)], 
                                                  state=paste("state", j-1, sep=""),
                                                  sampling=sampling,
                                                  run=runnumber))
        }
        
  
        
        param_index = which(params$run==runnumber)
        
        param_index = param_index[length(param_index)]
      
        # loop over the Ne0's
        for (j in seq(1, length(jsondat$I))){
          for (k in seq(1, 26)){
            name = paste("SkylineNe.state", j-1, ".", k, sep="")
            mean = mean(t[,name])
            median = mean(t[,name])
            hpd = HPDinterval(as.mcmc(t[,name]))
            value_name = strsplit(name, split="\\.")[[1]][[1]]
            data = rbind(data, data.frame(time=last_time - (k-1)/25*t$Tree.height[1], name = value_name, mean=mean, median=median, 
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                          state=paste("state", j-1, sep=""),
                                          sampling=sampling,
                                          run=runnumber))
            
          }
        }
        
        
        # loop over the Ne0's
        for (j in seq(2,length(param_labels))){
          mean = mean(t[,param_labels[j]])
          median = mean(t[,param_labels[j]])
          hpd = HPDinterval(as.mcmc(t[,param_labels[j]]))
          value_name = strsplit(param_labels[j], split="\\.")[[1]][[1]]
          true=params[param_index, param_labels[j]]
          
          mig = rbind(mig, data.frame(true=true, mean=mean, median=median, 
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                          sampling=sampling,
                                          run=runnumber, method="MASCOT-Skyline"))
        }
        
  
    }
  }
}




# 
p_skyline <- ggplot(data=data[data$name=="SkylineNe" & data$run<20,], aes(x=time, color=state, fill=state, group=state))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
  geom_line(aes(y=median,linetype = "median estimate"))+
  geom_line(data=incidence[incidence$run<20,], aes(x=time, y=log(y/50), linetype = "true value", color=state))+
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
  coord_cartesian(ylim=c(-4,8))

plot(p_skyline)
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
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/2)), ]
  # t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
  # t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]
  
  t = t1
  # t = rbind(t1, t2, t3)
  # calculate ess values
  ess <- effectiveSize(t)
  
  if (min(ess[2:3])<40){
    print(sprintf("ESS value is %f for file %s",min(ess[2:3]),filename1))
  }else{
    dfname <- data.frame(filename = filename1)
    
    # get the number of the run
    tmp = strsplit(filename1, split="_")
    tmp = strsplit(tmp[[1]][2], split="S")
    
    runnumber = as.numeric(tmp[[1]][2])

    param_index = which(params$run==runnumber)
    param_index = param_index[length(param_index)]
    
    if (runnumber>50){
      sampling="higher migration"
    }else{
      sampling="lower migration"
    }
    
    # loop over the Ne0's
    for (j in seq(2,length(param_labels))){
      label = gsub("f_migrationRatesSkyline", "geo.rates", param_labels[j])
      label = gsub("_to_", ".", label)
      mean = mean(t[,label]*t[,"geo.meanRate"])
      median = mean(t[,label]*t[,"geo.meanRate"])
      hpd = HPDinterval(as.mcmc(t[,label]*t[,"geo.meanRate"]))
      true=params[param_index, param_labels[j]]
      
      mig = rbind(mig, data.frame(true=true, mean=mean, median=median, 
                                  lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                  sampling = sampling,
                                  run=runnumber, method="DTA"))
    }
  }
}


sampls = unique(mig$sampling)
for (i in seq(1, length(sampls))){
  ratio = mean(truemig[truemig$sampling==sampls[[i]], "events" ]/truemig[truemig$sampling==sampls[[i]], "tips" ])
  
  mig[mig$sampling==sampls[[i]], "facet"] = paste(sampls[[i]], "\n",ratio, "\nmigration events per tip", sep="")
  print(ratio)
}


require("ggpubr")
p_migration <- ggplot(mig) +
  geom_point(aes(x=true, y=median, color=method))+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(color="red") +
  geom_smooth(aes(x=true, y=median, color=method, fill=method), method="lm")+
  stat_cor(aes(x=true, y=median, color=method, fill=method, label = ..r.label..), label.x = -1.5, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_grid(facet~.)+
  coord_cartesian(xlim=c(0.01,100), ylim=c(0.01,100))+
  ylab("estimated") + xlab("true") + theme_minimal() +
  scale_color_OkabeIto()+
  scale_fill_OkabeIto()
plot(p_migration)


ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/migrationrates.pdf", width = 9, height=6)



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




