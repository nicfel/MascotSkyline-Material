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

# get the names of all SISCO first (of three) run log files
log <- list.files(path="../out", pattern="BD.*1mascot.log", full.names = TRUE)

# define the number of states in the params file
states <- 3

# initialize the data frames
data = data.frame()
mig = data.frame()
incidence = data.frame()
samples = data.frame()
truemig = data.frame()

# get the full trees to later compute the number of migration events
truetrees <- list.files(path="./master", pattern="*.single.tree", full.names = TRUE)

# initialize the plots
plots = list()
pl = 1

# loop over the trees files
for (tr in seq(1, length(truetrees))){
  print(tr) # print the number of the tree file
  trees = read.beast(truetrees[[tr]]) # read in the tree file
  edges = trees@phylo$edge # get the edges of the tree
  nodes = as.numeric(trees@data$node) # get the nodes of the tree
  
  tmp = strsplit(truetrees[[tr]], split="_") # get the run number
  tmp = strsplit(tmp[[1]][2], split="S") # get the run number  
  runnumber = as.numeric(tmp[[1]][2]) # get the run number

  if (runnumber>50){ # if the run number is larger than 50, it is a higher migration run
    sampling="lower migration" # set the sampling to higher migration
    if (runnumber<56){ 
      plots[[runnumber-45]] = ggtree(trees, aes(color=location), size=0.25) + scale_color_OkabeIto() +theme(legend.position = "none") # plot the tree
    }
  }else{
    sampling="higher migration"
    if (runnumber<6){
      plots[[runnumber]] = ggtree(trees, aes(color=location), size=0.25) + scale_color_OkabeIto() +theme(legend.position = "none")
    }
  }

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
  truemig = rbind(truemig, data.frame(events=events, tips = length(trees@phylo$tip.label), sampling=sampling, run=runnumber)) # add the number of migration events to the data frame
}
  
p.trees1 = ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]], ncol = 5, labels = c("F", "G", "H", "I", "J"))
p.trees2 = ggarrange(plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]], ncol = 5, labels = c("F", "G", "H", "I", "J"))

p.trees = ggarrange(p.trees1, p.trees2, ncol=1)

plot(p.trees)
ggsave(plot=p.trees, file="../../../MascotSkyline-Text/Figures/BD_trees.pdf", width = 9, height=6)



support = read.table("nodesupport.tsv", header=T, sep="\t") # read in the node support file
support$ratio=c() # initialize the ratio of migration events per tip

uni_runs = unique(support$run) # get the unique run numbers
for (i in seq(1, length(uni_runs))){ # loop over the unique run numbers
  ind = which(truemig$run==uni_runs[[i]]) # get the index of the run number in the truemig data frame
  ind2 = which(support$run==uni_runs[[i]]) # get the index of the run number in the support data frame  
  support[ind2, "ratio"] = truemig[ind, "events"]/truemig[ind, "tips"] # add the ratio of migration events per tip to the support data frame
}

support$sampling="higher migration"
support[support$run>50, "sampling"]="lower migration"

p.density = ggplot(support, aes(x=support, fill=method, color=method)) + geom_histogram(alpha=0.4, position = 'identity') + 
  theme_minimal() + scale_y_log10() + facet_grid(sampling~.)+
  scale_fill_OkabeIto() +  scale_color_OkabeIto() +xlab("posterior support for true node location") + ylab("count")

plot(p.density)
ggsave(plot=p.density, file="../../../MascotSkyline-Text/Figures/BD_nodesupport.pdf", width = 9, height=6)

# plot the mean node support per run on the y-axis to the ratio of migration events per tip on the x-axis
p.support = ggplot(support, aes(x=ratio, y=support, color=method)) + 
  scale_x_log10() +facet_grid(sampling~.)+theme_minimal() + geom_smooth() +
  scale_color_OkabeIto() +xlab("migration events per tip") + ylab("posterior support for true node location")
  

plot(p.support)


# = ggplot(support, aes(x=ratio, y=support)) + 
#   geom_point(alpha=0.1) +scale_x_log10() +facet_grid(method~.)+theme_minimal()



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
  if (length(t$posterior)>20){
    ess <- effectiveSize(t)
  
    if (min(ess[2:3])<20){
      print(sprintf("ESS value is %f for file %s",min(ess[2:3]),filename1))
    }else{
        dfname <- data.frame(filename = filename1)
        
        # get the number of the run
        tmp = strsplit(filename1, split="_")
        tmp = strsplit(tmp[[1]][2], split="S")
        
        runnumber = as.numeric(tmp[[1]][2])
        
        require("rjson")
        
        if (runnumber>50){
          sampling="lower migration"
        }else{
          sampling="higher migration"
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




log <- list.files(path="../out", pattern="DTABD.*1dta.log", full.names = TRUE)

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
      sampling="lower migration"
    }else{
      sampling="higher migration"
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




p_migration <- ggplot(mig[sample(1:nrow(mig)), ]) +
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


ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/BD_migrationrates.pdf", width = 9, height=6)





corr_comp = data.frame()

for (m in unique(mig$method)){
  for (r in unique(mig[mig$method==m, "run"])){
    ind1 = which(mig$method==m & mig$run==r)
    c = cor(mig[ind1, "true"], mig[ind1, "median"])
    corr_comp = rbind(corr_comp, data.frame(c=c, method=m, run=r, events=truemig[truemig$run==r, "events"]))

  }
}


p_migration <- ggplot(corr_comp, aes(x=events, y=c)) +
  geom_point(aes(color=method))+
  # geom_errorbar(aes(x=true, ymin=lower, ymax=upper, color=method), alpha=0.25)+
  scale_x_log10()+
  # scale_y_log10()+
  # geom_abline(color="red") +
  geom_smooth(aes(color=method, fill=method))+
  stat_cor(aes(olor=method, fill=method, label = ..r.label..), label.x = -1.5, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_grid(facet~.)+
  # coord_cartesian(xlim=c(0.01,100), ylim=c(0.01,100))+
  ylab("estimated") + xlab("true") + theme_minimal() +
  scale_color_OkabeIto()+
  scale_fill_OkabeIto()
plot(p_migration)


ggsave(plot=p_migration, file="../../../MascotSkyline-Text/Figures/BD_correlation.pdf", width = 9, height=6)


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




