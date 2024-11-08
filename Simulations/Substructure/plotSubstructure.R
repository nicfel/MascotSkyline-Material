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

log <- list.files(path="../out", pattern="Substructure.*1mascot.log", full.names = TRUE)

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

  events = c(0, 0); # initialize the number of migration events
  background_events = 0
  for (i in seq(1, length(edges[,1]))){ # loop over all the edges
    parent_edge = which(edges[,2]==edges[i,1]) # get the parent edge
    if (length(parent_edge)==1){ # if there is a parent edge
      from = which(nodes==edges[i,1]) # get the node number of the parent
      to = which(nodes==edges[i,2]) # get the node number of the child
      # state 0 is any location <= 19, state 1 ant location > 19
      l1 = as.numeric(as.numeric(trees@data$location[from]) > 19)
      l2 = as.numeric(as.numeric(trees@data$location[to]) > 19)
      if (l1!=l2){ # if the parent and child are in different locations
        events[l1+1]=events[l1+1]+1 # increase the number of migration events
      }
      
      if (l1==0){
        ln1 = as.numeric(as.numeric(trees@data$location[from]) ==0)
        ln2 = as.numeric(as.numeric(trees@data$location[to]) > 0)
        if (ln1!=ln2){
          background_events = background_events+1
        }
        
      }
      
    }
  }
  samples = c(0, 0);
  for (i in seq(1,length(trees@phylo$tip.label)+1)){
    index = which(edges[,2]==i)
    loc = as.numeric(as.numeric(trees@data$location[index]) > 19)+1
    samples[loc] = samples[loc]+1
  }
  coal_events = c(0, 0)
  for (i in seq(1,length(edges))){
    index = which(edges[,1]==i)
    if (length(index)==2){
      from = which(nodes==edges[index[1],2])
      loc = as.numeric(as.numeric(trees@data$location[from]) > 19)+1
      coal_events[loc] = coal_events[loc]+1
    }
  }

  truemig = rbind(truemig, 
                  data.frame(migrationEvents.state0_to_state1=events[1],
                             migrationEvents.state1_to_state0=events[2],
                             tips = length(trees@phylo$tip.label), 
                             run=runnumber, 
                             samples0=samples[1],
                             samples1=samples[2],
                             coal_sink = coal_events[2],
                             background_events = background_events,
                             rootTime=min(as.numeric(trees@data$time)))) # add the number of migration events to the data frame
}

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")
  # Read in the SISCO *.logs
  t1 <- read.table(filename1, header=TRUE, sep="\t")
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
  t=t1
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
        
        # jsondat <- fromJSON(file=gsub("1mascot.log","master.json", gsub("../out/","./master/",log[[i]])))
        # tot_incidence = c()
        # times = jsondat$t
        # for (j in seq(1, length(jsondat$I))){
        #   incidence = rbind(incidence, data.frame(time=times, y=jsondat$I[[j]],
        #                                           state=paste("state", j-1, sep=""),
        #                                           run=runnumber))
        #   tot_incidence = c(tot_incidence, sum(diff(jsondat$I[[j]])==1))
        # }
        # 
        
        param_index = which(params$run==runnumber)
        param_index = param_index[length(param_index)]
        
        last_time = truemig[which(truemig$run==runnumber), "rootTime"]
        
        
        param_labels = labels(truemig)[[2]]
        
        # loop over the migration events
        for (j in seq(1,2)){
          mean = mean(t[,param_labels[j]])
          median = mean(t[,param_labels[j]])
          hpd = HPDinterval(as.mcmc(t[,param_labels[j]]))
          true=truemig[truemig$run==runnumber, param_labels[j]]
          str = strsplit(param_labels[j], split="\\.")[[1]][[2]]
          from = strsplit(str, split="_to_")[[1]][1]
          to = strsplit(str, split="_to_")[[1]][2]
          from = gsub("state0", "events from source to sink", from)
          from = gsub("state1", "events sink to source", from)
          mig = rbind(mig, data.frame(true=true, mean=mean, median=median, 
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                          samples_source = truemig[which(truemig$run==runnumber),"samples0"],
                                          samples_sink = truemig[which(truemig$run==runnumber),"samples1"],
                                          coal_sink = truemig[which(truemig$run==runnumber),"coal_sink"],
                                          background_events = truemig[which(truemig$run==runnumber),"background_events"],
                                          from = from,
                                          to = to,
                                          run=runnumber))
        }
        
        if (mig[mig$run==runnumber & mig$from=="events sink to source", "mean"]>100){
          # loop over the Ne0's if the the run number is x01 x02 x03 x04 x05 (take modulo 100)
          for (j in seq(1, 2)){
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
        
    }
  }
}

p_skyline <- ggplot(data=mig, aes(x=coal_sink, y=median, color=from, group=from))+
  geom_errorbar(aes(ymin=lower, ymax=upper), alpha=0.2)+
  geom_point()+
  geom_smooth(aes(x=coal_sink, y=true,
                  color="true number of events (smooth average)"), fill=NA, linetype='dashed')+
  ylab("number of migration events") + 
  xlab("true number of coalescent events observed in the sink") +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  # facet_grid(.~order) +
  theme_minimal()  +
  scale_color_OkabeIto(name="")+
  scale_fill_OkabeIto(name="")+
  theme(legend.position="top")

plot(p_skyline)
ggsave(plot=p_skyline, file="../../../MascotSkyline-Text/Figures/Substructure_events.pdf", width = 6, height=4)


data[data$state=="state0", "state"] = "source"
data[data$state=="state1", "state"] = "sink"

p = ggplot(data, aes(x=time, y=median, ymin=lower, ymax=upper, fill=state, color=state))+
  geom_ribbon(color=NA, alpha=0.25)+
  geom_line()+
  theme_minimal()  +
  facet_wrap(.~run, ncol=5) +
  scale_color_OkabeIto(name="")+
  scale_fill_OkabeIto(name="")+
  theme(legend.position="top")

plot(p)

log <- list.files(path="../out", pattern="Substructure.*1mascot.log", full.names = TRUE)

# get the full trees to later compute the number of migration events
dtatrees <- list.files(path="./dtaout", pattern="*.trees", full.names = TRUE)
dtamig = data.frame()
# loop over the trees files
for (tr in seq(1, length(dtatrees))){
  print(dtatrees[[tr]]) # print the number of the tree file
  trees = read.beast(dtatrees[[tr]]) # read in the tree file
  tmp = strsplit(dtatrees[[tr]], split="_") # get the run number
  tmp = strsplit(tmp[[1]][2], split="S") # get the run number  
  runnumber = as.numeric(tmp[[1]][2]) # get the run number
  events = list(); # initialize the number of migration events
  c =1
  for (k in seq(50, length(trees), 5)){
    events[[c]] = c(0, 0)
    edges = trees[[k]]@phylo$edge # get the edges of the tree
    nodes = as.numeric(trees[[k]]@data$node) # get the nodes of the tree
    background_events = 0
    for (i in seq(1, length(edges[,1]))){ # loop over all the edges
      parent_edge = which(edges[,2]==edges[i,1]) # get the parent edge
      if (length(parent_edge)==1){ # if there is a parent edge
        from = which(nodes==edges[i,1]) # get the node number of the parent
        to = which(nodes==edges[i,2]) # get the node number of the child
        # state 0 is any location <= 19, state 1 ant location > 19
        l1 = as.numeric(trees[[k]]@data$geo[from]=="state1")
        l2 = as.numeric(trees[[k]]@data$geo[to]=="state1")
        if (l1!=l2){ # if the parent and child are in different locations
          events[[c]][l1+1]=events[[c]][l1+1]+1 # increase the number of migration events
        }
      }
    }
    c=c+1
  }
  # convert the list to a matrix
  events_matrix = do.call(rbind, events)
  
  from_name = c("events from source to sink", "events sink to source")
  for (k in seq(1, length(from_name))){
    mean = mean(events_matrix[,k])
    hpd = HPDinterval(as.mcmc(events_matrix[,k]))
    median = median(events_matrix[,k])
    true = 1
    to=""
    dtamig = rbind(dtamig, data.frame(true=true, mean=mean, median=median, 
                                lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                samples_source = truemig[which(truemig$run==runnumber),"samples0"],
                                samples_sink = truemig[which(truemig$run==runnumber),"samples1"],
                                coal_sink = truemig[which(truemig$run==runnumber),"coal_sink"],
                                background_events = truemig[which(truemig$run==runnumber),"background_events"],
                                from = from_name[k],
                                to = to,
                                run=runnumber))
  }
  
}

p_skyline <- ggplot(data=dtamig, aes(x=samples_sink, y=median, color=from, fill=from, group=from))+
  geom_errorbar(aes(ymin=lower, ymax=upper), alpha=0.2)+
  geom_point()+
  ylab("number of migration events") + 
  xlab("samples in sink population") +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  # facet_grid(.~order) +
  theme_minimal()  +
  scale_color_OkabeIto(name="")+
  scale_fill_OkabeIto(name="")+
  theme(legend.position="top")

plot(p_skyline)

ggsave(plot=p_skyline, file="../../../MascotSkyline-Text/Figures/SubstructureDTA_events.pdf", width = 6, height=4)
