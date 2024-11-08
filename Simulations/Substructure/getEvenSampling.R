# load the libraries
library(ggplot2)
library(coda)
library("methods")
library(colorblindr)
library(ggtree)
library(treeio)
library(ggpubr)
library(ape)


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
incidence = data.frame()
truemig = data.frame()
totmigdir = data.frame()

# get the full trees to later compute the number of migration events
truetrees <- list.files(path="./master", pattern="*.tree", full.names = TRUE)
# initialize the plots
plots = list()
pl = 1

# loop over the trees files
for (tr in seq(1, length(truetrees))){
  print(truetrees[[tr]]) # print the number of the tree file
  tmp = strsplit(truetrees[[tr]], split="_") # get the run number
  tmp = strsplit(tmp[[1]][2], split="S") # get the run number  
  runnumber = as.numeric(tmp[[1]][2]) # get the run number
  # if run number is not between 101 and 200 or 301 and 400, skip the run
  if ((runnumber<=800)){
    next
  }
  if (!grepl(".original.", truetrees[[tr]])){
    next
  }

  trees = read.beast(truetrees[[tr]]) # read in the tree file

  # get all the sampling times of the leafs from location 0 and 1
  times0 = c()
  times1 = c()
  index0 = c()
  index1 = c()
  for (i in seq(1, length(trees@data$reaction))){
    # if reaction is sampling
    if (trees@data$reaction[[i]] == "Sampling"){
      # change the corresponding tip label
      trees@phylo$tip.label[as.numeric(trees@data$node[[i]])] = paste("inv", 
                                                          trees@phylo$tip.label[as.numeric(trees@data$node[[i]])],
                                                          "loc_",
                                                          trees@data$location[[i]], 
                                                          "_time_",
                                                          trees@data$time[[i]], 
                                                          sep="")
      # if location is 0
      if (trees@data$location[[i]] == 0){
        times0 = c(times0, as.numeric(trees@data$time[[i]]))
        index0 = c(index0, i)
      }
      # if location is 1
      if (trees@data$location[[i]] == 1){
        times1 = c(times1, as.numeric(trees@data$time[[i]]))
        index1 = c(index1, i)
      }
    }
  }
  # Vector to keep track of the tips to keep
  keep = c()

  # split the times into 20 bins based on the minimum and maximum
  # value in times0
  bins = seq(min(times0), max(times0), length.out=20)
  # until there are 200 samples from times0, select a random bin,
  # then select a random time from that bin if it has not been used yet
  selected0 = c()
  while (length(selected0)<250){
    bin = sample(length(bins)-1, 1)
    tmp = times0[times0>=bins[bin] & times0<bins[bin+1]]
    if (length(tmp)>0){
      tmp1 = sample(tmp, 1)
       if (!any(selected0==tmp1)){
        selected0 = c(selected0, tmp1)
        keep = c(keep, index0[times0==tmp1])
       }
    }
  }
  # do the same for times1
  bins = seq(min(times1), max(times1), length.out=20)
  selected1 = c()
  while (length(selected1)<250){
    bin = sample(length(bins)-1, 1)
    tmp = times1[times1>=bins[bin] & times1<bins[bin+1]]
    if (length(tmp)>0){
      tmp1 = sample(tmp, 1)
      if (!any(selected1==tmp1)){
        selected1 = c(selected1, tmp1)
        keep = c(keep, index1[times1==tmp1])
      }
    }
  }

  # plot the density of both over time
  p = ggplot() +
    geom_density(aes(x=times0, color="0"), data=data.frame(times0)) +
    geom_density(aes(x=times1, color="1"), data=data.frame(times1)) +
    geom_density(aes(x=selected0, color="0"), data=data.frame(selected0)) +
    geom_density(aes(x=selected1, color="1"), data=data.frame(selected1)) +
    
    scale_color_manual(name="Location", values=c("0"="red", "1"="blue")) +
    labs(x="Time", y="Density") +
    theme_bw()
  plot(p)
  
  # Get the names of the tips to be removed
  tips_to_remove <- setdiff(trees@phylo$tip.label, trees@phylo$tip.label[as.numeric(trees@data$node[keep])])
  
  # Prune the tree
  pruned_tree <- drop.tip(trees@phylo, tips_to_remove)
  
  # make a new tree file name by adding 700 to the run number when the runnumber is 
  # between 101 and 200 or by adding 600 when it is between 301 and 400, but 
  newname = paste("./master/SIR_S", runnumber, "_master.tree", sep="")

  # write the tree to file in master folder
  write.tree(pruned_tree, file=newname)

}
  
