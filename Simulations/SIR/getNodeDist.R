getNodeDist <- function(file) {
  # load the libraries
  library(ggplot2)
  library(coda)
  library("methods")
  library(colorblindr)
  library(ggtree)
  library(treeio)
  library(ggpubr)
  require(ape)
  require(phytools)

  print(file)
  # read the tree as a beast tree
  require(treeio)
  trees <- read.beast(file)
  
  # also read in the tree in ./master that has the same S%d in its name
  S = strsplit(file, split="_")[[1]][2]
  master = paste0("./master/SIR_", S, "_master.tree")
  master.tree <- read.beast(master)
  
  # get the height of every node in master.tree
  master.heights <- master.tree@phylo$edge.length
  
  master.tree@data$location = paste("state", master.tree@data$location, sep="")
  
  difference = c()
  # loop over the trees in trees
  for (tree in trees) {
    # get the height of every node in the current tree
    tree.heights <- tree@phylo$edge.length
    diff = 0
    # loop through each node in the tree
    for (node in 1:length(tree.heights)) {
      # find the closest node in the master tree based on height
      closest_node <- which.min(abs(master.heights - tree.heights[node]))
      # print(min(abs(master.heights - tree.heights[node])))
      # Use closest_node as the corresponding node in master.tree
      # Perform your intended operations with closest_node here
      closest_node_nr = which(master.tree@data$node==master.tree@phylo$edge[closest_node,2])
      node_nr = which(tree@data$node==tree@phylo$edge[node,2])
      
      if (master.tree@data$location[closest_node_nr]!=tree@data$typeTrait[node_nr]){
        if (master.tree@phylo$edge[closest_node,2]<=500){
          print(master.tree@phylo$tip.label[closest_node_nr])
          print(master.tree@data$location[closest_node_nr])
          print(tree@data$typeTrait[node_nr])
        }
      }
      diff = diff + as.numeric(master.tree@data$location[closest_node_nr]!=tree@data$typeTrait[node_nr])
    }
    difference = c(difference, diff)
  }
  # make a new log file in the treeDist with the differences, the file is \t delimited and the header is difference
  write.table(difference, paste0("./treeDist/SIR", S, ".mascot.log"), sep="\t", row.names=FALSE, col.names=FALSE)

  # now, read it the corresponding dta tree file in ../out with the filename of type DTASIR_S%d_1dta.trees
  # first, check if the file exists
  if (!file.exists(paste0("../out/DTASIR_", S, "_1dta.trees"))){
    print(paste0("../out/DTASIR_", S, "_1dta.trees"))
    next
  }
  dta = paste0("../out/DTASIR_", S, "_1dta.trees")
  dta.trees <- read.beast(dta)
  difference = c()
  # loop over the trees in trees
  for (tree in dta.trees){
    # get the height of every node in the current tree
    tree.heights <- tree@phylo$edge.length
    diff = 0
    # loop through each node in the tree
    for (node in 1:length(tree.heights)) {
      # find the closest node in the master tree based on height
      closest_node <- which.min(abs(master.heights - tree.heights[node]))
      # print(min(abs(master.heights - tree.heights[node])))
      # Use closest_node as the corresponding node in master.tree
      # Perform your intended operations with closest_node here
      closest_node_nr = which(master.tree@data$node==master.tree@phylo$edge[closest_node,2])
      node_nr = which(tree@data$node==tree@phylo$edge[node,2])
      
      if (master.tree@data$location[closest_node_nr]!=tree@data$geo[node_nr]){
        if (master.tree@phylo$edge[closest_node,2]<=500){
          print(master.tree@phylo$tip.label[closest_node_nr])
          print(master.tree@data$location[closest_node_nr])
          print(tree@data$geo[node_nr])
        }
      }
      diff = diff + as.numeric(master.tree@data$location[closest_node_nr]!=tree@data$geo[node_nr])
    }
    difference = c(difference, diff)
  }
  # make a new log file in the treeDist with the differences, the file is \t delimited and the header is difference
  write.table(difference, paste0("./treeDist/SIR", S, ".dta.log"), sep="\t", row.names=FALSE, col.names=FALSE)
}
