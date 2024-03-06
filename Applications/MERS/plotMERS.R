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
library(ggtree)
library(treeio)
library(gridExtra)
library(ggpubr)



# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get the names of all SISCO first (of three) run log files
log <- list.files(path="./out", pattern=".*.trees", full.names = TRUE)
log = log[!grepl("it", log)]

for ( i in seq(1, length(log))){
  if (!grepl("events", log[[i]])){
    print(log[[i]])
    # system(paste("/Applications/BEAST\\ 2.7.4/bin/treeannotator -b 10",
    #              log[[i]], gsub(".trees", ".tree", log[[i]])),ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
}

cols = c("camel"="#E6A777", "human"="#448C82","c"="#E6A777", "h"="#448C82")


mascotlog <- list.files(path="./out", pattern="mascot_.*.log$", full.names = TRUE)
mascotlog = mascotlog[!grepl("it", mascotlog)]

headers = data.frame(number=seq(1,7), lables=c("10% of camel\n100% of human\nsamples", "40% of camel\n100% of human\nsamples", "70% of camel\n100% of human\nsamples", "complete dataset", "100% of camel\n70% of human\nsamples", "100% of camel\n40% of human\nsamples", "100% of camel\n10% of human\nsamples"));


data = data.frame()

for (i in seq(1, length(mascotlog), 1)){
  t <- read.table(mascotlog[[i]], header=TRUE, sep="\t")
  
  t = t[-seq(1, length(t$Sample)/10), ]
  
  tmp = strsplit(mascotlog[[i]], split="\\.")[[1]][[3]]
  ind = (as.numeric(tmp)-10)/30+1
  
  for (j in seq(1,11)){
    labels = paste("SkylineNe.camel.", j, sep="")
    mean = mean(t[,labels])
    median = mean(t[,labels])
    hpd = HPDinterval(as.mcmc(t[,labels]))
    value_name = j
    
    data = rbind(data, data.frame(time=j, mean=mean, median=median, 
                                lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                run=ind, method="MASCOT-Skyline", state="camel"))
  }
  
  for (j in seq(1,11)){
    labels = paste("SkylineNe.human.", j, sep="")
    mean = mean(t[,labels])
    median = mean(t[,labels])
    hpd = HPDinterval(as.mcmc(t[,labels]))
    value_name = j
    
    data = rbind(data, data.frame(time=j, mean=mean, median=median, 
                                  lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                  run=ind, method="MASCOT-Skyline", state="human"))
  }
  

}



mascot <- list.files(path="./out", pattern="mascot_.*.tree$", full.names = TRUE)
mascot = mascot[!grepl("it", mascot)]

locations = c("camel", "human")

pml = list()
titels = list()
trends = list()

for (i in seq(1, length(mascot))){
  tree = read.beast(mascot[[i]])
  
  tree@data$location = c()
  tree@data$entropy = c()
  
  for (k in seq(1,length(tree@data$camel))){
    prob_vals = c()
    entr = 0
    for (j in seq(1,length(locations))){
      prob_vals[[j]] = as.numeric(tree@data[k, locations[[j]]])
      if (prob_vals[[j]]>0){
        entr = entr + prob_vals[[j]]*log(prob_vals[[j]])
      }
    }
    tree@data$location[k] = gsub("_", " ", locations[[which.max(prob_vals)]])
    tree@data$entropy[k] = entr
  }  
  tmp = strsplit(mascot[[i]], split="\\.")[[1]][[3]]
  ind = (as.numeric(tmp)-10)/30+1  
  
  pml[[ind]] = ggtree(tree, aes(color=location)) + 
    # geom_nodepoint(aes(color=location),size=1)+
    geom_tippoint(aes(color=location),size=0.25)+
    scale_color_manual(values=cols)+
    # scale_alpha(range=c(0.2,1))+
    theme_minimal()+
    theme(legend.position = "none") +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())

  # make a figure with only the header label and otherwise completely empty
  titels[[ind]] = ggplot(data=headers[headers$number==ind, ]) + 
    geom_text(aes(x=0, y=0, label=lables), size=4)+
    theme_void()+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  trends[[ind]] = ggplot(data=data[which(data$run==ind),], aes(x=time, color=state, fill=state, group=state))+
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
    geom_line(aes(y=median,linetype = "median estimate"))+
    # facet_grid(.~run)+
    theme_minimal()+
    theme(legend.position = "none",
          panel.margin = unit(0, "cm"),  # set margin to zero
          panel.spacing = unit(0, "cm")  # set spacing to zero
    )+
    scale_fill_manual(values=cols)+
    scale_color_manual(values=cols)+
    scale_x_reverse()+
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
    
    
}

# create a title plot
p0 = do.call(grid.arrange, c(titels, nrow=1, ncol=7))
# create the plot grid
p1 = do.call(grid.arrange, c(pml, nrow=1, ncol=7))

ggsave(filename="mers_mascot.pdf", plot=p1, height=3, width=20)

dta <- list.files(path="./out", pattern="dta_.*.tree$", full.names = TRUE)
dta <- list.files(path="./out", pattern="dta_.*.tree$", full.names = TRUE)
dta = dta[!grepl("it", dta)]


pdl = list()

for (i in seq(1, length(dta))){
  tree = read.beast(dta[[i]])
  
  for (k in seq(1,length(tree@data$geo))){
    
    tree@data$geo[k] = strsplit(tree@data$geo[k], split="+")[[1]][[1]]
  }
  tmp = strsplit(dta[[i]], split="\\.")[[1]][[3]]
  ind = (as.numeric(tmp)-10)/30+1
  
  pdl[[ind]] = ggtree(tree, aes(color=geo)) + 
    # geom_nodepoint(aes(color=location),size=1)+
    geom_tippoint(aes(color=geo),size=0.25)+
    scale_color_manual(values=cols)+
    # scale_alpha(range=c(0.2,1))+
    theme_minimal()+
    theme(legend.position = "none") +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
}

p2 = do.call(grid.arrange, c(pdl, nrow=1, ncol=7))

ggsave(filename="mers_dta.pdf", plot=p2, height=3, width=20)

p_trends = do.call(grid.arrange, c(trends, nrow=1, ncol=7))

spacer_plot <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank(), panel.grid = element_blank(), panel.border = element_blank())



p_tot = ggarrange(p0, p1, p_trends,spacer_plot, p2, ncol = 1, heights=c(0.32,1,1,0.1,1))


plot(p_tot)

ggsave(filename="../../../MascotSkyline-Text/Figures/mers.pdf", plot=p_tot, height=6.5, width=10)


p_legend = trends[[1]] + theme(legend.position = "bottom")

ggsave(filename="../../../MascotSkyline-Text/Figures/mers_legend.pdf", plot=p_legend, height=1, width=4)



mascotlog <- list.files(path="./out", pattern="mascotnolo_.*.log$", full.names = TRUE)
mascotlog <- list.files(path="./out", pattern="mascotnolo_.*.log$", full.names = TRUE)


data = data.frame()

for (i in seq(1, length(mascotlog), 1)){
  t <- read.table(mascotlog[[i]], header=TRUE, sep="\t")
  
  t = t[-seq(1, length(t$Sample)/10), ]
  
  tmp = strsplit(mascotlog[[i]], split="\\.")[[1]][[3]]
  ind = (as.numeric(tmp)-10)/30+1
  
  for (j in seq(1,11)){
    labels = paste("SkylineNe.camel.", j, sep="")
    mean = mean(t[,labels])
    median = mean(t[,labels])
    hpd = HPDinterval(as.mcmc(t[,labels]))
    value_name = j
    
    data = rbind(data, data.frame(time=j, mean=mean, median=median, 
                                  lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                  run=ind, method="MASCOT-Skyline", state="camel"))
  }
  for (j in seq(1,11)){
    labels = paste("SkylineNe.human.", j, sep="")
    mean = mean(t[,labels])
    median = mean(t[,labels])
    hpd = HPDinterval(as.mcmc(t[,labels]))
    value_name = j
    
    data = rbind(data, data.frame(time=j, mean=mean, median=median, 
                                  lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                  run=ind, method="MASCOT-Skyline", state="human"))
  }
  
  
}


mascot <- list.files(path="./out", pattern="mascotnolo_.*.tree$", full.names = TRUE)

locations = c("camel", "human")

pml = list()

for (i in seq(1, length(mascot))){
  tree = read.beast(mascot[[i]])
  
  tree@data$location = c()
  tree@data$entropy = c()
  
  for (k in seq(1,length(tree@data$camel))){
    prob_vals = c()
    entr = 0
    for (j in seq(1,length(locations))){
      prob_vals[[j]] = as.numeric(tree@data[k, locations[[j]]])
      if (prob_vals[[j]]>0){
        entr = entr + prob_vals[[j]]*log(prob_vals[[j]])
      }
    }
    tree@data$location[k] = gsub("_", " ", locations[[which.max(prob_vals)]])
    tree@data$entropy[k] = entr
  }
  
  tmp = strsplit(mascot[[i]], split="\\.")[[1]][[3]]
  ind = (as.numeric(tmp)-10)/30+1
  
  
  pml[[ind]] = ggtree(tree, aes(color=location)) + 
    # geom_nodepoint(aes(color=location),size=1)+
    geom_tippoint(aes(color=location),size=0.25)+
    scale_color_manual(values=cols)+
    # scale_alpha(range=c(0.2,1))+
    theme_minimal()+
    theme(legend.position = "none") +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
  
  trends[[ind]] = ggplot(data=data[which(data$run==ind),], aes(x=time, color=state, fill=state, group=state))+
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
    geom_line(aes(y=median,linetype = "median estimate"))+
    # facet_grid(.~run)+
    theme_minimal()+
    theme(legend.position = "none",
          panel.margin = unit(0, "cm"),  # set margin to zero
          panel.spacing = unit(0, "cm")  # set spacing to zero
    )+
    scale_fill_manual(values=cols)+
    scale_color_manual(values=cols)+
    scale_x_reverse()+
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
  
}

library(gridExtra)
# create the plot grid
p1 = do.call(grid.arrange, c(pml, nrow=1, ncol=7))

ggsave(filename="mers_mascot.pdf", plot=p1, height=3, width=20)


dta <- list.files(path="./out", pattern="dtanolo_.*.tree$", full.names = TRUE)
dta <- list.files(path="./out", pattern="dtanolo_.*.tree$", full.names = TRUE)


pdl = list()

for (i in seq(1, length(dta))){
  tree = read.beast(dta[[i]])
  
  for (k in seq(1,length(tree@data$geo))){
    
    tree@data$geo[k] = strsplit(tree@data$geo[k], split="+")[[1]][[1]]
  }
  tmp = strsplit(dta[[i]], split="\\.")[[1]][[3]]
  ind = (as.numeric(tmp)-10)/30+1
  
  
  pdl[[ind]] = ggtree(tree, aes(color=geo)) + 
    # geom_nodepoint(aes(color=location),size=1)+
    geom_tippoint(aes(color=geo),size=0.25)+
    scale_color_manual(values=cols)+
    # scale_alpha(range=c(0.2,1))+
    theme_minimal()+
    theme(legend.position = "none") +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
}


ptrends = do.call(grid.arrange, c(trends, nrow=1, ncol=7))

p2 = do.call(grid.arrange, c(pdl, nrow=1, ncol=7))

ggsave(filename="mers_dta.pdf", plot=p2, height=3, width=20)



p_tot = ggarrange(p0, p1,ptrends,spacer_plot,p2, ncol = 1, heights=c(0.32,1,1,0.1,1))

ggsave(filename="../../../MascotSkyline-Text/Figures/mersnolo.pdf", plot=p_tot, height=6.5, width=10)

plot(p_tot)



