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

cols = c("camel"="#E6A777", "human"="#448C82","c"="#E6A777", "h"="#448C82")


mascotlog <- list.files(path="./out", pattern="mascot_.*.trees$", full.names = TRUE)
mascotlog = mascotlog[!grepl("it", mascotlog)]
mascotlog = mascotlog[!grepl("events", mascotlog)]
data = data.frame()
for (i in seq(1, length(mascotlog), 1)){
  print(mascotlog[[i]])
  # make a new tree file tmp.tree that is exactly the same as mascotlog[[i]],
  # but ads END; as the final line
  system(sprintf("cp %s tmp.tree", mascotlog[[i]]))
  system("echo 'END;' >> tmp.tree")
  # read in the trees using ggtree
  mascot_tree = read.beast("tmp.tree")
  
  tmp = strsplit(mascotlog[[i]], split="\\.")[[1]][[3]]
  ind = (as.numeric(tmp)-10)/30+1
  # read in the corresponding dta tree of type dta_mers_red.tmp.trees
  dta_tree = read.beast(paste("./out/dta_mers_red.", tmp, ".trees", sep=""))
  
  # get the height and length of all trees in mascot_tree
  lengths = c()
  height = c()
  for (j in seq(round(length(dta_tree)/10), length(dta_tree))){
    lengths = c(lengths, sum(dta_tree[[j]]@phylo$edge.length))
    height = c(height, max(node.depth.edgelength(dta_tree[[j]]@phylo)))
  }
  # compute HPD
  hpd_l = HPDinterval(as.mcmc(lengths))
  hpd_h = HPDinterval(as.mcmc(height))
  data = rbind(data, data.frame(median=median(height), 
                                lower=hpd_h[1,"lower"], upper=hpd_h[1,"upper"],
                                run=ind, method="DTA", stat="height"))
  data = rbind(data, data.frame(median=median(lengths),
                                lower=hpd_l[1,"lower"], upper=hpd_l[1,"upper"],
                                run=ind, method="DTA", stat="length"))
  # add mascot
  lengths = c()
  height = c()
  for (j in seq(round(length(mascot_tree)/10), length(mascot_tree))){
    lengths = c(lengths, sum(mascot_tree[[j]]@phylo$edge.length))
    height = c(height, max(node.depth.edgelength(mascot_tree[[j]]@phylo)))
  }
  # compute HPD
  hpd_l = HPDinterval(as.mcmc(lengths))
  hpd_h = HPDinterval(as.mcmc(height))
  data = rbind(data, data.frame(median=median(height), 
                                lower=hpd_h[1,"lower"], upper=hpd_h[1,"upper"],
                                run=ind, method="MASCOT-Skyline", stat="height"))
  data = rbind(data, data.frame(median=median(lengths),
                                lower=hpd_l[1,"lower"], upper=hpd_l[1,"upper"],
                                run=ind, method="MASCOT-Skyline", stat="length"))
}


headers = data.frame(number=seq(1,7), lables=c("10% of camel\n100% of human\nsamples", "40% of camel\n100% of human\nsamples", "70% of camel\n100% of human\nsamples", "complete dataset", "100% of camel\n70% of human\nsamples", "100% of camel\n40% of human\nsamples", "100% of camel\n10% of human\nsamples"));


# plot the dta and mascot tree heights and lengths slightly offset using run as the x-axis
p1 = ggplot(data[data$stat=="height", ], aes(x=run, y=median, ymin=lower, ymax=upper, 
                                             color=method)) +
  geom_point(position=position_dodge(width=0.2)) +
  geom_errorbar(position=position_dodge(width=0.2), width=0.1) +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_minimal() +
  theme(legend.position="top") +
  scale_x_continuous(breaks=seq(1,7,1), labels=headers$lables) +
  labs(x="Run", y="Height in years", color="Method", linetype="Statistic") +
  theme(legend.position = "none") +
  xlab("")
plot(p1)

p2 = ggplot(data[data$stat=="length", ], aes(x=run, y=median, ymin=lower, ymax=upper, 
                                             color=method)) +
  geom_point(position=position_dodge(width=0.2)) +
  geom_errorbar(position=position_dodge(width=0.2), width=0.1) +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_minimal() +
  theme(legend.position="top") +
  scale_x_continuous(breaks=seq(1,7,1), labels=headers$lables) +
  labs(x="Run", y="Length in years", color="Method", linetype="Statistic") +
  theme(legend.position = c(0.5, 0.25)) +
  xlab("")
plot(p2)

p = ggarrange(p1, p2, ncol=1)

ggsave(filename="../../../MascotSkyline-Text/Figures/mers_treecomp.pdf", plot=p, height=6.5, width=10)



