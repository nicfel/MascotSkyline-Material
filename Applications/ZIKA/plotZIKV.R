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
require(ggtree)
library(ggpubr)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get the names of all SISCO first (of three) run log files
log <- list.files(path="./out", pattern="*.log", full.names = TRUE)

locations = c("Polynesia", "SouthAmerica", "CentralAmerica", "Caribbean", "Brazil_North", "Brazil_Northeast", "Brazil_Southeast")
# cols = brewer.pal(n = 7, name = 'Dark2')

cols = c("Polynesia"="#1B9E77", "SouthAmerica"="#D95F02", "CentralAmerica"="#7570B3", "Caribbean"="#E7298A", "Brazil North"="#66A61E", "Brazil Northeast"="#E6AB02", "Brazil Southeast"="#A6761D")


# Read in the logs
t1 <- read.table(log[[1]], header=TRUE, sep="\t")
t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
t2 <- read.table(log[[2]], header=TRUE, sep="\t")
t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
t3 <- read.table(log[[3]], header=TRUE, sep="\t")
t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]

t.constant = rbind(t1,t2,t3)


t1 <- read.table(log[[4]], header=TRUE, sep="\t")
t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
t2 <- read.table(log[[5]], header=TRUE, sep="\t")
t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
t3 <- read.table(log[[6]], header=TRUE, sep="\t")
t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]

t.skygrid = rbind(t1,t2,t3)

dat = data.frame()
for (i in seq(1,length(t.skygrid))){
  if (grepl("Ne_", labels(t.skygrid)[[2]][[i]])){
    loc = strsplit(gsub("Ne_", "", labels(t.skygrid)[[2]][[i]]), split="\\.")[[1]]
    hpd.5 = HPDinterval(as.mcmc(exp(t.skygrid[, i])), prob=0.5)
    hpd.95 = HPDinterval(as.mcmc(exp(t.skygrid[, i])), prob=0.95)
    timestart = as.Date("2016-10-12")-0.025*(as.numeric(loc[[2]])-1)*365
    timeend= as.Date("2016-10-12")-0.025*(as.numeric(loc[[2]])-1)*365
    
    
    dat = rbind(dat, data.frame(time = timestart, 
                                l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                location=gsub("_", " ", loc[[1]]),
                                method="skygrid"))
    dat = rbind(dat, data.frame(time = timeend, 
                                l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                location=gsub("_", " ", loc[[1]]),
                                method="skygrid"))
  }
}

# for (i in seq(1,length(t.constant))){
#   if (grepl("Ne.", labels(t.constant)[[2]][[i]])){
#     loc = strsplit(gsub("Ne.", "", labels(t.constant)[[2]][[i]]), split="\\.")[[1]]
#     hpd.5 = HPDinterval(as.mcmc(t.constant[, i]), prob=0.5)
#     hpd.95 = HPDinterval(as.mcmc(t.constant[, i]), prob=0.95)
#     timestart = as.Date("2016-10-12")
#     timeend= as.Date("2016-10-12")-0.025*(158)*365
#     
#     
#     dat = rbind(dat, data.frame(time = timestart, 
#                                 l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
#                                 l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
#                                 location=gsub("_", " ", loc[[1]]),
#                                 method="constant"))
#     dat = rbind(dat, data.frame(time = timeend, 
#                                 l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
#                                 l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
#                                 location=gsub("_", " ", loc[[1]]),
#                                 method="constant"))
#     
#   }
# }


p_ne1 <- ggplot(dat[which(dat$location=="Polynesia" | dat$location=="SouthAmerica" | dat$location=="CentralAmerica"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  # facet_wrap(.~method, ncol=3)+
  xlab("")+
  coord_cartesian(ylim=c(0,75))+
  theme(legend.position = c(.3, .8))
plot(p_ne1)

p_ne2 <- ggplot(dat[which(dat$location=="Brazil North" | dat$location=="Brazil Northeast"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  # facet_wrap(.~method, ncol=3)+
  xlab("")+
  coord_cartesian(ylim=c(0,75))+
  theme(legend.position = c(.3, .8))
plot(p_ne2)

p_ne3 <- ggplot(dat[which(dat$location=="Caribbean" | dat$location=="Brazil Southeast"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  # facet_wrap(.~method, ncol=3)+
  xlab("")+
  coord_cartesian(ylim=c(0,75))+
  theme(legend.position = c(.3, .8))
plot(p_ne3)

# ggsave(plot=p_ne, filename = "../../../MascotSkyline-Text/Figures/zikv_ne.pdf", height=15, width=15)



system(" /Applications/BEAST\\ 2.6.3/bin/logcombiner -b 10 -log out/zikv_skygrid_rep0.trees out/zikv_skygrid_rep1.trees out/zikv_skygrid_rep2.trees -o out/zikv_skygrid.trees")
system(" /Applications/BEAST\\ 2.6.3/bin/logcombiner -b 10 -log out/zikv_constant_rep0.trees out/zikv_constant_rep1.trees out/zikv_constant_rep2.trees -o out/zikv_constant.trees")


system(" /Applications/BEAST\\ 2.6.3/bin/treeannotator -b 0 out/zikv_skygrid.trees out/zikv_skygrid.tree")
system(" /Applications/BEAST\\ 2.6.3/bin/treeannotator -b 0 out/zikv_constant.trees out/zikv_constant.tree")

tree = read.beast("./out/zikv_skygrid.tree")

tree@data$location = c()
tree@data$entropy = c()

for (i in seq(1,length(tree@data$Brazil_North))){
  prob_vals = c()
  entr = 0
  for (j in seq(1,length(locations))){
    prob_vals[[j]] = as.numeric(tree@data[i, locations[[j]]])
    entr = entr + prob_vals[[j]]*log(prob_vals[[j]])
  }
  tree@data$location[[i]] = gsub("_", " ", locations[[which.max(prob_vals)]])
  tree@data$entropy[[i]] = -entr
}

p_skygrid = ggtree(tree, aes(color=location, alpha=entropy)) + geom_tree() + theme_tree() +
  scale_color_manual(values=cols)+
  theme(legend.position = "none") + coord_flip() + scale_x_reverse()

plot(p_skygrid)


tree = read.beast("./out/zikv_constant.tree")

tree@data$location = c()
tree@data$entropy = c()

for (i in seq(1,length(tree@data$Brazil_North))){
  prob_vals = c()
  entr = 0
  for (j in seq(1,length(locations))){
    prob_vals[[j]] = as.numeric(tree@data[i, locations[[j]]])
    entr = entr + prob_vals[[j]]*log(prob_vals[[j]])
  }
  tree@data$location[[i]] = gsub("_", " ", locations[[which.max(prob_vals)]])
  tree@data$entropy[[i]] = -entr
}

p_constant = ggtree(tree, aes(color=location, alpha=entropy)) + geom_tree() + theme_tree() +
  scale_color_manual(values=cols)+
  theme(legend.position = "none") + coord_flip() + scale_x_reverse()

plot(p_constant)


library("gridExtra")

p1 = ggarrange(p_ne1, p_ne2, p_ne3, ncol = 3, labels = c("A", "B", "C"))
p2 = ggarrange(p_skygrid, p_constant, ncol = 2, labels = c("D", "E"))

p3 = ggarrange(p1, p2, ncol = 1)

plot(p3)


ggsave(plot=p3, filename = "../../../MascotSkyline-Text/Figures/zikv.pdf", height=6, width=10)
