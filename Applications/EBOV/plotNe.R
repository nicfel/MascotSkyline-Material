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
log <- list.files(path="./out", pattern="*rep0.log", full.names = TRUE)

locations = c("SLE", "GIN", "LBR")

dat = data.frame()
for (j in seq (1,length(log))){
  t1 <- read.table(log[[j]], header=TRUE, sep="\t")
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
  t2 <- read.table(gsub("rep0","rep1",log[[j]]), header=TRUE, sep="\t")
  t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
  t3 <- read.table(gsub("rep0","rep2",log[[j]]), header=TRUE, sep="\t")
  t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]
  t = rbind(t1,t2,t3)

  sampling =strsplit(log[j], split="_")[[1]][[2]]
  struct = strsplit(log[j], split="_")[[1]][[3]]
  for (i in seq(1,length(t))){
    if (grepl("Ne_", labels(t)[[2]][[i]])){
      loc = strsplit(gsub("logNe_", "", labels(t)[[2]][[i]]), split="\\.")[[1]]
      hpd.5 = HPDinterval(as.mcmc(exp(t[, i])), prob=0.5)
      hpd.95 = HPDinterval(as.mcmc(exp(t[, i])), prob=0.95)
      timestart = as.Date("2015-09-15")-0.01*(as.numeric(loc[[2]])-1)*365
      timeend= as.Date("2015-09-15")-0.01*(as.numeric(loc[[2]])-1)*365
      dat = rbind(dat, data.frame(time = timestart, 
                                  l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                  l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                  location=gsub("_", " ", loc[[1]]),
                                  sampling=sampling, structure=struct))
      dat = rbind(dat, data.frame(time = timeend, 
                                  l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                  l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                  location=gsub("_", " ", loc[[1]]),
                                  sampling=sampling, structure=struct))
      
    }
  }
}

p_ne1 <- ggplot(dat[which(dat$location=="SLE" | dat$location=="GIN" | dat$location=="LBR"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.75) +
  theme_minimal()+
  scale_color_OkabeIto(breaks=c("state0","state1"))+
  scale_fill_OkabeIto()+
  scale_x_date(limits=c(as.Date("2014-01-01"), as.Date("2015-07-01")))+
  # facet_wrap(.~location, ncol=3)+
  xlab("")
plot(p_ne1)


das

p_ne2 <- ggplot(dat[which(dat$location=="Brazil North" | dat$location=="Brazil Northeast"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.5) +
  theme_minimal()+
  scale_color_OkabeIto(breaks=c("state0","state1"))+
  scale_fill_OkabeIto()+
  scale_x_date()+
  # facet_wrap(.~location, ncol=3)+
  xlab("")+
  theme(legend.position = "none")
plot(p_ne2)

p_ne2 <- ggplot(dat[which(dat$location=="Caribbean" | dat$location=="Brazil Southeast"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.5) +
  theme_minimal()+
  scale_color_OkabeIto(breaks=c("state0","state1"))+
  scale_fill_OkabeIto()+
  scale_x_date()+
  # facet_wrap(.~location, ncol=3)+
  xlab("")+
  theme(legend.position = "none")
plot(p_ne2)




ggsave(plot=p_ne, filename = "../../../MascotSkyline-Text/Figures/zikv_ne.pdf", height=15, width=15)

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

p = ggtree(tree, aes(color=location, alpha=entropy)) + geom_tree() + theme_tree() +
  scale_color_OkabeIto() + theme(legend.position = "left") + coord_flip() + scale_x_reverse()

ggsave(plot=p, filename = "../../../MascotSkyline-Text/Figures/zikv_tree.pdf", height=6, width=6)


