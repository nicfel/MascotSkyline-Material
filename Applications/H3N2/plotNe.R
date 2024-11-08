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

locations = c("Oceania", "Other")

mrsi.structured = as.Date("2005-12-30")
mrsi.unstructured = as.Date("2005-09-02")


dat = data.frame()
mig = data.frame()
for (j in seq (1,length(log))){
  t1 <- read.table(log[[j]], header=TRUE, sep="\t")
  t1 <- t1[-seq(1,ceiling(length(t1$Sample)/2)), ]
  t2 <- read.table(gsub("rep0","rep1",log[[j]]), header=TRUE, sep="\t")
  t2 <- t2[-seq(1,ceiling(length(t2$Sample)/2)), ]
  t3 <- read.table(gsub("rep0","rep2",log[[j]]), header=TRUE, sep="\t")
  t3 <- t3[-seq(1,ceiling(length(t3$Sample)/2)), ]
  t = rbind(t1,t2,t3)

  sampling = strsplit(log[j], split="_")[[1]][[2]]
  outsidesamples = as.numeric(strsplit(log[j], split="_")[[1]][[3]])

    # print(sampling)
  if (outsidesamples==50){
    mrsi = mrsi.structured
  }else{
    mrsi = mrsi.unstructured
  }
  struct = sampling
  print(struct)
  print(mrsi)
  for (i in seq(1,length(t))){
    if (grepl("Ne_", labels(t)[[2]][[i]])){
      loc = strsplit(gsub("Ne_", "", labels(t)[[2]][[i]]), split="\\.")[[1]]
      hpd.5 = HPDinterval(as.mcmc(exp(t[, i])), prob=0.5)
      hpd.95 = HPDinterval(as.mcmc(exp(t[, i])), prob=0.95)
      timestart = mrsi-0.01*(as.numeric(loc[[2]])-1)*365
      timeend = mrsi-0.01*(as.numeric(loc[[2]])-1)*365
      dat = rbind(dat, data.frame(time = timestart, 
                                  l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                  l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                  location=gsub("_", " ", loc[[1]]),
                                  sampling=struct,
                                  outsidesamples=outsidesamples))
      dat = rbind(dat, data.frame(time = timeend, 
                                  l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                  l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                  location=gsub("_", " ", loc[[1]]),
                                  sampling=struct,
                                  outsidesamples=outsidesamples))
    }else if (grepl("migrationEvents.", labels(t)[[2]][[i]])){
      loc = gsub("_"," ",gsub("migrationEvents.", "", labels(t)[[2]][[i]]))
      
      hpd.5 = HPDinterval(as.mcmc(t[, i]), prob=0.5)
      hpd.95 = HPDinterval(as.mcmc(t[, i]), prob=0.95)
      mig = rbind(mig, data.frame(l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                  l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                  sampling=struct,
                                  outsidesamples=outsidesamples,
                                  location=loc))
    }
  }
}

dat$sampling <- relevel(dat$sampling, "structured")

p_ne1 <- ggplot(dat[which(dat$location=="Oceania"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=interaction(sampling,location)), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=interaction(sampling,location)), alpha=0.75) +
  theme_minimal()+
  scale_fill_OkabeIto(name="",order = c(1,2), labels=c("Ne in Oceania when\naccounting for structure","Ne in Oceania when\nnot accounting for structure"))+
  scale_y_continuous(breaks=seq(0,15,5))+
  # facet_wrap(.~outsidesamples, ncol=2)+
  coord_cartesian(ylim=c(0,20))+
  scale_x_date(limits=c(as.Date("2000-01-01"), as.Date("2006-01-01")))+
  # facet_wrap(.~location, ncol=3)+
  xlab("")+ 
  ylab("Ne")+ 
  theme(legend.position = c(0.5, 0.7))

dat2 = rbind(dat[which(dat$location=="Oceania" & dat$sampling=="unstructured"),], 
             dat[which(dat$location=="Other" & dat$sampling=="structured"),])
dat2$sampling <- relevel(dat2$sampling, "structured")
dat2$location <- relevel(dat2$location, "Other")

p_ne2 <- ggplot(dat2) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=interaction(sampling, location)), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=interaction(sampling, location)), alpha=0.75) +
  theme_minimal()+
  scale_fill_OkabeIto(name="",order = c(2,3), labels=c("Ne outside Oceania when\naccounting for structure","Ne in Oceania when\nnot accounting for structure"))+
  scale_y_continuous(breaks=seq(0,15,5))+
  # facet_wrap(.~outsidesamples, ncol=2)+
  coord_cartesian(ylim=c(0,20))+
  scale_x_date(limits=c(as.Date("2000-01-01"), as.Date("2006-01-01")))+
  # facet_wrap(.~location, ncol=3)+
  xlab("") +
  ylab("Ne")+ 
  theme(legend.position = c(0.5, 0.7))
plot(p_ne2)

system(" /Applications/BEAST\\ 2.6.3/bin/logcombiner -b 20 -log out/h3n2_structured_0_rep0.trees out/h3n2_structured_0_rep1.trees out/h3n2_structured_0_rep2.trees -o out/h3n2.trees")
system(" /Applications/BEAST\\ 2.6.3/bin/treeannotator -b 0 out/h3n2.trees out/h3n2.tree")

require(ggtree)
library(treeio)

tree = read.beast("./out/h3n2.tree")

tree@data$location = list()
tree@data$location2 = c()

tree@data$entropy = list()

for (i in seq(1,length(tree@data$Oceania))){
  prob_vals = c()
  entr = 0
  for (j in seq(1,length(locations))){
    prob_vals[[j]] = as.numeric(tree@data[i, locations[[j]]])
    if (prob_vals[[j]]!=0){
      entr = entr + prob_vals[[j]]*log(prob_vals[[j]])
    }
    
  }

  tree@data$location[[i]] = prob_vals[[1]]
  tree@data$entropy[[i]] = -entr
}

p = ggtree(tree, aes(colour=as.numeric(location)),  mrsd=mrsi.unstructured) + 
  geom_nodepoint(aes(color=as.numeric(location)), size=2)+
  geom_tippoint(aes(color=as.numeric(location)), size=1)+
  theme_tree() +
  theme_minimal()+
  scale_color_continuous(name="Probability of\nbeing in Oceania\n",high="#E69F00", low="#009E73") +
  theme(legend.position = "left") +
  coord_flip() +
  scale_x_reverse()+
  scale_y_discrete(breaks=c()) +
  theme(legend.position = c(0.8, 0.7))
plot(p)
# ggsave(plot=p, filename = "../../../MascotSkyline-Text/Figures/zikv_tree.pdf", height=6, width=6)


library("gridExtra")

p2 = ggarrange(p_ne1, p_ne2, ncol = 1, labels = c("B", "C"))
p3 = ggarrange(p,p2, ncol = 2, labels = c("A", ""))

plot(p3)

ggsave(plot=p3, filename = "../../../MascotSkyline-Text/Figures/h3n2.pdf", height=6, width=12)

