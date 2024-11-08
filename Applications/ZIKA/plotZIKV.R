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
library(treeio)
library(plyr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get the names of all SISCO first (of three) run log files
log <- list.files(path="./out", pattern="*.log", full.names = TRUE)

locations = c("Polynesia", "SouthAmerica", "CentralAmerica", "Caribbean", "Brazil_North", "Brazil_Northeast", "Brazil_Southeast")
# cols = brewer.pal(n = 7, name = 'Dark2')

cols = c("Polynesia"="#A4A9AA", "SouthAmerica"="#03A381", "CentralAmerica"="#B39B31", "Caribbean"="#F3756D", "Brazil North"="#4A578B", "Brazil Northeast"="#5594C6", "Brazil Southeast"="#D19BC4")


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

param.labels = labels(t.skygrid)[[2]]
tmp = strsplit(param.labels[grepl("SkylineNe", param.labels)], split="\\.")
max_val = 0
for (i in seq(1, length(tmp))){
  if (length(tmp[[i]])==3){
    max_val = max(max_val, as.numeric(tmp[[i]][[3]]))
  }
}


dat = data.frame()
mig = data.frame()

time_points = seq(0,3.5, 0.05)
change_points = seq(0,max_val-1,1)/(max_val-1)


for (i in seq(1, length(locations))){
  params = param.labels[startsWith(param.labels,paste("SkylineNe.", locations[[i]], ".",sep=""))]
  tmp.ne = c()
  for (j in seq(1, length(t.skygrid$Sample))){
    times = change_points*t.skygrid[j, "Tree.height"]
    vals = t.skygrid[j, params]
    tmp.ne = c(tmp.ne, approx(times, vals, xout=time_points, method = "linear")$y)
  }
  for (j in seq(1, length(time_points))){
    vals=exp(as.mcmc(tmp.ne[seq(j, length(tmp.ne), length(time_points))]))
    hpd.5 = HPDinterval(vals, prob=0.5)
    hpd.95 = HPDinterval(vals, prob=0.95)
    timestart = as.Date("2016-10-12")-time_points[[j]]*365
    dat = rbind(dat, data.frame(time = timestart,
                                            l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                            l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                            location=gsub("_"," ",locations[[i]]),
                                            method="skygrid"))
  }
}



for (i in seq(1,length(t.skygrid))){
  if (grepl("migrationEvents.", labels(t.skygrid)[[2]][[i]])){
    print(labels(t.skygrid)[[2]][[i]])
    loc = strsplit(gsub("migrationEvents.", "", labels(t.skygrid)[[2]][[i]]), split="_to_")[[1]]
    mig = rbind(mig, data.frame(events=median(t.skygrid[, i]),
                                from=gsub("_", " ", loc[[1]]),
                                to=gsub("_", " ", loc[[2]]),
                                method="skygrid"))
    
  }

}

dat.const = data.frame()
mig.constant = data.frame()

for (i in seq(1,length(t.constant))){
  if (grepl("ConstantNe.", labels(t.constant)[[2]][[i]])){
    loc = strsplit(gsub("ConstantNe.", "", labels(t.constant)[[2]][[i]]), split="\\.")[[1]]
    hpd.5 = HPDinterval(as.mcmc(t.constant[, i]), prob=0.5)
    hpd.95 = HPDinterval(as.mcmc(t.constant[, i]), prob=0.95)
    timestart = as.Date("2013-10-12")
    timeend= as.Date("2013-10-12")-0.025*(30)*365

    index = which(gsub("_", " ", loc[[1]])==labels(cols))
    dat.const = rbind(dat.const, data.frame(time = index-0.3,
                                l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                location=gsub("_", " ", loc[[1]]),
                                method="constant"))
    dat.const = rbind(dat.const, data.frame(time = index+0.3,
                                l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                location=gsub("_", " ", loc[[1]]),
                                method="constant"))

  }else if (grepl("migrationEvents.", labels(t.constant)[[2]][[i]])){
    # print(labels(t.skygrid)[[2]][[i]])
    loc = strsplit(gsub("migrationEvents.", "", labels(t.constant)[[2]][[i]]), split="_to_")[[1]]
    mig.constant = rbind(mig.constant, data.frame(events=median(t.constant[, i]),
                                from=gsub("_", " ", loc[[1]]),
                                to=gsub("_", " ", loc[[2]]),
                                method="skygrid"))

  }

}


p_ne1 <- ggplot(dat[which(dat$location=="Polynesia" | dat$location=="SouthAmerica" | dat$location=="CentralAmerica" | dat$location=="Caribbean"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location, group=interaction(method,location)), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location, group=interaction(method,location)), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  xlab("")+
  theme(legend.position = "top", 
        plot.margin=unit(c(0,0,0,0), "mm"))+
  coord_cartesian(ylim=c(0,100), xlim=c(as.Date("2013-01-01"),as.Date("2017-01-01")))
plot(p_ne1)

p_ne2 <- ggplot(dat[which(dat$location=="Brazil North" | dat$location=="Brazil Northeast" | dat$location=="Brazil Southeast"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  xlab("")+
  coord_cartesian(ylim=c(0,100), xlim=c(as.Date("2013-01-01"),as.Date("2017-01-01")))+
  theme(legend.position = "top", 
        plot.margin=unit(c(0,0,0,0), "mm"))
  plot(p_ne2)


p_ne3 <- ggplot(dat.const) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_continuous(breaks=seq(1,7),labels=labels(cols))+
  xlab("")+
  coord_cartesian(ylim=c(0,20)) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
plot(p_ne3)


# ggsave(plot=p_ne, filename = "../../../MascotSkyline-Text/Figures/zikv_ne.pdf", height=15, width=15)



system(" /Applications/BEAST\\ 2.7.4/bin/logcombiner -b 10 -log out/zikv_skyline_rep0.basta_alignment.trees out/zikv_skyline_rep1.basta_alignment.trees out/zikv_skyline_rep2.basta_alignment.trees -o out/zikv_skygrid.trees")
system(" /Applications/BEAST\\ 2.7.4/bin/logcombiner -b 10 -log out/zikv_constant_rep0.basta_alignment.trees out/zikv_constant_rep1.basta_alignment.trees out/zikv_constant_rep2.basta_alignment.trees -o out/zikv_constant.trees")


system(" /Applications/BEAST\\ 2.7.4/bin/treeannotator -b 0 out/zikv_skygrid.trees out/zikv_skygrid.tree")
system(" /Applications/BEAST\\ 2.7.4/bin/treeannotator -b 0 out/zikv_constant.trees out/zikv_constant.tree")

tree = read.beast("./out/zikv_skygrid.tree")

tree@data$location = c()
tree@data$entropy = c()

for (i in seq(1,length(tree@data$Brazil_North))){
  prob_vals = c()
  entr = 0
  for (j in seq(1,length(locations))){
    prob_vals[[j]] = as.numeric(tree@data[i, locations[[j]]])
    if (prob_vals[[j]]>0){
      entr = entr + prob_vals[[j]]*log(prob_vals[[j]])
    }
  }
  tree@data$location[i] = gsub("_", " ", locations[[which.max(prob_vals)]])
  tree@data$entropy[i] = entr
}

p_skygrid = ggtree(tree, aes(color=location), mrsd="2016-10-12") + 
  geom_nodepoint(aes(color=location),size=1)+
  geom_tippoint(aes(color=location),size=1)+
  theme_tree2() +
  scale_color_manual(values=cols)+
  # scale_alpha(range=c(0.2,1))+
  theme_minimal()+
  ggtitle("MASCOT-Skyline") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  coord_cartesian(xlim=c(2013,2017)) +
  coord_flip()  + 
  scale_x_reverse()

plot(p_skygrid)


tree = read.beast("./out/zikv_constant.tree")

tree@data$location = character(length(tree@data$Brazil_North))
tree@data$entropy = numeric(length(tree@data$Brazil_North))

for (i in seq(1,length(tree@data$Brazil_North))){
  entr = 0
  prob_vals = numeric(length(locations))
  for (j in seq(1, length(locations))) {
    prob_vals[j] = as.numeric(tree@data[i, locations[j]])
    entr = entr + prob_vals[j] * log(prob_vals[j])
  }
  tree@data$location[i] = gsub("_", " ", locations[which.max(prob_vals)])
  tree@data$entropy[i] = -entr
}

tree@data$location <- as.factor(tree@data$location)


library(dplyr)
require(maps)
require(viridis)
theme_set(
  theme_void()
)
p_constant = ggtree(tree, aes(color=location), mrsd="2016-10-12") +
  geom_nodepoint(aes(color=location),size=1)+
  geom_tippoint(aes(color=location),size=1)+
  theme_tree2() +
  scale_color_manual(values=cols)+
  theme_minimal()+
  ggtitle("MASCOT-Constant") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  coord_cartesian(xlim=c(2013,2017)) +
  coord_flip()  + scale_x_reverse()
plot(p_constant)




worldmap <- map_data("world")

# code from https://web.stanford.edu/~cengel/cgi-bin/anthrospace/great-circles-on-a-recentered-worldmap-in-ggplot
center = 260
worldmap$long.recenter <-  ifelse(worldmap$long  < center - 180 , worldmap$long + 360, worldmap$long)

RegroupElements <- function(df, longcol, idcol){  
  g <- rep(1, length(df[,longcol]))
  if (diff(range(df[,longcol])) > 300) {          # check if longitude within group differs more than 300 deg, ie if element was split
    d <- df[,longcol] > mean(range(df[,longcol])) # we use the mean to help us separate the extreme values
    g[!d] <- 1     # some marker for parts that stay in place (we cheat here a little, as we do not take into account concave polygons)
    g[d] <- 2      # parts that are moved
  }
  g <-  paste(df[, idcol], g, sep=".") # attach to id to create unique group variable for the dataset
  df$group.regroup <- g
  df
}

### Function to close regrouped polygons
# takes dataframe, checks if 1st and last longitude value are the same, if not, inserts first as last and reassigns order variable
ClosePolygons <- function(df, longcol, ordercol){
  if (df[1,longcol] != df[nrow(df),longcol]) {
    tmp <- df[1,]
    df <- rbind(df,tmp)
  }
  o <- c(1: nrow(df))  # rassign the order variable
  df[,ordercol] <- o
  df
}


worldmap.rg <- ddply(worldmap, .(group), RegroupElements, "long.recenter", "group")

# close polys
worldmap.cp <- ddply(worldmap.rg, .(group.regroup), ClosePolygons, "long.recenter", "order")  # use the new grouping var

locations = c("Polynesia", "SouthAmerica", "CentralAmerica", "Caribbean", "Brazil_North", "Brazil_Northeast", "Brazil_Southeast")
coord = data.frame(x=c(220,293,264,290,310,322,315),y=c(-17,7,18,20,0,-5,-19),locations=labels(cols))

for (i in seq(1,length(mig$events))){
  mig[i,c("xfrom","yfrom")] = coord[which(coord$locations==mig[i,"from"]),c("x","y")]
  mig[i,c("xto","yto")] = coord[which(coord$locations==mig[i,"to"]),c("x","y")]
  mig.constant[i,c("xfrom","yfrom")] = coord[which(coord$locations==mig.constant[i,"from"]),c("x","y")]
  mig.constant[i,c("xto","yto")] = coord[which(coord$locations==mig.constant[i,"to"]),c("x","y")]
  
}

mig = mig[-which(mig$events==0),]
mig.constant = mig.constant[-which(mig.constant$events==0),]

p.events.skl = ggplot() +
  geom_polygon(aes(long.recenter,lat,group=group.regroup), size = 0.5, fill="#f9f9f9", colour = "grey65", data=worldmap.cp) +
  coord_cartesian(xlim=c(210,320), ylim=c(-40,30)) +
  geom_point(data=coord, aes(x=x,y=y,color=locations), size=2) +
  geom_curve(data=mig, aes(x=xfrom,y=yfrom,xend=xto,yend=yto,color=from,size=events),
             arrow = arrow(length = unit(0.03, "npc"),type = "closed"), curvature=0.3)+
  scale_color_manual(values=cols)+
  # theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin=unit(c(8,8,8,8), "mm"))+
  ggtitle("MASCOT-Skyline") +
  scale_size_continuous(name="number of\nevents",limits=c(1,60),breaks=c(1,20,40),range = c(0.25, 2)) +
  guides(color=FALSE) +
  theme(legend.position = c(0.2, -0.25))
plot(p.events.skl)


p.events.con = ggplot() +
  geom_polygon(aes(long.recenter,lat,group=group.regroup), size = 0.5, fill="#f9f9f9", colour = "grey65", data=worldmap.cp) +
  coord_cartesian(xlim=c(210,320), ylim=c(-40,30)) +
  geom_point(data=coord, aes(x=x,y=y,color=locations),size=2) +
  geom_curve(data=mig.constant, aes(x=xfrom,y=yfrom,xend=xto,yend=yto,color=from,size=events),
             arrow = arrow(length = unit(0.03, "npc"),type = "closed"), curvature=0.3)+
  scale_color_manual(values=cols) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin=unit(c(8,8,8,8), "mm"))+
  ggtitle("MASCOT-Constant") +
  scale_size_continuous(limits=c(1,60),range = c(0.25, 2))


plot(p.events.con)

  
library("gridExtra")



p1 = ggarrange(p_skygrid, p_constant, ncol = 2, labels = c("A", "B"))
p2 = ggarrange(p_ne1, p_ne2, ncol = 1, labels = c("C", "E"))
p3 = ggarrange(p.events.skl, p.events.con, ncol = 1, labels = c("D", "F"))

p4 = ggarrange(p2, p3, ncol = 2, widths =c(2,1))
p_complete = ggarrange(p1, p4, ncol = 1, heights=c(0.8,1))
plot(p_complete)

ggsave(plot=p1, filename = "../../../MascotSkyline-Text/Figures/zikv_trees.pdf", height=4, width=10)
ggsave(plot=p2, filename = "../../../MascotSkyline-Text/Figures/zikv_ne.pdf", height=4, width=7)
ggsave(plot=p_complete, filename = "../../../MascotSkyline-Text/Figures/zikv.pdf", height=9, width=12)


