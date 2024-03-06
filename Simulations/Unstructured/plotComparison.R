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

rootHeight=3.6771
timepoints = seq(0,1,0.01)

first = TRUE

names = c("constant", "exponential")
p = list()

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  # Read in the logs
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$m1)/10)), ]

  # calculate ess values
  ess <- effectiveSize(t)

  # compute trajectories
  for (s in seq(1,2)){
    for (j in seq(1,length(timepoints))){
      Ne = c()
      for (k in seq(1,length(t$Sample))){
        Ne = append(Ne, t[k,paste("Nenull",s,sep=".")]-t[k,paste("growth",s,sep=".")]*timepoints[j]*rootHeight)
      }
      Ne = exp(Ne)
      hpd = HPDinterval(as.mcmc(Ne))
      new.dat = data.frame(time=timepoints[j]*rootHeight, Ne=mean(Ne), lower=hpd[1,"lower"], upper=hpd[1,"upper"], state=paste("state",s-1, sep=""), run=names[[i]])
      if (first){
        first=F
        data = new.dat
      }else{
        data = rbind(data, new.dat)
      }
    }
  }
  system(paste("/Applications/BEAST\\ 2.6.4/bin/treeannotator -burnin 10", gsub(".log", ".trees", log[i]), "./out/summary.tree" ))
  x <- read.beast("./out/summary.tree")
  p[[i]] = ggtree(x, right=TRUE, aes(color=max)) +
    scale_color_OkabeIto() +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(breaks=NULL) +
    ylab("") +
    xlab("time")
  
  
  new.arrow = data.frame(x=1.1,xend=1.9,y=1.1,yend=1.1,size=mean(t$f_migration.state1_to_state0), state="state1", run=names[[i]])
  new.arrow = rbind(new.arrow, data.frame(x=1.9,xend=1.1,y=0.9,yend=0.9,size=mean(t$f_migration.state0_to_state1), state="state0", run=names[[i]]))
  if(i==1){
    arrow = new.arrow
  }else{
    arrow = rbind(arrow, new.arrow)
  }
}


# make the Ne plots
red.data1 = data[which(data$run=="constant"), ]
red.data2 = data[which(data$run=="exponential" & data$state=="state0"), ]

p_ne1 <- ggplot(red.data1) +
  geom_line(aes(x=abs(time-rootHeight), y=Ne, color=state, group=state)) +
  geom_ribbon(aes(x=abs(time-rootHeight), ymin=lower,ymax=upper, color=state, group=state, fill=state), alpha=0.1) +
  theme_minimal()+
  scale_color_OkabeIto(breaks=c("state0","state1"))+
  scale_fill_OkabeIto()+
  scale_y_log10(limits=c(0.5,120)) +
  xlab("")+
  theme(legend.position = "none")

p_ne2 <- ggplot(red.data2) +
  geom_line(aes(x=abs(time-rootHeight), y=Ne, color=state, group=state)) +
  geom_ribbon(aes(x=abs(time-rootHeight), ymin=lower,ymax=upper, color=state, group=state, fill=state), alpha=0.1) +
  theme_minimal()+
  scale_color_OkabeIto(breaks=c("state0","state1"))+
  scale_fill_OkabeIto()+
  scale_y_log10(limits=c(0.5,120)) +
  xlab("")+
  theme(legend.position = "none")



red.arrow1 = arrow[which(arrow$run=="constant"), ]
red.arrow2 = arrow[which(arrow$run=="exponential"), ]


p_arrow1 <- ggplot(red.arrow1) +
  geom_point(aes(x=round(x),y=round(y), color=state), size=10) +
  geom_curve(aes(x=x,y=y,xend=xend,yend=yend, color=state, size=size),curvature=-0.5,arrow = arrow(length = unit(0.03, "npc"))) +
  scale_color_OkabeIto(order=c(2,1))+
  scale_y_continuous(limits=c(0.5,1.5))+
  scale_x_continuous(limits=c(0.5,2.5))+
  scale_size_continuous(limits=c(0,max(arrow$size)), range=c(0,3)) +
  theme_void()+
  theme(legend.position = "none")

p_arrow2 <- ggplot(red.arrow2) +
  geom_point(aes(x=round(x),y=round(y), color=state), size=10) +
  geom_curve(aes(x=x,y=y,xend=xend,yend=yend, color=state, size=size),curvature=-0.5,arrow = arrow(length = unit(0.03, "npc"))) +
  scale_color_OkabeIto(order=c(2,1))+
  scale_y_continuous(limits=c(0.5,1.5))+
  scale_x_continuous(limits=c(0.5,2.5))+
  scale_size_continuous(limits=c(0,max(arrow$size)), range=c(0,3)) +
  theme_void()+
  theme(legend.position = "none")

theme_set(theme_pubr())



figure <- ggarrange(p[[1]], p_ne1, p_arrow1, p[[2]], p_ne2, p_arrow2,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 3, nrow = 2)

plot(figure)

ggsave(plot=figure, paste("../../../MascotSkyline-Text/Figures/expo_unstructured.pdf", sep=""),width=10, height=3.5)


