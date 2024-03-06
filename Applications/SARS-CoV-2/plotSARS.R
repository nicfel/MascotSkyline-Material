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

locations = c("WesternWashington", "EasternWashington",
              "NorthAmerica", "restOfWorld")

cols = c("WesternWashington"="#A4A9AA", 
"EasternWashington"="#03A381", "NorthAmerica"="#B39B31", "restOfWorld"="#F3756D")


# read in all *.log files in the out folder that have *mascot* in the pattern 
log <- list.files(  # nolint
  path = "./out",
  pattern = "*mascot*",
  full.names = TRUE
)
# remove all files that don't end in log
log <- log[endsWith(log, "log")]
log <- log[!grepl("it", log)]


# Read in the logs
t1 <- read.table(log[[1]], header=TRUE, sep="\t")
t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
t2 <- read.table(log[[2]], header=TRUE, sep="\t")
t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
t3 <- read.table(log[[3]], header=TRUE, sep="\t")
t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]

t.skyline = rbind(t1,t2,t3)

# read in all *.log files in the out folder that have *dta* in the pattern
log <- list.files(  # nolint
  path = "./out",
  pattern = "*dta*",
  full.names = TRUE
)
# remove all files that don't end in log
log <- log[endsWith(log, "log")]
log <- log[!grepl("it", log)]

t1 <- read.table(log[[1]], header=TRUE, sep="\t")
t1 <- t1[-seq(1,ceiling(length(t1$Sample)/10)), ]
t2 <- read.table(log[[2]], header=TRUE, sep="\t")
t2 <- t2[-seq(1,ceiling(length(t2$Sample)/10)), ]
t3 <- read.table(log[[3]], header=TRUE, sep="\t")
t3 <- t3[-seq(1,ceiling(length(t3$Sample)/10)), ]

t.dta = rbind(t1,t2,t3)



# get all the labels of t.skyline that start with SkylineNe.WesternWashington
param.labels = labels(t.skyline)[[2]]
wa.labels = param.labels[startsWith(param.labels,"SkylineNe.WesternWashington.")]

time_points = seq(0,1.1,1.1/(length(wa.labels)-1))  


dat = data.frame()
mig = data.frame()

change_points = seq(0,1.1, 1.1/(length(wa.labels)-1))
time_points = seq(0,1.1, 1.1/100)


for (i in seq(1, length(locations))){
  params = param.labels[startsWith(param.labels,paste("SkylineNe.", locations[[i]], ".",sep=""))]
  tmp.ne = c()
  for (j in seq(1, length(t.skyline$Sample))){
    times = change_points*t.skyline[j, "Tree.height"]
    vals = t.skyline[j, params]
    tmp.ne = c(tmp.ne, approx(times, vals, xout=time_points, method = "linear")$y)
  }
  
  for (j in seq(1, length(time_points))){
    vals=exp(as.mcmc(tmp.ne[seq(j, length(tmp.ne), length(time_points))]))
    hpd.5 = HPDinterval(vals, prob=0.5)
    hpd.95 = HPDinterval(vals, prob=0.95)
    timestart = as.Date("2020-10-14")-time_points[[j]]*365
    dat = rbind(dat, data.frame(time = timestart,
                                            l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                            l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                            location=gsub("_"," ",locations[[i]]),
                                            method="MASCOT-Skyline"))
  }
}





# dat.const = data.frame()
# mig.constant = data.frame()

# for (i in seq(1,length(t.constant))){
#   if (grepl("Ne.", labels(t.constant)[[2]][[i]])){
#     loc = strsplit(gsub("Ne.", "", labels(t.constant)[[2]][[i]]), split="\\.")[[1]]
#     hpd.5 = HPDinterval(as.mcmc(t.constant[, i]), prob=0.5)
#     hpd.95 = HPDinterval(as.mcmc(t.constant[, i]), prob=0.95)
#     timestart = as.Date("2013-10-12")
#     timeend= as.Date("2013-10-12")-0.025*(30)*365

#     index = which(gsub("_", " ", loc[[1]])==labels(cols))
#     dat.const = rbind(dat.const, data.frame(time = index-0.3,
#                                 l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
#                                 l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
#                                 location=gsub("_", " ", loc[[1]]),
#                                 method="constant"))
#     dat.const = rbind(dat.const, data.frame(time = index+0.3,
#                                 l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
#                                 l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
#                                 location=gsub("_", " ", loc[[1]]),
#                                 method="constant"))

#   }else if (grepl("migrationEvents.", labels(t.constant)[[2]][[i]])){
#     # print(labels(t.skygrid)[[2]][[i]])
#     loc = strsplit(gsub("migrationEvents.", "", labels(t.constant)[[2]][[i]]), split="_to_")[[1]]
#     mig.constant = rbind(mig.constant, data.frame(events=median(t.constant[, i]),
#                                 from=gsub("_", " ", loc[[1]]),
#                                 to=gsub("_", " ", loc[[2]]),
#                                 method="skygrid"))

#   }

# }


p_ne1 <- ggplot(dat[which(dat$location=="WesternWashington" | dat$location=="EasternWashington"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location, group=interaction(method,location)), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location, group=interaction(method,location)), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  xlab("")+
  theme(legend.position = "top")+
  coord_cartesian(ylim=c(0,5), xlim=c(as.Date("2020-02-01"),as.Date("2020-10-14")))
plot(p_ne1)

p_ne1 <- ggplot(dat[which(dat$location=="NorthAmerica" | dat$location=="EasternWashington"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location, group=interaction(method,location)), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location, group=interaction(method,location)), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  xlab("")+
  theme(legend.position = "top")+
  coord_cartesian(ylim=c(0,10), xlim=c(as.Date("2020-02-01"),as.Date("2020-10-14")))
plot(p_ne1)

p_ne1 <- ggplot(dat[which(dat$location=="NorthAmerica" | dat$location=="restOfWorld"),]) +
  geom_ribbon(aes(x=time, ymin=l.95,ymax=u.95, fill=location, group=interaction(method,location)), alpha=0.25) +
  geom_ribbon(aes(x=time, ymin=l.5,ymax=u.5, fill=location, group=interaction(method,location)), alpha=0.9) +
  theme_minimal()+
  scale_fill_manual(name="",values=cols)+
  scale_x_date()+
  xlab("")+
  theme(legend.position = "top")+
  coord_cartesian(ylim=c(0,50), xlim=c(as.Date("2020-02-01"),as.Date("2020-10-14")))
plot(p_ne1)









for (a in seq(1, length(locations))){
  for (b in seq(1, length(locations))){
    if (a!=b){
      name = paste("migrationEvents.", locations[[a]], "_to_", locations[[b]], sep="")
      name.rates = paste("f_migrationRatesSkyline.", locations[[a]], "_to_", locations[[b]], sep="")
      
      # put the vals and rates into vectors sorted 
      vals = t.skyline[,name]
      rates = t.skyline[,name.rates]
      # vals = vals[order(vals)]
      # rates = rates[order(rates)]

      # remove the first and last 2.5% of entries from vals and rates
      # vals = vals[seq(ceiling(length(vals)*0.025), length(vals)-ceiling(length(vals)*0.025))]
      # rates = rates[seq(ceiling(length(rates)*0.025), length(rates)-ceiling(length(rates)*0.025))]
      
      
      mig = rbind(mig, data.frame(events=vals,
                                  rates = rates,
                                  from=locations[[a]],
                                  to=locations[[b]],
                                  method="MASCOT-Skyline"))
    }
  }
}

for (a in seq(1, length(locations))){
  for (b in seq(1, length(locations))){
    if (a!=b){
      name = paste("c_migrationEvents.", locations[[a]], "_to_", locations[[b]], ".1.", sep="")
      name.rates = paste("geo.rates.", locations[[a]], ".", locations[[b]], sep="")
      
      # put the vals and rates into vectors sorted 
      vals = t.dta[,name]
      rates = t.dta[,name.rates]*t.dta[, "geo.clock.rate"]
      # vals = vals[order(vals)]
      # rates = rates[order(rates)]

      # remove the first and last 2.5% of entries from vals and rates
      # vals = vals[seq(ceiling(length(vals)*0.025), length(vals)-ceiling(length(vals)*0.025))]
      # rates = rates[seq(ceiling(length(rates)*0.025), length(rates)-ceiling(length(rates)*0.025))]
      
      mig = rbind(mig, data.frame(events=vals,
                                  rates=rates,
                                  from=locations[[a]],
                                  to=locations[[b]],
                                  method="DTA"))
    }
  }
}

# change the order of from and to facets to the same as locations
mig$from = factor(mig$from, levels=locations)
mig$to = factor(mig$to, levels=locations)

p_ne3 <- ggplot(mig) +
  # plot as a violin plot, but remove the smoothing, remove outliers and adjust the width of the violins
  geom_violin(aes(x=method, y=events, fill=method, color=method), trim=T, adjust=0.01, size=1)+
  facet_grid(from~to, scales="free", space="free", switch = "y")+
  xlab("")+
  # scale_y_log10()+
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  # remove vertical lines in plots
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  # move the facet headers outside of the y labels
  theme(strip.placement = "outside", 
        strip.background = element_blank(),
        legend.position = "none")+
  scale_y_continuous(position = "right")
plot(p_ne3)


p_ne3 <- ggplot(mig) +
  geom_violin(aes(x=method, y=rates, fill=method)) +
  facet_grid(from~to, scales="free", space="free")+
  xlab("")+
  scale_y_log10()+
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  # remove vertical lines in plots
  theme_minimal()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
plot(p_ne3)

# vector to shorten names
locations.short = c("WW"="WesternWashington", "EW"="EasternWashington",
              "NA"="NorthAmerica", "RoW"="restOfWorld")



# start a second data frame for mig that contains as the rates the ratios between for and to in both directions,
# i.e. from A to B and from B to A
mig.ratios = data.frame()
for (a in seq(1, length(locations)-1)){
  for (b in seq(a+1, length(locations))){
    # also loop over both methods
    for (c in c("MASCOT-Skyline", "DTA")){
      # get the rates from A to B and from B to A
      rates_from = mig[which(mig$from==locations[[a]] & mig$to==locations[[b]] & mig$method==c), "rates"]
      rates_to = mig[which(mig$from==locations[[b]] & mig$to==locations[[a]] & mig$method==c), "rates"]
      events_from = mig[which(mig$from==locations[[a]] & mig$to==locations[[b]] & mig$method==c), "events"]
      events_to = mig[which(mig$from==locations[[b]] & mig$to==locations[[a]] & mig$method==c), "events"]

      from.short = names(locations.short)[locations.short==locations[[a]]]
      to.short = names(locations.short)[locations.short==locations[[b]]]
      
      # calculate the ratio of the rates
      mig.ratios = rbind(mig.ratios, data.frame(rates=rates_from/rates_to,
                                                events=(events_from+0.1) /(events_to+0.1),
                                                direction = paste(from.short, " & ", to.short, sep=""),
                                                from=locations[[a]],
                                                to=locations[[b]],
                                                method=c))
      


    }
  }
}


# Convert interaction to a numeric variable
mig.ratios$x_numeric <- as.numeric(interaction(mig.ratios$method, mig.ratios$direction))

# Find unique labels and calculate midpoints
unique_labels <- unique(gsub(".*\\.", "", interaction(mig.ratios$method, mig.ratios$direction)))
midpoints <- sapply(unique(gsub(".*\\.", "", mig.ratios$direction)), function(label) {
  mean(mig.ratios$x_numeric[gsub(".*\\.", "", interaction(mig.ratios$method, mig.ratios$direction)) == label])
})

# Plot
p_mig_ratio <- ggplot(mig.ratios, aes(x = x_numeric, y = rates, fill = method, color = method, group=x_numeric)) +
  geom_violin() +
  # stat_summary(
  #   aes(label = sprintf("%.2f", ..y..)),
  #   fun = median,
  #   geom = "text",
  #   vjust = 0,
  #   color = "black"
  # ) +
  theme_minimal() +
  scale_y_log10() +
  scale_x_continuous(breaks = midpoints, labels = unique(gsub(".*\\.", "", unique_labels))) +
  xlab("")+
  ylab("rate ratio") +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))

# Display the plot
plot(p_mig_ratio)


sampled_mig <- mig[sample(1:nrow(mig), 0.05 * nrow(mig)), ]

p_dir_corr_all = ggplot(mig, aes(x=events+0.1, y=rates, color=method, fill=method)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point(data=sampled_mig, alpha=0.1) +
  geom_smooth(method="lm")+
  stat_cor(aes(label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  facet_grid(from~to)+
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  theme_minimal()

ggsave(p_dir_corr_all, filename = "../../../MascotSkyline-Text/Figures/sars_dircorr_all.pdf", height=8, width=10)

mig.ratios2 = data.frame()
for (a in seq(1, length(locations)-1)){
  for (b in seq(1, length(locations))){
    if (a!=b){
      # also loop over both methods
      for (c in c("MASCOT-Skyline", "DTA")){
        # get the rates from A to B and from B to A
        rates = mig[which(mig$from==locations[[a]] & mig$to==locations[[b]] & mig$method==c), "rates"]
        events = mig[which(mig$from==locations[[a]] & mig$to==locations[[b]] & mig$method==c), "events"]+0.1
        # compute the 95 % HPD for rates and events
        hpd.95.rates = HPDinterval(as.mcmc(rates), prob=0.95)
        hpd.95.events = HPDinterval(as.mcmc(events), prob=0.95)        
        # calculate the ratio of the rates
        mig.ratios2 = rbind(mig.ratios2, data.frame(rates=mean(rates),
                                                  rates.lower=hpd.95.rates[1,"lower"],
                                                  rates.upper=hpd.95.rates[1,"upper"],
                                                  events.lower=hpd.95.events[1,"lower"],
                                                  events.upper=hpd.95.events[1,"upper"],
                                                  events=mean(events),
                                                  from=locations[[a]],
                                                  to=locations[[b]],
                                                  method=c))
      }
    }
  }
}

p_dir_corr = ggplot(mig.ratios2, aes(x=events, y=rates, color=method, fill=method)) +
  scale_y_log10() +
  scale_x_log10() +
  # geom_errorbar(aes(ymin=rates.lower, ymax=rates.upper), width=0.1)+
  # geom_errorbarh(aes(xmin=events.lower, xmax=events.upper), height=0.1)+
  geom_smooth(method="lm")+
  stat_cor(aes(label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  geom_point() +theme_minimal()+
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510")) +
  theme(legend.position = "top")


ggsave(p_dir_corr, filename = "../../../MascotSkyline-Text/Figures/sars_dircorr.pdf", height=3, width=4)

  # coord_cartesian(xlim=c(0.1, 50), ylim=c(0.001,10))

events.corr = data.frame()
for (from in unique(mig$from)){
  for (to in unique(mig$to)){
    if (from!=to){
      val.mascot=mig[mig$from==from & mig$to==to & mig$method=="MASCOT-Skyline", "events"]
      val.dta=mig[mig$from==from & mig$to==to & mig$method=="DTA", "events"]
      
      median.mascot = mean(val.mascot)
      median.dta = mean(val.dta)
      
      hpd.m = HPDinterval(as.mcmc(val.mascot))
      hpd.d = HPDinterval(as.mcmc(val.dta))
      
      events.corr = rbind(events.corr, data.frame(x=median.mascot, y=median.dta,
                                                  xmin=hpd.m[1,"lower"],xmax=hpd.m[1,"upper"],
                                                  ymin=hpd.d[1,"lower"],ymax=hpd.d[1,"upper"]))
    }
  }
}


p_events_corr = ggplot(events.corr, aes(x=x, y=y,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
  geom_errorbar(alpha=0.2)+
  geom_errorbarh(alpha=0.2)+
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method="lm")+
  stat_cor(aes(label = ..r.label..), label.x = -1, p.accuracy = 0.01, r.accuracy = 0.01)+
  geom_point() +theme_minimal() +
  xlab("mean number of events MASCOT-Skyline")+
  ylab("mean number of events DTA")
plot(p_events_corr)
ggsave(p_events_corr, filename = "../../../MascotSkyline-Text/Figures/sars_events_corr.pdf", height=4, width=6)


p_migevents_ratio <- ggplot(mig.ratios, aes(x = x_numeric, y = events, fill = method, color = method, group=x_numeric)) +
  geom_violin() +
  # stat_summary(
  #   aes(label = sprintf("%.2f", ..y..)),
  #   fun = median,
  #   geom = "text",
  #   vjust = 0,
  #   color = "black"
  # ) +
  theme_minimal() +
  scale_y_log10() +
  scale_x_continuous(breaks = midpoints, labels = unique(gsub(".*\\.", "", unique_labels))) +
  xlab("")+
  ylab("events ratio") +
  scale_color_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))+
  scale_fill_manual(values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510"))


# Display the plot
plot(p_migevents_ratio)


library(gridExtra)
p.ratios = grid.arrange(p_mig_ratio, p_migevents_ratio, ncol = 1)
plot(p.ratios)
ggsave(p.ratios, filename = "../../../MascotSkyline-Text/Figures/sars_ratios.pdf", height=5, width=7)







# system(" /Applications/BEAST\\ 2.7.4/bin/logcombiner -resample 500000 -b 10 -log out/SARS2_mascot_rep0.SARSCoV2_WA.trees out/SARS2_mascot_rep1.SARSCoV2_WA.trees out/SARS2_mascot_rep2.SARSCoV2_WA.trees -o out/SARS2_mascot.trees")
# system(" /Applications/BEAST\\ 2.7.4/bin/logcombiner -resample 500000 -b 10 -log out/SARS2_dta_rep0.trees out/SARS2_dta_rep1.trees out/SARS2_dta_rep2.trees -o out/SARS2_dta.trees")
# 
# 
# system(" /Applications/BEAST\\ 2.7.4/bin/treeannotator -b 0 -heights mean out/SARS2_mascot.trees out/SARS2_mascot.tree")
# system(" /Applications/BEAST\\ 2.7.4/bin/treeannotator -b 0 -heights mean out/SARS2_dta.trees out/SARS2_dta.tree")

tree = read.beast("./out/SARS2_mascot.tree")

tree@data$location = c()
tree@data$entropy = c()

for (i in seq(1,length(tree@data$WesternWashington))){
  prob_vals = c()
  entr = 0
  for (j in seq(1,length(locations))){
    prob_vals[[j]] = as.numeric(tree@data[i, locations[[j]]])
    if (prob_vals[[j]]>0){
      entr = entr + prob_vals[[j]]*log(prob_vals[[j]])
    }
  }
  tree@data$location[i] = locations[[which.max(prob_vals)]]
  tree@data$entropy[i] = entr
}

p_mascot_tree = ggtree(tree, aes(color=location), mrsd="2020-10-13") + 
  # geom_nodepoint(aes(color=location),size=1)+
  geom_tippoint(aes(color=location),size=2)+
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
  # coord_cartesian(xlim=c(2013,2017)) +
  coord_flip()  + 
  scale_x_reverse()

plot(p_mascot_tree)


tree.dta = read.beast("./out/SARS2_dta.tree")


library(dplyr)
require(maps)
require(viridis)
theme_set(
  theme_void()
)

p_dta = ggtree(tree.dta, aes(color=geo), mrsd="2020-10-13") +
  # geom_nodepoint(aes(color=geo),size=1)+
  geom_tippoint(aes(color=geo),size=2)+
  theme_tree2() +
  scale_color_manual(values=cols)+
  theme_minimal()+
  ggtitle("DTA") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  # coord_cartesian(xlim=c(2013,2017)) +
  coord_flip()  + scale_x_reverse()

plot(p_dta)

p = grid.arrange(p_mascot_tree, p_dta, ncol = 1)

ggsave(plot=p, filename = "../../../MascotSkyline-Text/Figures/sars_trees.pdf", height=6, width=7.5)



p_mascot_tree_legend_customized <- p_mascot_tree + 
  theme(legend.position = "bottom") +
  scale_color_manual(values=cols, labels=c("Eastern Washington (EW)", "North America (NA)", "Rest of World (ROW)", "Western Washington (WW)"))+
  guides(
    fill = guide_legend(ncol = 2),
    color = guide_legend(ncol = 2)
  )


plot(p_mascot_tree_legend_customized)
ggsave(plot=p_mascot_tree_legend_customized, filename = "../../../MascotSkyline-Text/Figures/sars_trees_legend.pdf", height=2, width=5)




