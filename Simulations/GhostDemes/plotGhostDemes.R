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
library(ggpubr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the parameters
params = read.table("rates.txt", header=TRUE, sep="\t", nrows = 0)
param_labels = labels(params)[[2]]

# get the names of all SISCO first (of three) run log files
log <- list.files(path="..//out", pattern="GhostDemes.*.log", full.names = TRUE)

data = data.frame()

# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")

  # Read in the SISCO *.logs
  t <- read.table(filename1, header=TRUE, sep="\t")
  t <- t[-seq(1,ceiling(length(t$Sample)/10)), ]
  
  if (length(t$Sample)>10){
    # calculate ess values
    ess <- effectiveSize(t)
    if (min(ess[2:3])<100){
      print("masco ESS value to low")
      print(sprintf("ESS value is %f for file %s",min(ess[2:6]),filename1))
    }else{
      dfname <- data.frame(filename = filename1)
      
      # get the number of the run
      tmp = strsplit(filename1, split="_")
      tmp = strsplit(tmp[[1]][2], split="S")
      
      runnumber = as.numeric(tmp[[1]][2])
      param_index = which(params$run==runnumber)
      # loop over the Ne0's
      for (j in seq(2,length(param_labels))){
        if (grepl("mascot", log[[i]])){
          mean = mean(t[,param_labels[j]])
          median = mean(t[,param_labels[j]])
          hpd = HPDinterval(as.mcmc(t[,param_labels[j]]))
          value_name = strsplit(param_labels[j], split="\\.")[[1]][[1]]
          true=params[param_index, param_labels[j]]
          
          if (grepl("Ne", value_name)){
            data = rbind(data, data.frame(name = value_name, true=true, mean=mean, median=median, 
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                          time = as.numeric(strsplit(param_labels[j], split="\\.")[[1]][[3]]),
                                          state = strsplit(param_labels[j], split="\\.")[[1]][[2]],
                                          run=runnumber))
          }else{
            data = rbind(data, data.frame(name = value_name, true=true, mean=mean, median=median, 
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"], 
                                          time = NA,
                                          state=strsplit(strsplit(param_labels[j], split="\\.")[[1]][[2]], split="_to_")[[1]][[1]],
                                          run=runnumber))
          }
        }else{
          if (grepl("SkylineNe.state0", param_labels[j])){
            mean = mean(t[,param_labels[j]])
            median = mean(t[,param_labels[j]])
            hpd = HPDinterval(as.mcmc(t[,param_labels[j]]))
            true=params[param_index, param_labels[j]]
            value_name = strsplit(param_labels[j], split="\\.")[[1]][[1]]
            data = rbind(data, data.frame(name = value_name, true=NA, mean=mean, median=median,
                                          lower=hpd[1,"lower"], upper=hpd[1,"upper"],
                                          time = as.numeric(strsplit(param_labels[j], split="\\.")[[1]][[3]]),
                                          state = "combined",
                                          run=runnumber))
          }
        }
      }
      
      # read in the single child trees
      # das
      # treefname = gsub("1mascot.log",".single.tree", gsub("../out/","./master/",log[[i]]))
      # trees = read.beast(truetrees[[tr]]) # read in the tree file
      # edges = trees@phylo$edge # get the edges of the tree
      # nodes = as.numeric(trees@data$node) # get the nodes of the tree
      
      
    }
  }
}


subset = data[data$name=="SkylineNe" & data$state=="state0",]
tot.ne.sampled = subset[,"lower"]<subset[,"true"] &
  subset[,"upper"]>subset[,"true"]
cov.ne.sampled = sum(tot.ne.sampled[!is.na(tot.ne.sampled)])/sum(!is.na(tot.ne.sampled))
subset = data[data$name=="SkylineNe" & data$state=="state1",]
tot.ne = subset[,"lower"]<subset[,"true"] &
  subset[,"upper"]>subset[,"true"]
cov.ne = sum(tot.ne[!is.na(tot.ne)])/sum(!is.na(tot.ne))

subset = data[data$name=="f_migrationRatesSkyline" & data$state=="state1",]
tot.mig_sampled = subset[,"lower"]<subset[,"true"] &
  subset[,"upper"]>subset[,"true"]
cov.mig.sampled = sum(tot.mig_sampled[!is.na(tot.mig_sampled)])/sum(!is.na(tot.mig_sampled))
subset = data[data$name=="f_migrationRatesSkyline" & data$state=="state0",]
tot.mig = subset[,"lower"]<subset[,"true"] &
  subset[,"upper"]>subset[,"true"]
cov.mig = sum(tot.mig[!is.na(tot.mig)])/sum(!is.na(tot.mig))

p_ne_sampled <- ggplot(data[data$name == "SkylineNe" & data$state=="state0",]) +
  geom_point(aes(x = true, y = median)) +
  geom_abline(color = "red") +
  geom_errorbar(aes(x = true, ymin = lower, ymax = upper), alpha = 0.5) +
  ylab("estimated") + xlab("true") +
  ggtitle(paste("log(Ne) in sampled\ndeme ( cov =", round(cov.ne.sampled, 2), ")")) +
  theme_minimal() +
  theme(legend.position = "none")
plot(p_ne_sampled)

p_ne_ghost <- ggplot(data[data$name == "SkylineNe" & data$state=="state1",]) +
  geom_point(aes(x = true, y = median)) +
  geom_abline(color = "red") +
  geom_errorbar(aes(x = true, ymin = lower, ymax = upper), alpha = 0.5) +
  ylab("estimated") + xlab("true") +
  ggtitle(paste("log(Ne) in ghost\ndeme ( cov =", round(cov.ne, 2), ")")) +
  theme_minimal() +
  theme(legend.position = "none")
plot(p_ne_ghost)

p_migration_sampled <- ggplot(data[data$name == "f_migrationRatesSkyline" & data$state=="state1",]) +
  geom_point(aes(x = true, y = median)) +
  geom_errorbar(aes(x = true, ymin = lower, ymax = upper), alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(color = "red") +
  ylab("estimated") + xlab("true") +
  ggtitle(paste("Migration rates from\nghost deme ( cov =", round(cov.mig.sampled, 2), ")")) +
  theme_minimal() +
  theme(legend.position = "none")
plot(p_migration_sampled)

p_migration_ghost <- ggplot(data[data$name == "f_migrationRatesSkyline" & data$state=="state0",]) +
  geom_point(aes(x = true, y = median)) +
  geom_errorbar(aes(x = true, ymin = lower, ymax = upper), alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(color = "red") +
  ylab("estimated") + xlab("true") +
  ggtitle(paste("Migration rates into\nghost deme ( cov =", round(cov.mig, 2), ")")) +
  theme_minimal() +
  theme(legend.position = "none")
plot(p_migration_ghost)


combined_plot <- ggarrange(p_ne_sampled, p_migration_sampled,p_ne_ghost,p_migration_ghost, labels = c("A", "B", "C", "D"), ncol = 2, nrow=2)

# Display the combined plot
print(combined_plot)


ggsave(plot=combined_plot,paste("../../../MascotSkyline-Text/Figures/ghostdeme_params.pdf", sep=""),width=8, height=8)

set.seed(1223)

runs=unique(data[data$state!="combined", "run"])
selected_runs <- sample(runs, 4)

# plot the skyline graphs
p_skyline <- ggplot(data=data[data$name=="SkylineNe" & data$run %in% selected_runs & data$state!="combined",], aes(x=time, color=state, fill=state, group=state))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
  geom_line(aes(y=median,linetype = "median estimate"))+
  geom_line(aes(y=true, linetype = "true value"))+
  facet_wrap(run~., ncol=2) +
  ylab("log Ne") + 
  xlab("time") + 
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme_minimal()  +
  scale_color_OkabeIto(labels=c("sampled deme", "ghost deme"))+
  scale_fill_OkabeIto(labels=c("sampled deme", "ghost deme"))
labels <- setNames(paste("random replicate", 1:length(selected_runs)), sort(selected_runs))

# Use the labels in facet_wrap
p_skyline <- p_skyline + 
  facet_wrap(~ run, ncol=2, labeller = labeller(run = labels))

plot(p_skyline)

ggsave(plot=p_skyline,paste("../../../MascotSkyline-Text/Figures/ghostdeme_trends.pdf", sep=""),width=8, height=8)



keep_runs = c()
for (i in runs){
  if (sum(data$run==i)>30){
    keep_runs = c(keep_runs, i)
  }
}

selected_runs <- sample(keep_runs, 20)

data$state = factor(data$state, levels=c("state0", "state1", "combined"))

p_skyline <- ggplot(data=data[data$name=="SkylineNe" & data$run %in% selected_runs,], aes(x=time, color=state, fill=state, group=state))+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,color=NA)+
  geom_line(aes(y=median,linetype = "median estimate"))+
  geom_line(aes(y=true, linetype = "true value"))+
  facet_wrap(run~., ncol=5) +
  ylab("log Ne") + 
  xlab("time") + 
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank()
  ) +
  theme_minimal()  +
  scale_color_OkabeIto(labels=c("sampled deme", "ghost deme", "one state"))+
  scale_fill_OkabeIto(labels=c("sampled deme", "ghost deme", "one state"))
plot(p_skyline)
labels <- setNames(paste("random replicate", 1:length(selected_runs)), sort(selected_runs))
p_skyline <- p_skyline + 
  facet_wrap(~ run, ncol=5, labeller = labeller(run = labels))

ggsave(plot=p_skyline,paste("../../../MascotSkyline-Text/Figures/ghostdeme_trends_onestate.pdf", sep=""),width=12, height=9)

