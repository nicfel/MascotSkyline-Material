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

# define the states
states = c("EasternWashington", "NorthAmerica", "WesternWashington", "restOfWorld")

# get all files in the out directory that contain "mascot" and "it" and end in log
log <- list.files(  # nolint
  path = "./out",
  pattern = "*mascot*",
  full.names = TRUE
) 
tip_state = data.frame()
# remove all files that don't end in log
log <- log[endsWith(log, "log")]
# remove all files that don't contain "it" and rep0
log <- log[grepl("it", log) & grepl("rep0", log)]
for (i in seq(1, length(log))){
# for (i in seq(1, 1)){
    # get the number in the file name directly after it
    run = gsub(".*it", "", log[i])
    run = strsplit(run, split="_")[[1]][[1]]
    # read in the log file as a dataframe and discard the first 10% of rows
    df1 <- read.table(log[i], header = TRUE, check.names=F)
    df3 <- read.table(gsub("rep0", "rep1", log[i]), header = TRUE, check.names=F)
    df4 <- read.table(gsub("rep0", "rep2", log[i]), header = TRUE, check.names=F)
    # combine the three dataframes after removing the first 10% of rows from each df
    df <- rbind(df1[-c(1:round(nrow(df1)*0.1)),], df3[-c(1:round(nrow(df3)*0.1)),], df4[-c(1:round(nrow(df4)*0.1)),])
    # get the parameter names
    param_names <- colnames(df)
    # loop over the parameter names
    for (param in param_names){
        # if the param name contains a | then it is for a sample that was tip sampled
        if (grepl("sampledState", param)){
          true_state = strsplit(param, "\\|")[[1]][4]
          # replace .sampledState with nothing
          true_state = gsub(".sampledState", "", true_state)
          state = which(states == true_state)
          # get the values for the parameter
          values <- df[[param]]
          # add the values to the dataframe
          tip_state <- rbind(tip_state, 
              data.frame(name = gsub(".sampledState","",param), state = state, value.mascot = mean(values==(state-1)), run = as.numeric(run))
          )       
        }
    }
}

# make an additional row in tip_state for the dta estimates
tip_state$value.dta = NA

# read in all the dta *.trees files
trees <- list.files(  # nolint
  path = "./out",
  pattern = "*.trees",
  full.names = TRUE
)
# remove all files that contain "it" and rep0 or dta
trees <- trees[grepl("it", trees) & grepl("rep0", trees) & grepl("dta", trees)]
# loop over the files
for (i in seq(1, length(trees))){
    print(i)
    # get the number in the file name directly after it
    run = gsub(".*it", "", trees[i])
    run = strsplit(run, split="_")[[1]][[1]]

    # read in the tree files line by line, discard any lines that don't start with tree STATE_
    tree1 <- readLines(trees[i])
    # get all the lines after the line that contains Translate and before tree_STATE_0 and save them into a new vector
    tip_names = data.frame()
    start = FALSE
    for (j in seq(1, length(tree1))){
      if (grepl("Translate", tree1[[j]])){
        start = TRUE
        j=j+1;
      }
      if (start & grepl(";", tree1[[j]])){
        break
      }
      if (start){
        tmp = strsplit(tree1[[j]], split=" ")[[1]]
        tip_names = rbind(tip_names, data.frame(number=as.numeric(tmp[[length(tmp)-1]]), name = gsub(",","", tmp[[length(tmp)]])))
      }
    }

    tree1 <- tree1[grepl("tree STATE_", tree1)]
    # remove the first 10% of lines
    tree1 <- tree1[-c(1:round(length(tree1)*0.1))]
    # do the same for the files ending in rep1 and rep2
    tree2 <- readLines(gsub("rep0", "rep1", trees[i]))
    tree2 <- tree2[grepl("tree STATE_", tree2)]
    tree2 <- tree2[-c(1:round(length(tree2)*0.1))]
    tree3 <- readLines(gsub("rep0", "rep2", trees[i]))
    tree3 <- tree3[grepl("tree STATE_", tree3)]
    tree3 <- tree3[-c(1:round(length(tree3)*0.1))]
    # combine the three trees
    tree <- c(tree1, tree2, tree3)
    # get all the unique names in tip_state for which run = run
    names = unique(tip_state$name[tip_state$run == run])
    # print(names)
    # print("")
    # print("")
    # print("")
    # print(i)
    
    # start a vector with names columns and length(trees) rows
    tip_state2 = matrix(NA, nrow=length(tree), ncol=length(names))

    # loop over the trees
    for (j in seq(1, length(tree))){
        # get all the positions of "]" in tree[[j]]
        positions = gregexpr("]", tree[[j]])[[1]]
        # loop over the names
        for (k in seq(1, length(names))){
            # in tree[[j]] find the position of what starts with a "," or "(" and then tip_names[tip_names$name==names[[k]], "number"] + "["
            block_start = gregexpr(paste0("[(,]", tip_names[tip_names$name==names[[k]], "number"], "\\["), tree[[j]])[[1]]        
            # block_end is the first value in positions > that block_start
            block_end = positions[positions > block_start][1]
            # get the block
            # print("")
            block = substr(tree[[j]], block_start, block_end)
            # print(paste(block, states[as.numeric(tip_state[tip_state$name==names[[k]], "state"][[1]])]))
            # if the block and names[[k]] both contain the same element from states, then add a 1 to the vector, otherwise, add a 0
            if (grepl(states[as.numeric(tip_state[tip_state$name==names[[k]], "state"])][[1]], block)){
                tip_state2[j,k] = 1
                # print(1)
            } else {
                tip_state2[j,k] = 0
                # print(0)
            }
        }
    }
    # add the mean support for the true value to tip_state
    for (k in seq(1, length(names))){
        # get the index in tip_state for which name = names[[k]] and run = run
        index = which(tip_state$name == names[[k]] & tip_state$run == run)
        # add the value.dta for this index
        tip_state$value.dta[index] = mean(tip_state2[,k])
    }
    print('done')
}
# for each of the states, get the frequency of the state in tip_names$name
count = c()
for (i in seq(1, length(states))){
  count = c(count, sum(grepl(states[[i]], tip_names$name)))
}

# make count a dataframe, such that it can be plotted with tip_state
count = data.frame(statename = states, freq = count/sum(count))

tip_state$statename = states[tip_state$state]

tip_state$runname = factor(as.character(tip_state$run), levels = as.character(seq(1,max(tip_state$run))))

p1 = ggplot(tip_state, aes(x = value.mascot, y = value.dta, color=runname)) +
  geom_point() +
  # add a dashed line
  geom_hline(data=count, aes(yintercept=freq), color="black", linetype="dashed")+
  geom_vline(data=count, aes(xintercept=freq), color="black", linetype="dashed")+
  facet_wrap(statename~., ncol=2)+
  geom_abline(intercept = 0, slope = 1, color = "red") +  # Add diagonal line
  scale_color_discrete_qualitative(name="Iteration")+
  xlab("posterior support for true tip state using MASCOT-Skyline")+
  ylab("posterior support for true tip state using DTA")+
  theme_minimal()
plot(p1)
ggsave(plot=p1, filename = "../../../MascotSkyline-Text/Figures/sars_comp_tipsampling.pdf", height=6, width=8)


library(tidyverse)

# Reshaping the data to long format
long_tip_state <- tip_state %>% 
  gather(key = "method", value = "value", value.mascot, value.dta) %>%
  mutate(method = ifelse(method == "value.mascot", "MASCOT-Skyline", "DTA"))

p2 = ggplot(long_tip_state) +
  geom_vline(data=count, aes(xintercept=freq), color="black", linetype="dashed") +
  geom_histogram(aes(x=value, y=..density.., fill=method, color=method), bins=10, alpha=0.6, position="dodge") +
  geom_density(aes(x=value, fill=method, color=method), alpha=0.2) +
  scale_fill_manual(name="", values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510")) +
  scale_color_manual(name="", values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510")) +
  facet_wrap(statename~., ncol=2) +
  theme_minimal()
plot(p2)
ggsave(plot=p2, filename = "../../../MascotSkyline-Text/Figures/sars_compdist_tipsampling.pdf", height=6, width=8)


p3 = ggplot(tip_state) +
  geom_vline(aes(xintercept=0), color="black", linetype="dashed")+
  stat_density(aes(x=value.mascot-value.dta), adjust=0.5,alpha=0.5, fill="black")+
  facet_wrap(statename~., ncol=2)+
  theme_minimal() + xlab("difference in posterior support for true state")
plot(p3)

ggsave(plot=p3, filename = "../../../MascotSkyline-Text/Figures/sars_md_tipsampling.pdf", height=6, width=8)


# get all log files  that contain it
log <- list.files(  # nolint
  path = "./out",
  pattern = "*log",
  full.names = TRUE
)
# remove all files that do contain "it" and rep0
log <- log[grepl("it", log) & grepl("rep0", log)]
mig=data.frame()

# loop over the files
for (i in seq(1, length(log))){
    # get the number in the file name directly after it
    run = gsub(".*it", "", log[i])
    run = strsplit(run, split="_")[[1]][[1]]
    # read in the log file as a dataframe and discard the first 10% of rows
    df1 <- read.table(log[i], header = TRUE, check.names=F)
    df3 <- read.table(gsub("rep0", "rep1", log[i]), header = TRUE, check.names=F)
    df4 <- read.table(gsub("rep0", "rep2", log[i]), header = TRUE, check.names=F)
    # combine the three dataframes after removing the first 10% of rows from each df
    df <- rbind(df1[-c(1:round(nrow(df1)*0.1)),], df3[-c(1:round(nrow(df3)*0.1)),], df4[-c(1:round(nrow(df4)*0.1)),])
    # get the parameter names
    param_names <- colnames(df)
    # loop over the parameter names
    for (param in param_names){
        # if the log file name contains mascot and the param starts with f_migration
        if (grepl("mascot", log[i]) & grepl("f_migration", param)){
            # get the from and to values by splitting param on "_" and "."
            from_to = strsplit(param, split="[_.]")[[1]] 
            # get the values for the parameter
            values <- df[[param]]
            # get the 95% HPD 
            hpd = HPDinterval(as.mcmc(values), prob=0.95)
            # add the values to the dataframe
            mig <- rbind(mig, 
                data.frame(value = median(values), from=from_to[3], to=from_to[5], hpd.lower=hpd[1], hpd.upper=hpd[2], run = as.numeric(run), method="MASCOT-Skyline")
            )       
        }
        # do the same if the log file name contains dta and the param starts with c_migrationEvents.
        if (grepl("dta", log[i]) & grepl("geo.rates.", param)){
          from_to = strsplit(param, split="[_.]")[[1]] 
          values <- df[[param]]
          hpd = HPDinterval(as.mcmc(values), prob=0.95)
          mig <- rbind(mig, 
                       data.frame(value = median(values), from=from_to[3], to=from_to[4], hpd.lower=hpd[1], hpd.upper=hpd[2], run = as.numeric(run), method="DTA")
          )       
        }
    }
}

# plot the results of the 10 runs and MASCOT-Skyline and DTA next to each other using a facet_grid(from~to, ncol=4)
# for the x values, use the run number -0.15 if Mascot and +0.15 if DTA
# for the y values, use the value and the upper and lower
# for the color, use the method
mig$x = mig$run
mig$x[mig$method=="DTA"] = mig$x[mig$method=="DTA"] + 0.15
mig$x[mig$method=="MASCOT-Skyline"] = mig$x[mig$method=="MASCOT-Skyline"] - 0.15
p4 <- ggplot(mig, aes(x=x,y=value, color=method)) +
  geom_point() +
  geom_errorbar(aes(ymin=hpd.lower, ymax=hpd.upper)) +
  scale_y_log10()+
  facet_grid(from~to, scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_color_manual(name="", values = c("MASCOT-Skyline" = "#228B45", "DTA" = "#8B4510")) +
  xlab("")+ylab("migration rate")
plot(p4)
ggsave(plot=p4, filename = "../../../MascotSkyline-Text/Figures/sars_migrationrates_tipsampling.pdf", height=9, width=12)
