###
# NOTE: THIS SCRIPT TAKES MULTIPLE HOURS TO RUN!
# (choosing fewer and smaller m's would shorten run time)
# 
# This script implements pdc for 15 subjects over a range of embedding 
# dimensions and plots MDS plots for each run.
# see https://cran.r-project.org/web/packages/pdc/pdc.pdf
# 
# Basic purpose here is illustration.
#
# The plots and the clusterings are NOT saved to disk.
# 
###

library(R.matlab)
library(pdc)

filelist <- list.files("data")
n_subj = length(filelist)

# set the list of embedding dimensions (m) to try
m_list = 3:7

for (i in 1:n_subj){
    # load raw data from matlab, will be a very ugly very nested list
    raw_dat_subj <- readMat(paste0("data/", filelist[i]))
    
    # now get the data as a 3D-array with (trials x channels x timepoints)
    dat_subj <- raw_dat_subj[[1]][[6]]
    
    # get number of trials for this subject / condition
    n_trials <- dim(dat_subj)[1]
    
    # reduce time-window to between 0 (stimulus onset) and 2.5 seconds
    # first find the indices of the desired timeframe
    time <- raw_dat_subj[[1]][[3]]
    tw_start <- which(time == 0)
    tw_end <- which(time == 2.5)
    
    # then slice the array
    dat_subj <- dat_subj[,,tw_start:tw_end]
    
    # now get the outcomes
    outcome <- unlist(raw_dat_subj[[1]][[11]]) # ugly vector, need only last values
    labels <- as.numeric(outcome[(length(outcome)-n_trials+1) : length(outcome)])
    
    # subject id
    subj <- strsplit(filelist[i], "_")[[1]][1]
    
    # prepare for pdc - change order of array dimensions
    dat_subj <- aperm(dat_subj, c(3,1,2)) # time X trials/stimuli X channel
    
    # true label colors
    true_col <- labels
    true_col[true_col == 0] = 'red'
    true_col[true_col == 1 | true_col == 2] = 'blue'
    
    # clustering over the list of embedding dimensions, m
    # a series of plots will be created, but not saved to disk
    for (m in m_list){
        clustering <- pdc::pdclust(dat_subj, m=m)
        
        mdsPlot(clustering, labels, true_col)
        title(sub = paste0("subject ", subj, ", m = ", m))
    }
    
}

