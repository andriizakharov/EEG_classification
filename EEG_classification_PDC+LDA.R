###
# This is a script that reads matlab EEG data, does some light pre-processing,
# and attempts to predict the outcome class based on Permutation Distribution 
# Clustering method (Brandmaier 2015) and a Linear Discriminant Analysis classifier.
# A permutation distribution with a given embedding dimension is calculated for
# each EEG channel, and classification features are created based on the
# per-channel entropy of this distribution.
#
# The data is upsampled to create balanced classes.
# 
# Note that although there were three outcome classes in the data (0,1,2), here
# the classification is binary, i.e. 0 vs 1+2. The results don't change
# much if run as 0 vs 2 or 0 vs 1. 
#
# The hardcoded list indices here come from viewing the "struct" object 
# in Matlab.
#
# The script is not parallelized and takes ~10 mins to run on my i5-6300u machine.
###

library(R.matlab)
library(pdc)
library(ggplot2)
source("PDC_utils.R")

# make reproducible
set.seed(423)

# list data files, there's one file per subject
filelist <- list.files("data")
n_subj = length(filelist)

# set the list of embedding dimensions (m) to try
m_list = 3:7

# initialize results data frame
results_df <- data.frame()

for (i in 1:n_subj){
    # load raw data from matlab, will be a very ugly very nested list
    raw_dat_subj <- readMat(paste0("data/", filelist[i]))
    
    # now get the data as a 3D-array with (trials x channels x timepoints)
    dat_subj <- raw_dat_subj[[1]][[6]]
    
    # get electrode labels for reference
    el_labels <- unlist(raw_dat_subj[[1]][[5]])
    
    # get number of trials and channels for this subject / condition
    n_trials <- dim(dat_subj)[1]
    n_channels <- dim(dat_subj)[2]
    
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
    
    # put classes 1 and 2 together to make the problem binary
    labels[labels == 2] <- 1
    
    # subject id
    subj <- strsplit(filelist[i], "_")[[1]][1]
    
    # loop over m's
    for (m in m_list){
        # initialize a feature matrix, features will be per-channel entropy
        # over all trials
        feature_mat <- matrix(NA, nrow = n_trials, ncol = n_channels)
        
        # for each channel for each trial
        for (t in 1:n_trials){
            for (c in 1:n_channels){
                # calculate the permutation distribution entropy
                feature_mat[t, c] <- codebook.entropy(dat_subj[t, c, ], m = m)
            }
        }
        
        # the data frame for classification
        df_subj <- data.frame(feature_mat, factor(labels))
        names(df_subj) <- c(el_labels, "class")
        
        # upsample the data to create balanced classes
        # rows from the less frequent class will be sampled with substitution
        df_subj <- upSample(df_subj[, -ncol(df_subj)], df_subj$class, yname = "class")
        
        # define the training parameters - 10-fold CV
        fit_control <- trainControl(method = "cv", number = 10)
        
        # run lda, get predictions
        mod <- train(class ~ ., data = df_subj,
                     method = "lda", trControl = fit_control)
        
        # get accuracy from CV
        acc <- mod$results$Accuracy
        
        # save subject, m, and acc to the results data frame
        results_df <- rbind.data.frame(results_df,
                            list(subj, m, acc), stringsAsFactors = FALSE)
    }
}

names(results_df) <- c("subj", "m", "acc")

# plot the results
acc_plot <- ggplot(data = results_df, aes(x = m, y = acc)) +
    geom_col(position = "dodge") +
    geom_hline(aes(yintercept = 0.5), col = "red", size = 1) +
    facet_wrap(~ subj) +
    labs(x = "Embedding dimension, m",
         y = "accuracy",
         title = "Classification accuracy for 15 subjects (3101-3116)",
         subtitle = "(the red horizontal line represents baseline accuracy of 0.5)",
         caption = "based on upsampled data")

acc_plot

# save plot to disk
ggsave("PDC_results.png", plot = acc_plot, width = 12, height = 8, dpi = 300)
    