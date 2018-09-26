###
# This is a script that reads matlab EEG data, does some light pre-processing,
# and attempts to predict the outcome class based on Common Spatial Patterns
# method (Ramoser et al., 2000) and a Linear Discriminant Analysis classifier.
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
library(ggplot2)
library(MASS)
library(caret)
source("CSP_utils.R")

# make reproducible
set.seed(423)

# list data files, there's one file per subject
filelist <- list.files("data")
n_subj = length(filelist)

# set parameter list for k to go over
# k is the number of eigenvector pairs to be used after CSP
k_list <- 1:5

# initialize results data frame
results_df <- data.frame(subj = character(n_subj*length(k_list)), 
                         k = numeric(n_subj*length(k_list)), 
                         acc = numeric(n_subj*length(k_list)),
                         stringsAsFactors = FALSE)

# loop over all subjects
for (i in 1:n_subj) {
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
    
    # get electrode labels for reference (not used in the analysis)
    el_labels <- unlist(raw_dat_subj[[1]][[5]])
    
    # now get the outcomes
    outcome <- unlist(raw_dat_subj[[1]][[11]]) # ugly vector, need only last values
    outcome <- as.numeric(outcome[(length(outcome)-n_trials+1) : length(outcome)])
    
    # indices for classes 0 (didn't remember) and 1 or 2 (somewhat or fully remembered)
    cl_1_idx <- which(outcome == 0)
    cl_2_idx <- which(outcome == 1 | outcome == 2)

    # subset the data into classes 1 (0 in the data) and 2 (1 and 2 in the data)
    dat_subj_cl_1 <- dat_subj[cl_1_idx, , ]
    dat_subj_cl_2 <- dat_subj[cl_2_idx, , ]
    
    # calculate average spatial covariance matrices for classes 1 and 2
    cov_subj_cl_1 <- avg_cov_mat(dat_subj_cl_1)
    cov_subj_cl_2 <- avg_cov_mat(dat_subj_cl_2)
    
    # get the projection matrix
    W <- proj_mat(cov_subj_cl_1, cov_subj_cl_2)
    
    # calculate base accuracy for this subject, i.e. the accuracy one would 
    # get naively when always predicting the more frequent class
    # -- used previously without upsampling
    # base_acc <- max(nrow(dat_subj_cl_1), nrow(dat_subj_cl_2)) / 
    #     (nrow(dat_subj_cl_1) + nrow(dat_subj_cl_2))
    
    # loop over different k's 
    for (k in k_list){
        # calculate feature matrices taking k pairs of projection vectors (first and last)
        feat_mat_cl_1 <- feature_mat(dat_subj_cl_1, W, k)
        feat_mat_cl_2 <- feature_mat(dat_subj_cl_2, W, k)
        
        # create new data frame for this subject
        df_subj <- data.frame(v1 = c(feat_mat_cl_1[, 1], feat_mat_cl_2[, 1]),
                              v2 = c(feat_mat_cl_1[, 2], feat_mat_cl_2[, 2]),
                              class = factor(c(rep(1, nrow(feat_mat_cl_1)), 
                                        rep(2, nrow(feat_mat_cl_2)))))
        
        # upsample the data to create balanced classes
        # rows from the less frequent class will be sampled with substitution
        df_subj <- upSample(df_subj[1:2], df_subj$class, yname = "class")
        
        # define the training parameters - 10-fold CV
        fit_control <- trainControl(method = "cv", number = 10)
        
        # run lda, get predictions
        mod <- train(class ~ ., data = df_subj,
                     method = "lda", trControl = fit_control)

        # get accuracy from CV
        acc <- mod$results$Accuracy
        
        # save subject, k, and acc to the results data frame
        subj <- strsplit(filelist[i], "_")[[1]][1]
        results_df[k+(length(k_list)*(i-1)), "subj"] <- subj
        results_df[k+(length(k_list)*(i-1)), "k"] <- k
        results_df[k+(length(k_list)*(i-1)), "acc"] <- acc
        #results_df[k+(length(k_list)*(i-1)), "base_acc"] <- base_acc
    }
}

# plot the results
acc_plot <- ggplot(data = results_df, aes(x = k, y = acc)) +
    geom_col(position = "dodge") +
    geom_hline(aes(yintercept = 0.5), col = "red", size = 1) +
    facet_wrap(~ subj) +
    labs(x = "# of CSP eigenvector pairs, k",
         y = "accuracy",
         title = "Classification accuracy for 15 subjects (3101-3116)",
         subtitle = "(the red horizontal line represents baseline accuracy of 0.5)",
         caption = "based on upsampled data")

acc_plot

# save plot to disk
ggsave("CSP_results.png", plot = acc_plot, width = 12, height = 8, dpi = 300)