library(R.matlab)
library(ggplot2)
library(MASS)
source("utils.R")

###
# This is a very hacky script to make the Matlab EEG data somewhat tidy and 
# analyzable in R.
#
# The hardcoded list indices here come from me viewing the "struct" object 
# in Matlab.
###

set.seed(42)

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
    
    #### maybe main.R should start here?
    # indices for classes 0 (didn't remember) and 1 or 2 (somewhat or fully remembered)
    cl_1_idx <- which(outcome == 0)
    cl_2_idx <- which(outcome == 1 | outcome == 2)
    
    # subset the data into classes 1 and 2
    dat_subj_cl_1 <- dat_subj[cl_1_idx, , ]
    dat_subj_cl_2 <- dat_subj[cl_2_idx, , ]
    
    # calculate average spatial covariance matrices for classes 1 and 2
    cov_subj_cl_1 <- avg_cov_mat(dat_subj_cl_1)
    cov_subj_cl_2 <- avg_cov_mat(dat_subj_cl_2)
    
    # loop over 5 different k's 
    for (k in k_list){
        # get the projection matrix using k pairs of eigenvectors (first and last)
        W <- proj_mat(cov_subj_cl_1, cov_subj_cl_2)
        
        # TBD: 1) http://www.dsp.toronto.edu/~haiping/Publication/R-CSP_EMBC2009.pdf formula 13
        # implement that over all trials, getting a matrix of new features for each class
        # 2) put both classes into one data frame with class labels
        # 3) split into train / test set
        # 4) run MASS::lda over this data frame, 
        # see https://tgmstat.wordpress.com/2014/01/15/computing-and-visualizing-lda-in-r/
        
        # calculate feature matrices
        feat_mat_cl_1 <- feature_mat(dat_subj_cl_1, W, k)
        feat_mat_cl_2 <- feature_mat(dat_subj_cl_2, W, k)
        
        # create new data frame for this subject
        df_subj <- data.frame(v1 = c(feat_mat_cl_1[, 1], feat_mat_cl_2[, 1]),
                              v2 = c(feat_mat_cl_1[, 2], feat_mat_cl_2[, 2]),
                              class = c(rep(1, nrow(feat_mat_cl_1)), 
                                        rep(2, nrow(feat_mat_cl_2))))
        
        # shuffle it row-wise
        df_subj <- df_subj[sample(nrow(df_subj)), ]
        
        # split into 80% train, 20% test
        df_subj_train <- df_subj[1:(round(0.8*nrow(df_subj))), ]
        df_subj_test <- df_subj[(round(0.8*nrow(df_subj))+1):nrow(df_subj), ]
        
        # run lda, get predictions
        mod <- MASS::lda(class ~ ., data = df_subj_train)
        preds <- predict(mod, df_subj_test)
        #preds_cl <- ifelse(preds$posterior[, 1] > 0.5, 1, 2)
        
        # calculate accuracy
        acc <- mean(as.numeric(preds$class) == df_subj_test["class"])
        
        # save subject, k, and acc to the results data frame
        subj <- strsplit(filelist[i], "_")[[1]][1]
        results_df[k+(length(k_list)*(i-1)), "subj"] <- subj
        results_df[k+(length(k_list)*(i-1)), "k"] <- k
        results_df[k+(length(k_list)*(i-1)), "acc"] <- acc
    }
}

acc_plot <- ggplot(data = results_df, aes(x = subj, y = acc, fill = k)) +
    facet_grid(k ~ .) +
    geom_col()
acc_plot

acc_plot2 <- ggplot(data = results_df, aes(x = subj, y = acc, fill = as.factor(k))) +
    geom_col(position = "dodge") +
    scale_fill_brewer("k", palette = "Set1") +
    geom_hline(aes(yintercept = mean(acc[k == 1])), col = "red", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 2])), col = "blue", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 3])), col = "green", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 4])), col = "violet", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 5])), col = "orange", size = 1)

acc_plot2

acc_plot3 <- ggplot(data = results_df, aes(x = as.numeric(subj), y = acc, col = as.factor(k))) +
    geom_line() +
    scale_color_brewer("k", palette = "Set1") +
    geom_hline(aes(yintercept = mean(acc[k == 1])), col = "red", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 2])), col = "blue", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 3])), col = "green", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 4])), col = "violet", size = 1) +
    geom_hline(aes(yintercept = mean(acc[k == 5])), col = "orange", size = 1)

acc_plot3
