### This script contains functions for Common Spatial Pattern (CSP)
### dimensionality reduction
### see http://www.dsp.toronto.edu/~haiping/Publication/R-CSP_EMBC2009.pdf

# calculate average spatial covariance matrix for a class
# over all trials for a subject
avg_cov_mat <- function(dat_class){
    # dat_class is trial matrix MxNxT, M trials, N channels, and T timepoints
    n_trials <- dim(dat_class)[1]
    sum_mat <- 0 # initialize zero matrix
    for (i in 1:n_trials){
        m <- dat_class[i, , ] # this trial
        S <- m %*% t(m) / sum(diag(m %*% t(m))) # covariance matrix
        sum_mat <- sum_mat + S # add up
    }
    avg_mat <- sum_mat / n_trials
    avg_mat
}

# # calculate projection matrix from CSP, take k eigenvector pairs
# proj_mat <- function(cov_mat_1, cov_mat_2, k){
#     decomp_mat <- solve(cov_mat_2) %*% cov_mat_1
#     ev <- eigen(decomp_mat)
#     W <- cbind(ev$vectors[, 1:k], ev$vectors[, (ncol(ev$vectors)-k+1):ncol(ev$vectors)])
#     W
# }
# 
# # calculate new feature matrix based on CSP for a class
# # over all trials for a subject
# feature_mat <- function(dat_class, proj_mat){
#     # dat_class is trial matrix MxNxT, M trials, N channels, and T timepoints
#     # proj_mat is projection matrix consisting of the first and last eigenvector
#     n_trials <- dim(dat_class)[1]
#     feat_mat <- matrix(NA, nrow = n_trials, ncol = ncol(proj_mat)) # initialize zero matrix
#     for (i in 1:n_trials){
#         m <- dat_class[i, , ] # this trial
#         Z <- t(proj_mat) %*% m # project
#         # calculate row-wise variance 
#         var_Z <- 0 
#         for (j in 1:nrow(Z)){
#             var_Z <- var_Z + var(Z[j, ])
#         }
#         #calculate feature vector
#         feat_vec <- rep(NA, nrow(Z))
#         for (k in 1:nrow(Z)){
#             y_q <- log( var(Z[k, ]) / var_Z)
#             feat_vec[k] <- y_q
#         }
#         feat_mat[i, ] <- feat_vec # add to feature matrix
#     }
#     feat_mat
# }


##########################################################

# calculate projection matrix from CSP, take k eigenvector pairs
proj_mat <- function(avg_cov_mat_1, avg_cov_mat_2){
    C <- avg_cov_mat_1 + avg_cov_mat_2 # composite matrix
    ev <- eigen(C) # get eigenvalues and eigenvectors
    P <- sqrt(solve(diag(ev$values))) %*% t(ev$vectors) # whitening
    S_1 <- P %*% avg_cov_mat_1 %*% t(P) # transform e.g. cov mat for class 1
    B <- eigen(S_1)$vectors # get eigenvectors we're looking for
    W <- t(t(B) %*% P)
    W
}

# calculate new feature matrix based on CSP for a class
# over all trials for a subject
feature_mat <- function(dat_class, proj_mat, k){
    # dat_class is trial matrix MxNxT, M trials, N channels, and T timepoints
    # proj_mat is projection matrix consisting of the first and last eigenvector
    n_trials <- dim(dat_class)[1]
    feat_mat <- matrix(NA, nrow = n_trials, ncol = ncol(proj_mat)) # initialize zero matrix
    for (i in 1:n_trials){
        m <- dat_class[i, , ] # this trial
        Z <- proj_mat %*% m # decomposition (mapping) for this trial
        Z <- rbind(Z[1:k, ], Z[nrow(Z)-k:nrow(Z), ]) # get only k first and last rows
        # calculate row-wise variance 
        var_Z <- 0 
        for (j in 1:nrow(Z)){
            var_Z <- var_Z + var(Z[j, ])
        }
        #calculate feature vector
        feat_vec <- rep(NA, nrow(Z))
        for (k in 1:nrow(Z)){
            y_q <- log( var(Z[k, ]) / var_Z)
            feat_vec[k] <- y_q
        }
        feat_mat[i, ] <- feat_vec # add to feature matrix
    }
    feat_mat
}