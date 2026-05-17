# Sparse Reduced Rank Regression Simulation Study
# Author: Lucy McLaughlin
# This script performs simulation experiments comparing
# OLS, RRR, and SRRR in multivariate linear regression settings.
# The study evaluates prediction error, estimation accuracy,
# and variable selection performance under varying sparsity
# and correlation structures.

# Required package:
# install.packages("rrpack")

library(rrpack) 


one_run <- function(n,p,p0, q,q0, nrank, s2n,rho_X,rho_E, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  sim <- rrr.sim2( ##rrr.sim2 is most similar to the model in Chen and Huang (2012)
    n= n, #sample size 
    p=p, #number of predictors
    p0 = p0, #number of relevant predictors (where sparsity can be induced)
    q = q, #number of responses
    q0= q0, #number of relevant responses (In this study: q0 = q)
    nrank = nrank, #model rank
    s2n = s2n, #signal to noise ratio 
    sigma = NULL, # error variance (not specified)
    rho_X = rho_X, #correlation parameter for generating predictor variables
    rho_E = rho_E #correlation parameter for generating random errors
  )
  
  C <- sim$C #coefficient matrix
  X <- sim$X #predictor matrix
  Y <- sim$Y #response matrix
 
##############   TESTING VS TRAINING SET   ##############
   
#the following splits the data into 70% training and 30% testing sets  
  n_obs <- nrow(X)  
  n_train <- floor(0.7*n_obs)
  train_id <- sample(1:n_obs, n_train) 
  test_id <- setdiff(1:n_obs, train_id) 

  X_train <- X[train_id,] #predictor matrix for the training set
  Y_train <- Y[train_id,] #response matrix for the training set
  
  X_test <- X[test_id,] #predictor matrix for the testing set
  Y_test <- Y[test_id, ] #response matrix for the testing set 
  
##############   OLS   ##############
  
#the following computes the OLS and Moore-Penrose Pseudoinverse estimate
#for the coefficient matrix
#the mean squared prediction error, signal prediction error and estimation error
#is also defined
#the if-else ensures that OLS is used for the scenarios where p < n,
#and the Moore-Penrose Pseudoinverse is used for p > n
  if(p <= n_train){
    
    C_ols <- solve(t(X_train)%*%X_train)%*%t(X_train)%*%Y_train
    
    MSPE_ols <- mean(((X_test %*% C_ols) - Y_test)^2)
    SP_ols <- mean((X_test %*% C - X_test %*% C_ols)^2)
    Est_ols <- norm(C_ols - C, "F")^2
  }else{
    C_ols <- MASS::ginv(X_train) %*% Y_train
    
    MSPE_ols <- mean(((X_test %*% C_ols) - Y_test)^2)
    SP_ols <- mean((X_test %*% C - X_test %*% C_ols)^2)
    Est_ols <- norm(C_ols - C, "F")^2
  }
  
##############   RRR   ##############  
  
  fit_rrr <- rrr.fit(Y_train,X_train, nrank = nrank) #fits RRR for the specified rank (nrank)
  C_rrr <- coef(fit_rrr) #uses RRR to estimate the coefficient matrix

#the following defines the mean squared prediction error, signal prediction error 
#and estimation error for RRR  
  MSPE_rrr <- mean(((X_test %*% C_rrr)- Y_test)^2)
  SP_rrr <- mean((X_test %*% C - X_test %*% C_rrr)^2)
  Est_rrr <- norm(C_rrr - C, "F")^2

##############   SRRR   ##############
    
#the same process is repeated for SRRR:  
  fit_srrr <- srrr(Y_train, X_train, nrank = nrank, method = "glasso")
#method = "glasso" ensures group lasso penalty is used as the method to induce sparsity
#induces row-wise sparsity  
  C_srrr <- coef(fit_srrr)
  MSPE_srrr <- mean(((X_test %*% C_srrr)- Y_test)^2)
  SP_srrr <- mean((X_test %*% C - X_test %*% C_srrr)^2)
  Est_srrr <- norm(C_srrr - C, "F")^2
  
##############   SPARSITY MEASURES  ##############
  
#the following calculates the number of rows that are zero and non zero
#for an analysis of Sparse Reduced Rank Regression's variable selection performance
  zero_rows_srrr <- sum(rowSums(abs(C_srrr)) < 1e-6)
  nonzero_rows_srrr <- nrow(C_srrr) - zero_rows_srrr
  
  true_rows <- which(rowSums(abs(C)) > 0)
  selected_rows <- which(rowSums(abs(C_srrr)) > 1e-6)
#the True Positive Rate (TPR) and False Positive Rate (FPR) are defined below
#TPR is used to measure SRRR's  ability to correctly identify relevant predictors
  TPR <- length(intersect(true_rows,selected_rows)) / p0
#FPR is used to measure how often SRRR includes irrelevant predictors 
  FPR <- if((p-p0) >0){
    length(setdiff(selected_rows, true_rows)) / (p-p0)
  }else{
    NA
  }
#setdiff(x,y) finds rows in x that aren't in y 
  
  return(c(MSPE_ols = MSPE_ols,
           SP_ols = SP_ols,
           Est_ols = Est_ols,
           MSPE_rrr = MSPE_rrr,
           SP_rrr = SP_rrr,
           Est_rrr = Est_rrr,
           MSPE_srrr = MSPE_srrr,
           SP_srrr = SP_srrr,
           Est_srrr = Est_srrr,
           zero_srrr = zero_rows_srrr,
           nonzero_srrr = nonzero_rows_srrr,
           True_pr =TPR,
           False_pr = FPR))
}

############  SCENARIO SET UP  ###########


run_scenario <- function(scenario_name,nrep, n, p, p0,q,q0,nrank, s2n, 
                         rho_X, rho_E){
  
  results <- replicate( 
    nrep, 
    one_run(
      n= n,
      p= p,
      p0 = p0,
      q = q,
      q0 = q0,
      nrank = nrank,
      s2n = s2n, 
      rho_X = rho_X, 
      rho_E = rho_E
    )
  )
  
  nonzero_rows_results <- results["nonzero_srrr", ]
  
  #boxplot of the log of the Mean Squared Error

  boxplot(
    data.frame(
      OLS = results["MSPE_ols",],
      RRR = results["MSPE_rrr",],
      SRRR = results["MSPE_srrr",]
    ),
    main = paste("MSPE\n", scenario_name),
    ylab = "MSPE (Log scale)",
    col = c("pink","skyblue","lavender"),
    log = "y",
    cex.lab = 1.5,
    cex.axis = 1.3,
    cex.main = 1.6
  )

  #boxplot of the log of the estimation error
  boxplot(
    data.frame(
      OLS = results["Est_ols",],
      RRR = results["Est_rrr",],
      SRRR = results["Est_srrr",]
    ),
    main = paste("Estimation Error\n", scenario_name),
    ylab = "Estimation Error (Log scale)",
    col = c("pink","skyblue","lavender"),
    log = "y",
    cex.lab = 1.5,
    cex.axis = 1.3,
    cex.main = 1.6
  )

  #boxplot of the log of the Signal prediction error:

  boxplot(
    data.frame(
      OLS = results["SP_ols",],
      RRR = results["SP_rrr",],
      SRRR = results["SP_srrr",]
    ),
    main = paste("Signal Prediction Error\n", scenario_name),
    ylab = "SPE (Log scale)",
    col = c("pink","skyblue","lavender"),
    log = "y",
    cex.lab = 1.5,
    cex.axis = 1.3,
    cex.main = 1.6
  )

  
#the following code creates a histogram showing the number of relevant predictors
#the SRRR selects in each simulation run compared to the true number of relevant
#predictors
  x_limits <- range(c(nonzero_rows_results, p0)) + c(-2,2)



  hist(nonzero_rows_results,
       main = paste("Predictors Selected by SRRR\n", scenario_name),
       xlab = "Selected Predictors",
       xlim = x_limits,
       border = "white",
       cex.main = 0.9)

  abline(v = p0, col = "red", lwd = 2, lty = 2)
  text(p0, par("usr")[4]*0.9, expression(p[0]), col="red", pos=4)

#the following code defines the mean, median and interquartile range (IQR) for
#the errors calculated in each simulation scenario for each method
  avg <- rowMeans(results, na.rm = TRUE)
  avg[is.nan(avg)] <- NA
  
  median_vals <- apply(results, 1, median, na.rm = TRUE)
  
  q1_vals <- apply(results, 1, function(x) quantile(x, 0.25, na.rm = TRUE))
  q3_vals <- apply(results, 1, function(x) quantile(x, 0.75, na.rm = TRUE))
  IQR <- q3_vals - q1_vals
  
  results_table <- data.frame(
    Scenario = scenario_name,
    Method = c("OLS","RRR","SRRR"),
    
    MSPE_mean = unname(c(avg["MSPE_ols"], avg["MSPE_rrr"], avg["MSPE_srrr"])),
    MSPE_median = unname(c(median_vals["MSPE_ols"], median_vals["MSPE_rrr"], median_vals["MSPE_srrr"])),
    MSPE_IQR = unname(c(IQR["MSPE_ols"], IQR["MSPE_rrr"], IQR["MSPE_srrr"])),

    Signal_PE_mean = unname(c(avg["SP_ols"], avg["SP_rrr"], avg["SP_srrr"])),
    Signal_PE_median = unname(c(median_vals["SP_ols"], median_vals["SP_rrr"], median_vals["SP_srrr"])),
    Signal_PE_IQR = unname(c(IQR["SP_ols"], IQR["SP_rrr"], IQR["SP_srrr"])),

    Estimation_Error_mean = unname(c(avg["Est_ols"], avg["Est_rrr"], avg["Est_srrr"])),
    Estimation_Error_median = unname(c(median_vals["Est_ols"], median_vals["Est_rrr"], median_vals["Est_srrr"])),
    Estimation_Error_IQR = unname(c(IQR["Est_ols"], IQR["Est_rrr"], IQR["Est_srrr"])),

    Number_of_non_zero_rows = unname(c(NA, NA, avg["nonzero_srrr"])),
    True_positive_rate = unname(c(NA, NA, avg["True_pr"])),
    False_positive_rate = unname(c(NA, NA, avg["False_pr"])),
    
    row.names = NULL
  )
  results_table
  
}

####################  SCENARIOS #######################

#the following code defines each simulation scenario
set.seed(123)
par(mfrow = c(1,1))

## baseline scenario p < n ##
scenario1 <- run_scenario(
  scenario_name = "Scenario 1: p<n, p=65, p0=30",
  nrep=100,
  n= 100,
  p= 65,
  p0 = 30, #creates moderate sparsity in the true coefficient matrix
  q = 50,
  q0 = 50,
  nrank = 3,
  s2n = 1,
  rho_X = 0.5,
  rho_E = 0
)

## no sparsity scenario ##
scenario2 <- run_scenario(
  scenario_name = "Scenario 2: p<n, p=65, p0=65",
  nrep=100,
  n= 100,
  p= 65,
  p0 = 65, # the true coefficient matrix is not sparse in this scenario
  q = 50,
  q0 = 50,
  nrank = 3,
  s2n = 1,
  rho_X = 0.5,
  rho_E = 0
)

## high dimension ##
scenario3 <- run_scenario(
  scenario_name = "Scenario 3: p>n,p=200, p0=75",
  nrep=100,
  n= 100,
  p= 200, #p > n tests the methods in high-dimensions
  p0 = 75, #creates a moderately sparse scenario
  q = 50,
  q0 = 50,
  nrank = 3,
  s2n = 1,
  rho_X = 0.5,
  rho_E = 0
)

# ## strong sparsity ##
scenario4 <- run_scenario(
  scenario_name = "Scenario 4: p>n, p=200, p0=10",
  nrep=100,
  n= 100,
  p= 200,
  p0 = 10, #creates a highly sparse setting
  q = 50,
  q0 = 50,
  nrank = 3,
  s2n = 1,
  rho_X = 0.5,
  rho_E = 0
)


##  Added Scenarios: uncorrelated predictors  ##


scenario5 <- run_scenario(
  scenario_name = "Scenario 5: uncorrelated predictors, p=200, p0=75",
  nrep=100,
  n= 100,
  p= 200,
  p0 = 75,
  q = 50,
  q0 = 50,
  nrank = 3,
  s2n = 1,
  rho_X = 0,
  rho_E = 0
)

scenario6 <- run_scenario(
  scenario_name = "Scenario 6: uncorrelated predictors, p=200, p0=10",
  nrep=100,
  n= 100,
  p= 200,
  p0 = 10,
  q = 50,
  q0 = 50,
  nrank = 3,
  s2n = 1,
  rho_X = 0,
  rho_E = 0
)


################## RESULTS #####################


all_results <- rbind(scenario1, scenario2, scenario3, scenario4, scenario5,
                     scenario6)
View(all_results)

