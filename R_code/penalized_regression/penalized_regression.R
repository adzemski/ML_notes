library(ggplot2)
library(glmnet)
library(data.table)
library(xtable)
library(MASS)

# function that generates data for the "dense" design
generate_data_dense <- function() {
  # round up to next even number
  p_prime <- p + p %% 2
  Sigma <- kronecker(diag(p_prime/2), matrix(c(1, corr_par, corr_par, 1),
                                             nrow = 2, ncol = 2))
  Sigma <- Sigma[seq_len(p), seq_len(p)]
  
  X <- mvrnorm(n, mu = rep_len(0, p), Sigma = Sigma)
  
  y <- X %*% coef_true$b + rnorm(n, sd = 3)
  
  setNames(data.table(cbind(y, X)), c("y", sprintf("x%d", seq_len(p))))
}

# function that generates data for the "sparse" desgin w/ high correlation of 
# the first 2 regressors
generate_data_sparse_corr <- function() {
  # round up to next even number
  Sigma <- diag(p)
  
  Sigma[seq_len(2), seq_len(2)] <- 
    matrix(c(1, 0.95, 0.95, 1), nrow = 2, ncol = 2)
  
  X <- mvrnorm(n, mu = rep_len(0, p), Sigma = Sigma)
  
  y <- X %*% coef_true$b + rnorm(n, sd = 2)
  
  setNames(data.table(cbind(y, X)), c("y", sprintf("x%d", seq_len(p))))
}

# function that generates data for the "sparse" desgin w/o correlation of 
# the first 2 regressors
generate_data_sparse_uncorr <- function() {
  # round up to next even number
  Sigma <- diag(p)
  
  X <- mvrnorm(n, mu = rep_len(0, p), Sigma = Sigma)
  
  y <- X %*% coef_true$b + rnorm(n, sd = 2)
  
  setNames(data.table(cbind(y, X)), c("y", sprintf("x%d", seq_len(p))))
}

# draw random data and fit regression models once
simulate_once <- function(s, generate_data, 
                          lambda_ridge = "lambda.min", lambda_lasso = "lambda.min", 
                          always_include_first = FALSE) {  
  
  training_data <- generate_data()
  verification_data <- generate_data()
  
  f_true <- as.vector(as.matrix(verification_data[, 2:(p+1)]) %*% coef_true$b)
  
  # data contains regressors and one outcome (-> subtract 1)
  p <- ncol(training_data) - 1 
  
  if (always_include_first) {
    penalty_factor <- c(0, rep(1, p - 1))
  } else {
    penalty_factor <- rep(1, p)
  }
  
  ols_fit <- lm(training_data$y ~ as.matrix(training_data[, 2:(p+1)]))
  coef_ols <- coef(ols_fit)[2:(p+1)]
  predict_ols <- predict(ols_fit, x = as.matrix(verification_data[, 2:(p+1)]))
  
  ridge_fit <- 
    cv.glmnet(as.matrix(training_data[, 2:(p+1)]), training_data$y, 
              family = "gaussian", alpha = 0, standardize = FALSE, 
              penalty.factor = penalty_factor)
  coef_ridge <- as.matrix(coef(ridge_fit, s = lambda_ridge))[2:(p+1)]
  predict_ridge <- 
    as.vector(predict(ridge_fit, newx = as.matrix(verification_data[, 2:(p+1)]), 
                      s = lambda_ridge))
  
  lasso_fit <- 
    cv.glmnet(as.matrix(training_data[, 2:(p+1)]), training_data$y, 
              family = "gaussian", alpha = 1, standardize = FALSE, 
              penalty.factor = penalty_factor)
  coef_lasso <- as.matrix(coef(lasso_fit, s = lambda_lasso))[2:(p+1)]
  predict_lasso <- 
    as.vector(predict(lasso_fit, newx = as.matrix(verification_data[, 2:(p+1)]), 
                      s = lambda_lasso))
  
  selected <- abs(coef_lasso) > 0
  post_lasso_fit <- lm(training_data$y ~ as.matrix(training_data)[, c(FALSE, selected)])
  coef_post_lasso <- numeric(p)
  coef_post_lasso[selected] <- coef(post_lasso_fit)[-1]
  predict_post_lasso <- predict(post_lasso_fit, x = as.matrix(verification_data[, 2:(p+1)]))

  list(coef_ols=coef_ols, predict_ols = predict_ols,
       coef_ridge=coef_ridge, predict_ridge = predict_ridge, 
       coef_lasso=coef_lasso, predict_lasso = predict_lasso,
       coef_post_lasso = coef_post_lasso, predict_post_lasso = predict_post_lasso, 
       f_true = f_true, 
       lambda.min_ridge = ridge_fit$lambda.min, 
       lambda.min_lasso = lasso_fit$lambda.min)
}

# summarize coefficients estimates obtained in MC simulations
summarize_sim_coef <- function(method, sim_results) {
  coef <- 
    do.call(rbind, lapply(sim_results, function(x) x[[paste0("coef_", method)]]) )
  
  data.frame(method = rep_len(method, p),
             j = seq_len(p), b = as.vector(colMeans(coef)), 
             included = as.vector(colMeans(abs(coef) > 0)), 
             stringsAsFactors = FALSE)
}

# estimate empirical MSE from MC simulations
summarize_sim_mse <- function(method, sim_results) {

  prediction <- 
    do.call(c, lapply(sim_results, function(x) x[[paste0("predict_", method)]]) )
  
  f_true <- 
    do.call(c, lapply(sim_results, function(x) x$f_true) )
  
  data.frame(method = method, mse = mean((prediction - f_true)^2))
}

# draw data once and fit regression models and output results
simulate_design_once <- function(generate_data, file, methods = c("ols", "ridge", "lasso"), 
                                 always_include_first = FALSE) {
  
  coef_as_df <- function(method) {
    data.frame(method = rep_len(method, p), j = seq_len(p), 
               b = simulation_result[[paste0("coef_", method)]])
  }
  
  simulation_result <- simulate_once(1, generate_data = generate_data, 
                                     always_include_first = always_include_first)
  
  fit_once_coef <- 
    rbind(do.call(rbind, lapply(methods, coef_as_df)),
          coef_true[, c("method", "j", "b")])
  
  plot_coef <-   ggplot2::ggplot(fit_once_coef, aes(x = j, y = b, color = as.factor(method))) + 
    geom_point() + theme_bw() + labs(color = "Method") + xlab("regressor index") + ylab("estimate")
  
  pdf(file, width = 4.2, height = 2.5)
    print(plot_coef)
  dev.off()
  
  plot_coef
}

# run an MC simulation, generate summaries and save them
simulate_design_MC <- function(generate_data, file_prefix, num_sim = 50, 
                               methods = c("ols", "ridge", "lasso"), 
                               lambda.min_ridge = "lambda.min", lambda.min_lasso = "lambda.min", 
                               always_include_first = FALSE) {
  
  sim_results <- lapply(seq_len(num_sim), simulate_once, lambda_ridge = lambda.min_ridge, 
                        lambda_lasso = lambda.min_lasso,
                        generate_data = generate_data, 
                        always_include_first = always_include_first)
  
  lambda.min_ridge <- mean(do.call(c, lapply(sim_results, function(x) x$lambda.min_ridge) )) 
  lambda.min_lasso <- mean(do.call(c, lapply(sim_results, function(x) x$lambda.min_lasso) )) 

  summary_sim_coef <- 
    rbind(do.call(rbind, 
                  lapply(methods, summarize_sim_coef, sim_results = sim_results)), 
          cbind(coef_true[c("method", "j", "b")], data.frame(included = rep.int(NA, nrow(coef_true))))
          )
  
  summary_sim_coef <- data.table(summary_sim_coef)
  
  plot_included <- 
    ggplot(summary_sim_coef[method %in% c("ridge", "lasso"), ], 
           aes(x = j, y = included, color = as.factor(method))) + 
    geom_point() + theme_bw() + labs(color = "Method") + xlab("regressor index") + 
    ylab("probability of selection")
  
  pdf(paste0(file_prefix, "_plot_included.pdf"), width = 4.2, height = 2.5)
    print(plot_included)
  dev.off()
  
  sum_sim_mse <- 
    do.call(rbind, lapply(methods, summarize_sim_mse, sim_results = sim_results))
  sum_sim_mse
  print(xtable(sum_sim_mse), booktabs = TRUE, file = paste0(file_prefix, "_mse_table.tex"))
  
  plot_coef <- ggplot2::ggplot(summary_sim_coef, aes(x = j, y = b, color = as.factor(method))) + 
    geom_point() + theme_bw() + labs(color = "Method") + xlab("regressor index") + 
    ylab("mean estimate")
  
  pdf(paste0(file_prefix, "_plot_mean.pdf"), width = 4.2, height = 2.5)
    print(plot_coef)
  dev.off()
  
  list(lambda.min_ridge = lambda.min_ridge, lambda.min_lasso = lambda.min_lasso, 
       sum_sim_mse = sum_sim_mse, 
       plot_included = plot_included, 
       plot_coef = plot_coef, summary_sim_coef = summary_sim_coef)
}

# set working directory
wd <- "/Users/xdzean/Dropbox/teaching/machine_learning/ML_notes/R_code/penalized_regression"
# the folder "output_for_slides" must exist in wd
setwd(file.path(wd, "output_for_slides"))

### Simulate a model that is NOT sparse

# number of regressors/features 
p <- 50
# number of observations 
n <- 100
# parameterization of some correlations
corr_par <- 0.8

# define coefficients
b <- seq(0.1, 1, length.out = p)
b <- abs(b)^(4)
coef_true <- data.frame(method = rep_len("truth", p), 
                        j = seq_len(p), b = b, group = ((seq_len(p)-1) %/% 2) %% 2)

# plot coefficients
pdf("dense_plot_coefficients.pdf", width = 4.2, height = 2.5)
print(ggplot2::ggplot(coef_true, aes(x = j, y = b, color = as.factor(group))) + 
        geom_point() + theme_bw() + guides(color = FALSE) + ylab("True value of coefficient"))
dev.off()

# fit 1 model
set.seed(20190516)
plot_once <- simulate_design_once(generate_data_dense, "dense_plot_coef_once.pdf")

# run MC simulation
MC_sim_dense <- simulate_design_MC(generate_data_dense, "dense", num_sim = 200, lambda.min_ridge = 3.2,
                                   lambda.min_lasso = 0.4)


### Simulate a model that is IS sparse

# define coefficients
b <- rep(0, length.out = p)
b[1] <- 0.5
b[2:10] <- 1
coef_true <- data.frame(method = rep_len("truth", p),
                        j = seq_len(p), b = b)

# plot coefficients
pdf("sparse_plot_coefficients.pdf", width = 4.2, height = 2.5)
print(ggplot2::ggplot(coef_true, aes(x = j, y = b)) + 
        geom_point() + theme_bw() + guides(color = FALSE) + ylab("True value of coefficient"))
dev.off()

# fit 1 model
set.seed(20190516)
plot_once <-
  simulate_design_once(generate_data_sparse_uncorr, "sparse_uncorr_plot_coef_once.pdf")

## uncorrelated case

# run MC simulation 
MC_sim_sparse_uncorr <-
  simulate_design_MC(generate_data_sparse_uncorr, "sparse_uncorr", num_sim = 200, 
                     lambda.min_ridge = 0.9, lambda.min_lasso = 0.2)

## uncorrelated case, fix first regressor, demonstration of post-lasso

# run MC simulation
MC_sim_sparse_uncorr_post <-
  simulate_design_MC(generate_data_sparse_uncorr, "sparse_uncorr_post",
                     num_sim = 200, lambda.min_ridge = "lambda.min",
                     lambda.min_lasso = "lambda.min", methods = c("ols", "lasso", "post_lasso"), 
                     always_include_first = TRUE)

## correlated case, fix first regressor, demonstration of post-lasso

# now run design with strong correlation in 1st 2 regressors
# run MC simulation
MC_sim_sparse_corr_post <-
  simulate_design_MC(generate_data_sparse_corr, "sparse_corr_post",
                     num_sim = 200, lambda.min_ridge = "lambda.min",
                     lambda.min_lasso = "lambda.min", methods = c("ols", "lasso", "post_lasso"), 
                     always_include_first = TRUE)

## cross-validation example 

set.seed(20180524)
training_data <- generate_data_sparse_corr()

lasso_cv <- 
  cv.glmnet(as.matrix(training_data[, 2:(p+1)]), training_data$y, 
            family = "gaussian", alpha = 1, standardize = FALSE, nfolds = 10, 
            lambda = exp(seq(log(0.002), log(1.6), length.out = 30)))

df_lasso_cv <- data.frame(lambda=lasso_cv$lambda, cvm=lasso_cv$cvm, 
                          upper = lasso_cv$cvm + 2*lasso_cv$cvsd,
                          lower = lasso_cv$cvm - 2*lasso_cv$cvsd)

plot_lasso_cv <- ggplot(df_lasso_cv, aes(x=log(lambda), y=cvm)) +
  geom_point(color = "red") + theme_bw() + geom_errorbar(aes(ymin = lower, ymax = upper), color="grey") 

plot_lasso_cv

pdf("cv_example_lambda_path.pdf", width = 4.2, height = 2.5)
  plot_lasso_cv
dev.off()

plot_lasso_cv <- plot_lasso_cv + 
  geom_vline(xintercept = log(lasso_cv$lambda.min), linetype = "dashed") +
  geom_vline(xintercept = log(lasso_cv$lambda.1se), linetype = "dashed") 

pdf("cv_example_lambda_choice.pdf", width = 4.2, height = 2.5)
  plot_lasso_cv
dev.off()
