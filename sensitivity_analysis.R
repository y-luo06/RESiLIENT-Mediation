## Time-varying mediation analysis
## Parametric g-formula
## Sensitivity analysis: Complete-case analysis

#------------------------------
# Settings
#------------------------------
set.seed(0)

d <- isi_phq50.cca
a <- "C5"        # intervention
a_star <- "C10"  # reference
K <- 50          # number of permutation repeats
n <- nrow(d)

# Rename for convenience
# Timeline: C -> A -> M1 -> L -> M2 -> Y
d <- d %>% rename(A = assignment,
                  L = phq9no3w5,
                  Y = phq9no3w50,
                  M1 = isipoint_w3,
                  M2 = isipoint_w6)

# Baseline confounders
C <- c("age", "sex", "education", "working", "cohab", "marital",
       "mentaltx", "ace4orgreater", "cage2orgreater",
       "big5neuroticism", "big5extraversion", "big5openness",
       "big5agreeableness", "big5conscientiousness",
       "wsasbase", "isibase", "phq9base", "gad7base", "swemwbsbase",
       "asc_ssbase", "asc_dbase", "asc_mbase")

#------------------------------
# Function to simulate Q estimates: Qa_a, Qa_aStar, QaStar_a, QaStar_aStar
#------------------------------
simulate_Q <- function(data) {
  
  # 1. Fit parametric models: C -> A -> M1 -> L -> M2 -> Y
  fit_M1 <- lm(as.formula(paste("M1 ~ A +", paste(C, collapse = " + "))), data = data)
  fit_L  <- lm(as.formula(paste("L  ~ A + M1 +", paste(C, collapse = " + "))), data = data)
  fit_M2 <- lm(as.formula(paste("M2 ~ A + M1 + L +", paste(C, collapse = " + "))), data = data)
  fit_Y  <- lm(as.formula(paste("Y  ~ A + M1 + L + M2 +", paste(C, collapse = " + "))), data = data)
  
  ## variance of each model will be used to add noise in the simulation later
  sigma_M1 <- sigma(fit_M1)
  sigma_L  <- sigma(fit_L)
  sigma_M2 <- sigma(fit_M2)
  sigma_Y  <- sigma(fit_Y)
  
  # 2. Store Q estimates
  ## Qa_a = Yama (Y in the BI where M at BI level), similar logic to others
  Q_store <- matrix(NA, nrow = K, ncol = 4)
  colnames(Q_store) <- c("Qa_a", "Qa_aStar", "QaStar_a", "QaStar_aStar")
  
  # 3. Main loop over K permutations (Permutation is used for random draw of M_a and M_aStar)
  for (k in 1:K) {
    
    ## 3a: simulate mediators under A = a 
    newdat_a <- data
    newdat_a$A <- a
    
    ### C are as observed, then M1 based on A and C; L based A, C, M1; M2 based on A, C, M1, and L: 
    ### simulate these for each i based on the fitted models, and + noise
    sim_M1_a <- rnorm(n, mean = predict(fit_M1, newdata = newdat_a), sd = sigma_M1)
    sim_L_a  <- rnorm(n, mean = predict(fit_L, newdata = transform(newdat_a, M1 = sim_M1_a)), sd = sigma_L)
    sim_M2_a <- rnorm(n, mean = predict(fit_M2, newdata = transform(newdat_a, M1 = sim_M1_a, L = sim_L_a)), sd = sigma_M2)
    
    ### randomly permute n values of the joint mediators (permutation): interventional or random analogue E
    ### instead of the most likely M value, we randomly assign the M based on the distribution (MARGINAL) to each subject
    G_a_M1 <- sample(sim_M1_a)
    G_a_M2 <- sample(sim_M2_a)
    
    ## 3b: simulate mediators under A = a*
    newdat_as <- data
    newdat_as$A <- a_star
    
    sim_M1_as <- rnorm(n, mean = predict(fit_M1, newdata = newdat_as), sd = sigma_M1)
    sim_L_as  <- rnorm(n, mean = predict(fit_L, newdata = transform(newdat_as, M1 = sim_M1_as)), sd = sigma_L)
    sim_M2_as <- rnorm(n, mean = predict(fit_M2, newdata = transform(newdat_as, M1 = sim_M1_as, L = sim_L_as)), sd = sigma_M2)
    
    ### random sample M from its distribution to replace the "naturally observed Ms" --> INTERVENTIONAL
    G_as_M1 <- sample(sim_M1_as)
    G_as_M2 <- sample(sim_M2_as)
    
    ## 3c: estimate Qa,a, Qa,a*, Qa*,a, Qa*,a* 
    ### Qa,a
    #### Permutation: repeat M1 and M2 with the permuted G, then predict Y
    newdata_Y <- data.frame(A = a, M1 = G_a_M1, L = sim_L_a, M2 = G_a_M2)
    newdata_Y[, C] <- data[, C]
    Q_store[k, "Qa_a"] <- mean(predict(fit_Y, newdata = newdata_Y))
    
    ### Qa,a*
    newdata_Y <- data.frame(A = a, M1 = G_as_M1, L = sim_L_a, M2 = G_as_M2)
    newdata_Y[, C] <- data[, C]
    Q_store[k, "Qa_aStar"] <- mean(predict(fit_Y, newdata = newdata_Y))
    
    ### Qa*,a
    newdata_Y <- data.frame(A = a_star, M1 = G_a_M1, L = sim_L_as, M2 = G_a_M2)
    newdata_Y[, C] <- data[, C]
    Q_store[k, "QaStar_a"] <- mean(predict(fit_Y, newdata = newdata_Y))
    
    ### Qa*,a*
    newdata_Y <- data.frame(A = a_star, M1 = G_as_M1, L = sim_L_as, M2 = G_as_M2)
    newdata_Y[, C] <- data[, C]
    Q_store[k, "QaStar_aStar"] <- mean(predict(fit_Y, newdata = newdata_Y))
  }
  
  # 4. Aggregate over K permutations
  colMeans(Q_store)
}

#------------------------------
# Estimate effects on original data using the above function
#------------------------------
Q_hat <- simulate_Q(d)

effects <- list(
  rTE  = Q_hat["Qa_a"] - Q_hat["QaStar_aStar"],
  rNDE = Q_hat["Qa_aStar"] - Q_hat["QaStar_aStar"],
  rNIE = Q_hat["Qa_a"] - Q_hat["Qa_aStar"],
  rPM  = (Q_hat["Qa_a"] - Q_hat["Qa_aStar"]) / (Q_hat["Qa_a"] - Q_hat["QaStar_aStar"])
)

print(Q_hat)
print(effects)

#------------------------------
# Bootstrap 95% CI
#------------------------------
library(parallel)

## Parallelized approach for boot 95%CI: fast
B <- 1000
set.seed(1)

## Function for a single bootstrap iteration
do_bootstrap <- function(b) {
  d_boot <- d[sample(1:n, size = n, replace = TRUE), ]
  Q_hat_boot <- simulate_Q(d_boot)
  
  rTE   <- Q_hat_boot["Qa_a"] - Q_hat_boot["QaStar_aStar"]
  rNDE  <- Q_hat_boot["Qa_aStar"] - Q_hat_boot["QaStar_aStar"]
  rNIE  <- Q_hat_boot["Qa_a"] - Q_hat_boot["Qa_aStar"]
  rPM   <- rNIE / rTE
  
  c(rTE, rNDE, rNIE, rPM)
}

## Parallel execution
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("d", "n", "simulate_Q", "a", "a_star", "K", "C"))

boot_effects <- t(parSapply(cl, 1:B, do_bootstrap))
stopCluster(cl)

colnames(boot_effects) <- c("rTE", "rNDE", "rNIE", "rPM")

CI_lower <- apply(boot_effects, 2, quantile, probs = 0.025)
CI_upper <- apply(boot_effects, 2, quantile, probs = 0.975)
results_CI <- data.frame(
  Effect = c("rTE", "rNDE", "rNIE", "rPM"),
  Estimate = unlist(effects),
  CI_lower = CI_lower,
  CI_upper = CI_upper
)
print(results_CI, digits = 4)
