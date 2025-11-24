---
title: "GcompBART"
output: github_document
---

# Introduction
GcompBART is a framework for performing G-computation using Bayesian Additive Regression Trees. It leverages the modeling capabilities of SoftBart (<https://cran.r-project.org/package=SoftBart>) internally to fit flexible observed data models, with the optional use of sparsity-inducing priors for longitudinal data as described in "Long-term cognitive effects of an incremental blood pressure intervention in a mortal cohort". Josefsson, Karalija and Daniels, https://doi.org/10.48550/arXiv.2311.00357. The package follows the same interface style as the SoftBart package for ease of use.

# Installation
From CRAN (once published):
```r
install.packages("GcompBART")
```

## Using remotes
```r
remotes::install_github("m4ryjo/GcompBART")
```

## Or devtools
```r
devtools::install_github("m4ryjo/GcompBART")
```

```r
library(GcompBART)
```


# Example: Simulate Data and Perform G-Computation
```r
# Simulation and model fitting example
n_burn <- 10
n_thin <- 1
n_save <- 10
n_J <- 100

set.seed(1234)

# Simulate dataset
n <- 100 
times_sim <- rep(1:4, each = 5) # grouping variable
lambda <- c(0.25,0.5,0.75,1)
a1 <- 0.4
a2 <- 0.2
a3 <- 0.1
sig <- 10

Sig <- matrix(c(c(1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0,0,0,0),
              c(0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0,0,0),
              c(0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0,0),
              c(0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0),
              c(0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3),
              c(a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0),
              c(0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0),
              c(0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0),
              c(0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0),
              c(0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2),
              c(a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0),
              c(0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0),
              c(0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0),
              c(0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0),
              c(0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1),
              c(a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0),
              c(0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0),
              c(0,0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0),
              c(0,0,0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0),
              c(0,0,0,0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1)), 20, 20, byrow = TRUE)

vartype_bl <- c(rep("X0", 5), rep("X", 15), "Y")
tgroup <- c(rep(1:4,each=5),4)

training_data <- sim_lfried(n, 20, Sig, sig, lambda)

# Fit the model and perform g-computation
BM <- BMfits(training_data[,1:21],
             var.type = vartype_bl,
             opts = Opts(num_burn = n_burn, 
                         num_thin = n_thin, 
                         num_save = n_save,
                         update_s = TRUE,
                         update_alpha = TRUE,
                         update_tvp = FALSE,
                         update_alpha_vec = FALSE,
                         update_eta = FALSE,
                         update_phi = FALSE,
                         update_tau = TRUE,
                         update_sigma_mu = TRUE),
             tgroup = tgroup)

out_gcomp <- gcompbart(training_data[,1:21],
                       BModels = BM,
                       var.type = vartype_bl,
                       J = n_J,
                       opts = Opts(num_burn = n_burn, 
                                   num_thin = n_thin, 
                                   num_save = n_save,
                                   update_s = TRUE,
                                   update_alpha = TRUE,
                                   update_tvp = FALSE,
                                   update_alpha_vec = FALSE,
                                   update_eta = FALSE,
                                   update_phi = FALSE,
                                   update_tau = TRUE,
                                   update_sigma_mu = TRUE),
                       tgroup = tgroup)

# Results
out_gcomp$summary_out
```

## Extension: Truncation by Death
Here is an example of how you might implement truncation by death in your simulated dataset:

```r

# Simulation and model fitting example
n_burn <- 10
n_thin <- 1
n_save <- 10
n_J <- 100

set.seed(1234)

# Simulate dataset
n <- 100 
times_sim <- rep(1:4, each = 5) # grouping variable
lambda <- c(0.25,0.5,0.75,1)
a1 <- 0.4
a2 <- 0.2
a3 <- 0.1
sig <- 10

Sig <- matrix(c(c(1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0,0,0,0),
              c(0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0,0,0),
              c(0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0,0),
              c(0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3,0),
              c(0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0,a3),
              c(a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0),
              c(0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0),
              c(0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0),
              c(0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0),
              c(0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2),
              c(a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0),
              c(0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0),
              c(0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0),
              c(0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0),
              c(0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1),
              c(a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0),
              c(0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0),
              c(0,0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0),
              c(0,0,0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0),
              c(0,0,0,0,a3,0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1)), 20, 20, byrow = TRUE)

training_data <- sim_lfried(n, 20, Sig, sig, lambda)

# Add survival status for three follow-up time points
surv1 <- rbinom(n, 1, 0.9)
surv2 <- rbinom(n, 1, 0.85)
surv3 <- rbinom(n, 1, 0.8)

# Insert survival indicators into dataset at specified positions
training_data <- cbind(training_data[,1:5], surv1, training_data[,6:10], surv2, training_data[,11:15], surv3, training_data[,16:21])

# Update variable types and grouping
vartype_bl <- c(rep("X0", 5), "S", rep("X", 5), "S", rep("X", 5), "S", rep("X", 5), "Y")
tgroup <- c(rep(1:4, each = 6), 4)  # Each block now has 6 variables including survival indicator

# Fit the model and perform g-computation
BM <- BMfits(training_data,
             var.type = vartype_bl,
             opts = Opts(num_burn = n_burn, num_thin = n_thin, num_save = n_save,
                         update_s = TRUE, update_alpha = TRUE, update_tvp = FALSE,
                         update_alpha_vec = FALSE, update_eta = FALSE, update_phi = FALSE,
                         update_tau = TRUE, update_sigma_mu = TRUE),
             tgroup = tgroup)

out_gcomp <- gcompbart(training_data,
                       BModels = BM,
                       var.type = vartype_bl,
                       J = n_J,
                       opts = Opts(num_burn = n_burn, num_thin = n_thin, num_save = n_save,
                                   update_s = TRUE, update_alpha = TRUE, update_tvp = FALSE,
                                   update_alpha_vec = FALSE, update_eta = FALSE, update_phi = FALSE,
                                   update_tau = TRUE, update_sigma_mu = TRUE),
                       tgroup = tgroup)

# Results
out_gcomp$summary_out
```


## Documentation
See ?Opts and ?Hypers for details.

## License
GPL-3 
# GcompBART
