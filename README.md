# hiper_glm

Course Project for **140.778.01 — Statistical Computing, Algorithm, and Software Development**  

A flexible implementation of generalized linear model (GLM) estimation in R.  
Supports multiple optimization methods including QR decomposition, BFGS, and Newton-Raphson for logistic regression.

## Project Overview

This project implements a customizable GLM estimation function in R, supporting multiple optimization methods:
- QR decomposition (for linear regression)
- BFGS (quasi-Newton optimization)
- Newton-Raphson (for logistic regression)

The core function `hiper_glm()` automatically selects the appropriate method based on the model type and user input. Users can also specify convergence thresholds and iteration settings through the `option` argument.
