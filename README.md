<table>
<tr>
    <td>
        <img src="https://github.com/federicopirola/LUME-DBN/blob/main/images/LUME-DBN.svg" 
             alt="LUME-DBN logo" 
             width="100" 
             height="100">
    </td>
    <td>
        <h1>LUME-DBN: Latent Uncertainty Modeling via MCMC Estimation in Dynamic Bayesian Networks</h1>
    </td>
</tr>
</table>
Latent Uncertainty Modeling via MCMC Estimation in DBNs (LUME-DBN) is a novel Gibbs sampling-based method for learning Dynamic Bayesian Networks from incomplete data. which iteratively sample the parameters and the missing values from their Full Conditional Distribution and update the DBN structure via Metropolis Hasting moves.

LUME-DBN jointly learns network structure, regression parameters, and missing values within a unified Bayesian framework. The workflow can be described in three main phases, repeated iteratively within an MCMC loop:

1. **Parameter learning (Collapsed Gibbs sampling):**  
   For each variable, regression coefficients, noise variances, and uncertainty parameters are updated from their posterior distributions, conditional on the current parent set and data.
<img src="https://github.com/federicopirola/LUME-DBN/blob/main/images/parameters%20update.png" alt="Parameters Update" width="1000"/>


2. **Structure learning (Metropolis–Hastings search):**  
   Parent sets are modified by proposing local moves (addition, deletion, or exchange of an edge). Proposals are accepted or rejected based on marginal likelihood ratios and structural priors, encouraging parsimonious networks.
<img src="https://github.com/federicopirola/LUME-DBN/blob/main/images/structure%20update.png" alt="Structure Update" width="1000"/>


3. **Missing data imputation (Gibbs step):**  
   At regular intervals, missing values are resampled from their full conditional distributions, which are Gaussian under the model assumptions. Each update accounts simultaneously for a variable’s own regression and its influence on child nodes, ensuring consistency across the network.
<img src="https://github.com/federicopirola/LUME-DBN/blob/main/images/missing%20values%20update.png" alt="Missing Update" width="1000"/>

<h1>Code Overview</h1>

The `code` folder contains all the functions and scripts used for Bayesian Linear Regression (BLR), Dynamic Bayesian Networks (DBNs), and Non-Homogeneous DBNs (NH-DBNs) analysis, learning, simulation, and visualization. It is organized into the following subfolders:

## distributions
Contains functions for:
- Sampling and computing densities of classic distributions, including Gaussian, Inverse Gamma, Poisson, Negative Binomial, and Geometric.  
- Sampling parameters and missing values from the Full Conditional Distribution.  
- Computing conditional likelihoods, prior probabilities, and marginal likelihoods.

## learning
Contains functions and scripts for:
- Metropolis-Hastings structural moves.  
- Collapsed Gibbs sampling for parameter moves.  
- Gibbs sampling for missing values imputation.  
- Learning single Bayesian Linear Regression models, Dynamic Bayesian Networks (DBNs) with complete data, and the LUME-DBN algorithm for learning DBNs from incomplete data.  
- Extensions for learning non-homogeneous DBNs in the case of complete data.

## random_generation
Contains functions for:
- Randomly generating Bayesian Linear Regression parameters.  
- Randomly generating DBN and NH-DBN structures and parameters.  
- Generating datasets from BLRs and DBNs.  
- Creating random incomplete datasets from complete data.

## plots
Contains functions for:
- Visualizing BLR, DBN, and NH-DBN models.  
- Plotting posterior sample distributions.  
- Generating missing value credible intervals.  
- Performing convergence diagnostics and visualizations.

## utils
Contains additional utility functions frequently used across simulations, learning, and plotting scripts.


📧 Federico Pirola: [f.pirola17@campus.unimib.it](mailto:f.pirola17@campus.unimib.it)
