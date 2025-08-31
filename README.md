<img src="https://github.com/federicopirola/LUME-DBN/blob/main/LUME-DBN.svg" alt="LUME-DBN logo" width="227" height="227"/>
<h1>LUME-DBN: Latent Uncertainty Modeling via MCMC Estimation in Dynamic Bayesian Networks</h1>
Latent Uncertainty Modeling via MCMC Estimation in DBNs (LUME-DBN) is a novel Gibbs sampling-based method for learning Dynamic Bayesian Networks from incomplete data. which iteratively sample the parameters and the missing values from their Full Conditional Distribution and update the DBN structure via Metropolis Hasting moves.

LUME-DBN jointly learns network structure, regression parameters, and missing values within a unified Bayesian framework. The workflow can be described in three main phases, repeated iteratively within an MCMC loop:

1. **Parameter learning (Collapsed Gibbs sampling):**  
   For each variable, regression coefficients, noise variances, and uncertainty parameters are updated from their posterior distributions, conditional on the current parent set and data.
<img src="https://github.com/federicopirola/LUME-DBN/blob/main/parameters%20update.png" alt="Parameters Update" width="1000"/>
2. **Structure learning (Metropolis–Hastings search):**  
   Parent sets are modified by proposing local moves (addition, deletion, or exchange of an edge). Proposals are accepted or rejected based on marginal likelihood ratios and structural priors, encouraging parsimonious networks.
<img src="https://github.com/federicopirola/LUME-DBN/blob/main/structure%20update.png" alt="Structure Update" width="1000"/>
3. **Missing data imputation (Gibbs step):**  
   At regular intervals, missing values are resampled from their full conditional distributions, which are Gaussian under the model assumptions. Each update accounts simultaneously for a variable’s own regression and its influence on child nodes, ensuring consistency across the network.
<img src="https://github.com/federicopirola/LUME-DBN/blob/main/missing%20values%20update.png" alt="Missing Update" width="1000"/>

📧 Federico Pirola: [f.pirola17@campus.unimib.it](mailto:f.pirola17@campus.unimib.it)
