### Package changes from previous ei.Circles version 0.0.1-7

Two new functions `simula_BPF` and `simula_BPF_with_deviations` are included in the package to simulate, departing from the basic underlying model in `BPF`, both marginal election results and unit vote transfer matrices, using (i) exclusively the underlying model in `BPF` and (ii) this model including ecological fallacy effects.

When covariates are used, the model transition probabilities estimated in each unit depend on the values of the covariates in the unit. A new output has been included in the `BPF` function to account for this. When covariates are not NULL, the list output of `BPF` contains an array with the estimated model transition probabilities corresponding to each unit. 

A new option,`"hyper"`, has been added for the `local` argument of the `BPF` function, which is call `"hyper"`. When `local = "hyper"`, transition matrices are estimated for each unit, and the default global solution is obtained by aggregating these unit-level estimates. In this case, the estimate for each unit assumes a multi-hypergeometric distribution for the table's inner values, given the observed row and column margins. The maximum likelihood estimate for each unit is then determined by randomly searching in the vicinity of the translated initial estimated transition matrix for that unit. 

### Package changes from previous ei.Circles version 0.0.1-6

The algorithm to estimate vote transfers at polling units when `local = "lik"` has been changed. The package NlcOptim (>= 0.6) is now required for this option of the `local` argument of `BPF`. 
