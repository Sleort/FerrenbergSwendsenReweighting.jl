# FerrenbergSwendsenReweighting.jl

[![Build Status](https://travis-ci.org/Sleort/FerrenbergSwendsenReweighting.jl.svg?branch=master)](https://travis-ci.org/Sleort/FerrenbergSwendsenReweighting.jl)

[![Coverage Status](https://coveralls.io/repos/Sleort/FerrenbergSwendsenReweighting.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Sleort/FerrenbergSwendsenReweighting.jl?branch=master)

[![codecov.io](http://codecov.io/github/Sleort/FerrenbergSwendsenReweighting.jl/coverage.svg?branch=master)](http://codecov.io/github/Sleort/FerrenbergSwendsenReweighting.jl?branch=master)

A Julia package for single and multiple histogram[^1] reweighting à la Ferrenberg and Swendsen.

Calculates the reweighting weights - the "reweights".

[^1]: The nomenclature is kind misleading. No histogram is actually constructed in the reweighting procedure.

## Installation

## Usage/examples





####Input:

* Single histogram: 
  * *The logarithm* of the probability weight of a configuration $x$ at parameter $\lambda$ (NB: can be a vector!), $s: (\lambda, x) \to \mathbb{R}$. 
  * A series of data $\{x\}$ (in the form of a vector) sampled from the distribution given by $s$ at parameter $\lambda_0$. 
  * The parameter $\lambda$ which we want to "reweight to".
* Multiple histograms: 
  * Like above, but now the input is a set of series of data $\{\{x\}_i\}$ obtained from samplings at $\{\lambda_i\}$. 
  * We also need to provide information on the integrated autocorrelation time for an observable $O$ for each series $O$, $\tau_\text{int}(O, \{x\}_i)$. Alternatively, if no input is given, the autocorrelation time of $O = s$ is used.

####Output:

A vector of weights associated with the datapoints. (Or vector of vector of weights, in the case of multiple histograms)

## Equations

### Single histogram reweighting

Given a set of samples $\{x\}_{\lambda_0}$ obtained from some distribution parametrized by $\lambda_0$ (e.g. by a Monte Carlo simulation), we want the (best) set of associated weights $\{w\}_\lambda$ we have to use to get the correct expectation value of observable $O$ at (arbitrary) $\lambda$. In other words, we want to find $w(\lambda, x)$ such that
$$
\langle O \rangle_\lambda = \frac{\sum_x w(\lambda, x) O(x)}{\sum_x w(\lambda, x)}
$$
Note that by construction $w(\lambda_0, x) = \text{constant}$.

Let the probability weight of a state $x$ at parameter $\lambda$ be given by 
$$
p(\lambda, x) \propto \exp(s(\lambda, x))
$$
then
$$
w(\lambda, x) = \exp[s(\lambda, x) - s(\lambda_0,x)]
$$

### Multiple histogram reweighting

Now we have a set of set of samples $\{\{x\}_{\lambda_i}\}$, each obtained at some $\lambda_i$. According to Ferrenberg and Swendsen, the optimal weights can be found by solving the nonlinear set of equations given by
$$
1 = \sum_i \sum_{x_i } \frac{{g_i}^{-1}\exp[s(\lambda, x_i) - \log Z_\lambda ]}{\sum_j N_j {g_j}^{-1} \exp[s(\lambda_j, x_i) - \log Z_{\lambda_j}]} \quad \forall \lambda \in \{\lambda_i\}
$$
for all $\{\log Z_{\lambda_i}\}$ (up to an overall constant factor, which can be fixed arbitrarily). Here $g_i = 2\tau_{\text{int}, i} \ge 1$, with $\tau_{\text{int},i}$ being the integrated autocorrelation time of the simulation series $i$. (Uncorrelated data corresponds to $\tau_\text{int} = 1/2$, i.e. $g=1$.) $N_i$ is the number of samples in series $i$.

Once $\{\log Z_{\lambda_i}\}$ is found, the weights at an arbitrary $\lambda$ can be found by (ignoring irrelevant constant scaling)
$$
w(\lambda, x_i) = \frac{{g_i}^{-1}\exp[s(\lambda, x_i) ]}{\sum_j N_j {g_j}^{-1} \exp[s(\lambda_j, x_i) - \log Z_{\lambda_j}]}
$$

#### Technicalities

* In solving the nonlinear equation, the scale is fixed such that $Z_{\lambda_1} = 1$. 

* To avoid numerical overflow, the equation solved is actually 

* $$
  1 = \sum_i \sum_{x_i } \frac{{g_i}^{-1}}{\sum_j N_j {g_j}^{-1} \exp[s(\lambda_j, x_i) - s(\lambda, x_i) - \log Z_{\lambda_j}  + \log Z_\lambda]} \quad \forall \lambda \in \{\lambda_i\}
  $$

  This guarantees that at least the $i = j$ terms give a reasonably large, non-zero contribution to the sum in the denominator. On the other hand, an numerical "infinite" exponential in the denominator simply leads to a zero contribution to the overall sum, thus preventing any overflow instability.

* The output weights are scaled such that the largest weight is of unit size.




