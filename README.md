# FerrenbergSwendsenReweighting.jl

[![Build Status](https://travis-ci.org/Sleort/FerrenbergSwendsenReweighting.jl.svg?branch=master)](https://travis-ci.org/Sleort/FerrenbergSwendsenReweighting.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/lx4l0r1eu9w79oe3?svg=true)](https://ci.appveyor.com/project/Sleort/ferrenbergswendsenreweighting-jl)
[![Coverage Status](https://coveralls.io/repos/github/Sleort/FerrenbergSwendsenReweighting.jl/badge.svg?branch=master)](https://coveralls.io/github/Sleort/FerrenbergSwendsenReweighting.jl?branch=master)

A Julia package for single and multiple histogram<sup>[1](#footnote1)</sup> reweighting à la Ferrenberg and Swendsen.

Calculates the reweighting weights - the "reweights".

<a name="footnote1">1</a>: The nomenclature is kind of misleading. No histogram is actually constructed in the reweighting procedure.


1. [Alan M. Ferrenberg and Robert H. Swendsen "New Monte Carlo technique for studying phase transitions", Physical Review Letters **61** pp. 2635-2638 (1988)](http://dx.doi.org/10.1103/PhysRevLett.61.2635)
2. [Alan M. Ferrenberg and Robert H. Swendsen "Optimized Monte Carlo data analysis", Physical Review Letters **63** pp. 1195-1198 (1989)](http://dx.doi.org/10.1103/PhysRevLett.63.1195)

## Installation
Sice this package is not (yet) registered, you must install it by cloning. To clone it, use
```julia
Pkg.clone("https://github.com/Sleort/FerrenbergSwendsenReweighting.jl")
```

## Usage
```julia
using FerrenbergSwendsenReweighting #Load package
```
Say you have a function
```julia
logprob(parameter, state)
```
mapping from a `parameter` (e.g. temperature) and a `state` (e.g. energy, or spin configuration, or...) to the *logarithm* of the probability weight of the `state` at that `parameter`. You also have a set (vector) of `sampled_states` sampled according to the distribution given by `exp(logprob(p0, state))` at parameter `p0`. Then we can find the relative `weights` of the `sampled_states` at some other parameter `p` by
```julia
rw = ReweightObj(logprob, p0, sampled_states) #Make reweighting object
weights = evaluate(rw, p1) #Find the weights at some parameter p1
evaluate!(rw, p2, weights) #Find the weights at some parameter p2, evaluated in-place (overwriting `weights`)
```


If no `logprob::Function` is provided, i.e.

```julia
rw = ReweightObj(p0, sampled_states) #Make reweighting object assuming the Boltzmann distribution
```
the Boltzmann weight is assumed, i.e. `parameter = 1/T`, `state = energy` and `logprob(parameter, state) = -parameter*state` will be used.



If `p0` is a subtype of `AbstractVector` *and* `sampled_states` a subtype of `AbstractVector{<:AbstractVector}` ("vector of vectors"), it is assumed that we're dealing with several simulation series `(p0[1],sampled_states[1]),(p0[2],sampled_states[2]),... ` and *multiple histogram reweighting* will be employed. If not, *single histogram reweighting* is used.



In case of multiple histogram reweighting, we may provide a set of integrated autocorrelation times for each simulation series. If that is not done, the autocorrelation of the `autocorrelation_observable(parameter,state)` will be used. The default here is `autocorrelation_observable = logprob`:

```julia
typeof(sampled_states) <: AbstractVector{<:AbstractVector} #true
length(p0) == length(sampled_states) == 2 #true

rw = ReweightObj(p0, sampled_states) #Basic Boltzmann distribution reweighting
# or
rw = ReweightObj(logprob, p0, sampled_states) #logprob distribution reweighting
# or
rw = ReweightObj(logprob, p0, sampled_states; essfs=ones(2)) #logprob distribution reweighting, all effective sampling size (ESS) factors set to 1
# or
myobs(parameter,state) = parameter*state^2 #A custom observable
rw = ReweightObj(logprob, p0, sampled_states; autocorrelation_observable = myobs) #logprob distribution reweighting, ESS factor according to the myobs observable

#The rest proceeds as before...
weights = evaluate(rw, p1) #Find the weights at some parameter p1
evaluate!(rw, p2, weights) #Find the weights at some parameter p2, evaluated in-place (overwriting `weights`)
```



## Examples

*Work in progress*




## Old Readme + background (needs to be fixed)

***For some strange reason Github doesn't support Latex style math... The stuff below will be fixed later/moved to documentation...***


**Input:**

* Single histogram:
  * *The logarithm* of the probability weight of a configuration $x$ at parameter $\lambda$ (NB: can be a vector!), $s: (\lambda, x) \to \mathbb{R}$.
  * A series of data $\{x\}$ (in the form of a vector) sampled from the distribution given by $s$ at parameter $\lambda_0$.
  * The parameter $\lambda$ which we want to "reweight to".
* Multiple histograms:
  * Like above, but now the input is a set of series of data $\{\{x\}_i\}$ obtained from samplings at $\{\lambda_i\}$.
  * We also need to provide information on the effective sample size factor (ESS) (proportional to the inverse integrated autocorrelation time) for an observable $O$ for each series $O$, $\text{ESS}(O, \{x\}_i)$. Alternatively, if no input is given, the ESS of $O = s$ is used.

**Output:**

A vector of weights associated with the data points. (Or vector of vector of weights, in the case of multiple histograms)

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

Now we have a set of set of samples $\{\{x\}_{\lambda_i}\}​$, each obtained at some $\lambda_i​$. According to Ferrenberg and Swendsen, the optimal weights can be found by solving the nonlinear set of equations given by
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

* The output weights are scaled such `sum(weights) = length(weights)`.
