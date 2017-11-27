__precompile__()

#=
INPUT:
    x: Vector of data (may themselves be vectors) OR xs: A vector of such vectors of data
    λ0: Parameter at which data has been sampled  OR λs: A vector of such parameters
    logprob: λ,x[i] → R; log probability of datapoint at some parameter λ. Default is -λ*x[i]; the Boltzmann weight

    τint: Integrated autocorrelation time
OUTPUT:
    vector of weights of the data


    Options:
        * Truncate weight to zero if insignificant on some scale?
        * Type of weights (Float64 is standard)

=#

#=
    TODO:
        * Check order of arguments
        * Calculation of integrated autocorrelation time (own package?) - include it as a keyword argument
=#

##########################################
include("AutocorrelationTime.jl")
##########################################



module FerrenbergSwendsenReweighting

import StatsBase: AbstractWeights
export AbstractReweights, SingleHistogramReweights, MultipleHistogramReweights,
    reweights, reweights!

##################################
# Common for both single and multiple histograms:
##################################
boltzmann(β,E) = -β*E
"Default reweighting to Boltzmann distribution, i.e. assuming x0=E, λ=-β"
reweights(args...; kwargs...) = reweights(boltzmann, args...; kwargs...)
reweights!(args...; kwargs...) = reweights!(boltzmann, args...; kwargs...) #Default is Boltzmann distribution

include("SingleHistogramReweighting.jl")
include("MultipleHistogramReweighting.jl")
end # module
