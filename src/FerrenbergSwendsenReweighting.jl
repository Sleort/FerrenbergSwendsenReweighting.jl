__precompile__()

#=
INPUT:
    x: Vector of data (may themselves be vectors) OR xs: A vector of such vectors of data
    λ0: Parameter at which data has been sampled  OR λs: A vector of such parameters
    logprob: λ,x[i] → R; log probability of datapoint at some parameter λ. Default is -λ*x[i]; the Boltzmann weight

    τint: Integrated autocorrelation time
OUTPUT:
    reweighting object, which can be evaluated at a parameter of choice

    TODO:
        * Check order of arguments
=#

##########################################
include("AutocorrelationTime.jl")
##########################################




module FerrenbergSwendsenReweighting

using ArgCheck
export Reweights

struct Reweights{T1,T2}
    logprob::Function #exp(logprob(λ,x[i])) is the probability weight of sample x[i]
    δlogprob::T1 #Vector of shifts in logprob
    x::T2 #Input values/samples
    w::T1 #Probability weight storage (normalized!)

    function Reweights{T1,T2}(logprob::Function, δlogprob::T1, x::T2, w::T1) where {T1<:AbstractVector, T2}
        @argcheck length(δlogprob) == length(x) == length(w)
        new(logprob, δlogprob, x, w)
    end
end


function Reweights{T}(logprob::Function, δlogprob::T1, x::T2, w::T1) where {T<:Real, T1<:AbstractVector{T}, T2}
    Reweights{T1,T2}(logprob, δlogprob, x, w)
end


#Default type is Float64!
Reweights(args...; kwargs...) = Reweights{Float64}(args...; kwargs...)


"Default reweighting to Boltzmann distribution, i.e. assuming x0=E, λ=-β"
boltzmann(β,E) = -β*E
Reweights{T}(args...; kwargs...) where T = Reweights{T}(boltzmann, args...; kwargs...)


"""Return the "reweights" at parameter λ"""
function (rw::Reweights)(λ)
    for (i,xi) ∈ enumerate(rw.x)
        rw.w[i] = rw.logprob(λ, xi) + rw.δlogprob[i]
    end
    mw = maximum(rw.w)
    rw.w .= exp.(rw.w .- mw) #Make sure the weights are maximum 1 (avoid overflow)
    normalize!(rw.w, 1) #sum(w) = 1
    return rw.w
end


##################################
include("SingleHistogramReweighting.jl")
include("MultipleHistogramReweighting.jl")
##################################


end # module
