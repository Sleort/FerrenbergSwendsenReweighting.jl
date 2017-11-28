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

    function Reweights{T1,T2}(logprob::Function, δlogprob::T1, x::T2) where {T1<:AbstractVector, T2}
        @argcheck length(δlogprob) == length(x)
        new(logprob, δlogprob, x)
    end
end



"""Return the "Reweights" at parameter λ"""
#In place calculation HOW TO MAKE THIS NICE LOOKING WITH EXCLAMATION MARK??
function (rw::Reweights{T1,T2})(w::T1, λ) where {T1,T2}
    @argcheck length(rw.x) == length(w)
    for (i,xi) ∈ enumerate(rw.x)
        w[i] = rw.logprob(λ, xi) + rw.δlogprob[i]
    end
    mw = maximum(w)
    w .= exp.(w .- mw) #Make sure the weights are maximum 1 (avoid overflow)
    z = sum(w)/length(w)
    w ./= z #Normalized such that sum(w) = length(w)
    return w
end

#Out-of place:
function (rw::Reweights)(λ)
    w = similar(rw.δlogprob)
    rw(w, λ)
end


##################################
include("SingleHistogramReweighting.jl")
include("MultipleHistogramReweighting.jl")


##################################
#Default type is Float64!
Reweights(args...; kwargs...) = Reweights{Float64}(args...; kwargs...)

"Default reweighting to Boltzmann distribution, i.e. assuming x0=E, λ=-β"
boltzmann(β,E) = -β*E
Reweights{T}(args...; kwargs...) where T = Reweights{T}(boltzmann, args...; kwargs...)

end # module
