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



module FerrenbergSwendsenReweighting

using ArgCheck
export ReweightObj, evaluate, evaluate!

struct ReweightObj{F,T1,T2}
    logprob::F #The probability weight of sample x[i] is exp(logprob(λ,x[i]))
    δlogprob::T1 #Vector of shifts in logprob
    x::T2 #Input values/samples

    function ReweightObj{F,T1,T2}(logprob::F, δlogprob::T1, x::T2) where {F<:Function, T1<:AbstractVector, T2}
        @argcheck length(δlogprob) == length(x)
        new(logprob, δlogprob, x)
    end
end

Base.length(rw::ReweightObj) = length(rw.δlogprob)

function evaluate!(rw::ReweightObj{F,T1,T2}, λ, w::T1) where {F,T1,T2}
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

function evaluate(rw::ReweightObj, λ)
    w = similar(rw.δlogprob)
    evaluate!(rw, λ, w)
    return w
end



##################################
include("SingleHistogramReweighting.jl")
include("MultipleHistogramReweighting.jl")


##################################
#Default type is Float64!
ReweightObj(args...; kwargs...) = ReweightObj{Float64}(args...; kwargs...)

"Default reweighting to Boltzmann distribution, i.e. assuming x0=E, λ=-β"
boltzmann(β,E) = -β*E
ReweightObj{T}(args...; kwargs...) where T = ReweightObj{T}(boltzmann, args...; kwargs...)

end # module
