__precompile__()

#=
INPUT:
    x: Vector of data (may themselves be vectors) OR xs: A vector of such vectors of data
    λ0: Parameter at which data has been sampled  OR λs: A vector of such parameters
    logprob: λ,x[i] → R; log probability of datapoint at some parameter λ. Default is -λ*x[i]; the Boltzmann weight

    τint: Integrated autocorrelation time
    This is not used?
OUTPUT:
    reweighting object, which can be evaluated at a parameter of choice

    TODO:
        * Check order of arguments
=#



module FerrenbergSwendsenReweighting

using ArgCheck
using Statistics
#import IterTools
export ReweightObj, evaluate, evaluate!, reweight

struct ReweightObj{F,T1,T2}
    logprob::F #The probability weight of sample x[i] is exp(logprob(λ,x[i]))
    δlogprob::T1 #Vector of shifts in logprob
    x::T2 #Input values/samples
    ΔlogZs::AbstractVector

    function ReweightObj{F,T1,T2}(logprob::F, δlogprob::T1, x::T2, ΔlogZs) where {F<:Function, T1<:AbstractVector, T2}
        @argcheck length(δlogprob) == length(x)
        new(logprob, δlogprob, x, ΔlogZs)
    end
end

Base.length(rw::ReweightObj) = length(rw.δlogprob)

#Helping function to find the ranges the input vectors make up
#Not very elegant, but it works...
#function chain_ranges(chain)
#    ranges = []
#    i = 1
#    for l ∈ length.(chain.xss) #Length of each chained component
#        j = i + l - 1
#        push!(ranges, i:j)
#        i = j + 1
#    end
#    return ranges
#end


# Evaluate function for single histogram analysis.
function evaluate!(rw::ReweightObj{F,T1,T2}, λ, w::T1) where {F,T1,T2<:Vector}
    #@argcheck length(collect(rw.x)) == length(w)
    
    for (i,xi) ∈ enumerate(rw.x)
        w[i] = rw.logprob(λ, xi) + rw.δlogprob[i]
    end
    mw = maximum(w)
    w .= exp.(w .- mw) #Make sure the weights are maximum 1 (avoid overflow)
    z = sum(w)/length(w)
    w ./= z #Normalized such that sum(w) = length(w)

    #if rw.x isa Base.Iterators.Flatten
    #    w = [w[r] for r ∈ chain_ranges(rw.x)] #Split the weights to a vector of vectors to match the original data!
    #end

    nothing
end

# Evaluate function for multi-histogram analysis.
function evaluate!(rw::ReweightObj{F,Array{T1, 1},Array{Array{T2, 1},1}}, λ::R, w::Array{Array{T3, 1}, 1}) where {F,R<:Real, T1, T2, T3}
    @argcheck length.(rw.x) == length.(w)

    N₀ = length(rw.x)
    mws = Array{typeof(rw.δlogprob[1][1]), 1}(undef, N₀)
    for (k, xₖ) = enumerate(rw.x)
        for (i, xₖᵢ) = enumerate(xₖ)
            w[k][i] = rw.logprob(λ, xₖᵢ) + rw.δlogprob[k][i]
        end
        mws[k] = maximum(w[k])
    end
    mw = maximum(mws)

    sₖ = [w[1][1] for k = 1:N₀]
    # We make sure the weights are maximum 1 (avoid overflow)
    # and calculate the sum of w[k]
    for (k, wₖ) = enumerate(w)
        wₖ .= exp.(wₖ.-mw)
        sₖ[k] = sum(wₖ)
    end

    z = sum(sum.(w))/sum(length.(w))
    # Normalize w s.t. sum(sum.(w)) = sum(length.(w))
    w .= [[wₖᵢ/z for wₖᵢ in wₖ] for wₖ in w]

    nothing
end



function evaluate(rw::ReweightObj, λ)
    w = [copy(xₖ) for xₖ in rw.x]
    evaluate!(rw, λ, w)
    w
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
