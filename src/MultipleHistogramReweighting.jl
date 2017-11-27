import IterTools: chain
using NLsolve
import AutocorrelationTime: integrated_autocorrelation_time
const τint = integrated_autocorrelation_time


#Denominator of Ferrenberg Swendsen equations:
function FSdenom(logprobλx::Matrix, iλ::Integer, ix::Integer, nginvs::Vector, ΔlogZs::Vector) #All logZs are compared with logZ[λ[1]]
    δ = logprobλx[iλ,ix]-ΔlogZs[iλ-1] #iλ ≥ 2 #Shift
    d = nginvs[1]*exp(logprobλx[1,ix] - δ) #Contribution from λ1
    for i = 2:size(logprobλx,1)
        d += nginvs[i]*exp(logprobλx[i,ix] - ΔlogZs[i-1] - δ)
    end
    return d
end


#In-place calculating for all λs at the same time...
function FSexpression(logprobλx::Matrix, ginvx::Vector, nginvs::Vector, ΔlogZs::Vector)
    accumulated = zeros(ΔlogZs)
    nλ, nx = size(logprobλx)
    for ix = 1:nx
        for iλ = 2:nλ #NB: Starting from index 2! (First λ sets reference point)
            accumulated[iλ-1] += ginvx[ix]/FSdenom(logprobλx, iλ, ix, nginvs, ΔlogZs)
        end
    end
    return accumulated
end


# find_δlogprob returns a vector with the weight scaling necessary for each sample in the (flattened) xs
function find_δlogprob(logprob::Function, λs::AbstractVector, xs::AbstractVector{<:AbstractVector}, τints;  WeightType=Float64)
    length(xs) == length(τints) == length(λs) || error("The number of input dataseries, integrated autocorrelation times (τints), and input parameters (λs) must match!")

    #Setup for solving the multiple histogram equations:
    ns = length.(xs)
    nλ = length(λs)
    ginvs = 1./(2.*τints)
    nginvs = ginvs .* ns
    ginvx = vcat(fill.(ginvs, ns)...) #A vector of ginvs, each matching an observation x
    x = chain(xs...)
    logprobλx = WeightType[logprob(λi,xj) for λi ∈ λs, xj ∈ x]

    function f!(ΔlogZs, fδ)
        fδ .= FSexpression(logprobλx, ginvx, nginvs, ΔlogZs) .- 1
    end

    #Solving the equations, thus determining the relative ΔlogZs (compared to the first series)
    ΔlogZs = zeros(WeightType, nλ-1)
    ΔlogZs = nlsolve(f!, ΔlogZs, autodiff = true, show_trace=true).zero

    #Then we can calculate the approriate weight shift for a new, arbitrary λ:
    logprobλx[2:end, :] .-= ΔlogZs #Shift the sampled logprobs according to the calculated ΔlogZs. ΔlogZ = 0 for the first series.
    logprobλx .+= log.(nginvs) #... also include prefactors
    logprobλx .-= maximum(logprobλx) #Remove an (arbitrary) constant factor
    δlogprob = WeightType.(log.(ginvx ./ vec(sum(exp,logprobλx,1)))) #The "true weight" is now logprob(λ,x) + δlogprob[ix]

    return δlogprob
end



################################################
"""
If no integrated autocorrelation times are given:
    The autocorrelation times of the `autocorrelation_observable(λ,x)` function is used!
If no `autocorrelation_observable` is given, `autocorrelation_observable=logprob`
"""
function Reweights{T}(logprob::Function, λs::AbstractVector, xs::AbstractVector{<:AbstractVector};
                    autocorrelation_observable::Function = logprob,
                    τints::AbstractVector{<:Real} = [τint([autocorrelation_observable(λs[i],x) for x ∈ xs[i]]) for i=1:length(λs)]) where T<:Real
    δlogprob = find_δlogprob(logprob, λs, xs, τints; WeightType=T)
    x = chain(xs...)

    T1 = typeof(δlogprob)
    T2 = typeof(x)
    return Reweights{T1,T2}(logprob, δlogprob, x)
end
