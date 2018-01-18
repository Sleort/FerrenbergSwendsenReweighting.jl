using NLsolve
import IterTools: chain
import MCMCDiagnostics: ess_factor_estimate #Effective sample size estimate factor


#Denominator of Ferrenberg Swendsen equations:
function FSdenom(logprobλx::Matrix, iλ::Integer, ix::Integer, esss::Vector, ΔlogZs::Vector) #All logZs are compared with logZ[λ[1]]
    δ = logprobλx[iλ,ix]-ΔlogZs[iλ-1] #Input iλ ≥ 2 #Shift due to ΔlogZs should start at 2...
    d = esss[1]*exp(logprobλx[1,ix] - δ) #Contribution from λ1
    for i = 2:size(logprobλx,1)
        d += esss[i]*exp(logprobλx[i,ix] - ΔlogZs[i-1] - δ)
    end
    return d
end


#In-place calculating for all λs at the same time...
function FSexpression(logprobλx::Matrix, essfx::Vector, esss::Vector, ΔlogZs::Vector)
    accumulated = zeros(ΔlogZs)
    nλ, nx = size(logprobλx)
    for ix = 1:nx
        for iλ = 2:nλ #NB: Starting from index 2! (First λ sets reference point)
            accumulated[iλ-1] += essfx[ix]/FSdenom(logprobλx, iλ, ix, esss, ΔlogZs)
        end
    end
    return accumulated
end


# find_δlogprob returns a vector with the weight scaling necessary for each sample in the (flattened) xs
function find_δlogprob(logprob::Function, λs::AbstractVector, xs::AbstractVector{<:AbstractVector}, essfs::AbstractVector;  WeightType=Float64)
    length(xs) == length(essfs) == length(λs) || error("The number of input dataseries, effectice sample size factors (essfs), and input parameters (λs) must match!")

    #Setup for solving the multiple histogram equations:
    ns = length.(xs)
    nλ = length(λs)
    esss = ns .* essfs #Effective Sample SizeS
    essfx = vcat(fill.(essfs, ns)...) #A vector of essfs, each matching an observation x
    x = chain(xs...)
    logprobλx = WeightType[logprob(λi,xj) for λi ∈ λs, xj ∈ x]


    function f!(residual, ΔlogZs)
        residual .= FSexpression(logprobλx, essfx, esss, ΔlogZs) .- 1
    end


    #Solving the equations, thus determining the relative ΔlogZs (compared to the first series)
    ΔlogZs = logprobλx[2:end,1] .- logprobλx[1,1] #Approximate initialization of logZ difference
    ΔlogZs = nlsolve(f!, ΔlogZs, autodiff=:forward, show_trace=true).zero


    #Then we can calculate the approriate weight shift for a new, arbitrary λ:
    logprobλx[2:end, :] .-= ΔlogZs #Shift the sampled logprobs according to the calculated ΔlogZs. ΔlogZ = 0 for the first series.
    logprobλx .+= log.(esss) #... also include prefactors
    logprobλx .-= maximum(logprobλx) #Remove an (arbitrary) constant factor
    δlogprob = WeightType.(log.(essfx ./ vec(sum(exp,logprobλx,1)))) #The "true weight" is now logprob(λ,x) + δlogprob[ix]

    return δlogprob
end



################################################
function essf_estimate(x::AbstractVector)
    mx = mean(x)
    δx = sqrt(eps())*randn(size(x)) #Small perturbation to avoid failure when all values of x are equal (e.g. when β = 0...). Shouldn't be a problem except when MC the sampling completely freezes... (In which case we have bigger problems!)
    ess_factor_estimate((x .- mx) + δx)[1]
end

"""
If no effectice sample size factors (essfs) are given:
    The essfs of the `autocorrelation_observable(λ,x)` function is used!
If no `autocorrelation_observable` is given, `autocorrelation_observable=logprob`
"""
function ReweightObj{T}(logprob::Function, λs::AbstractVector, xs::AbstractVector{<:AbstractVector};
                    autocorrelation_observable::Function = logprob,
                    essfs::AbstractVector{<:Real} = [essf_estimate([autocorrelation_observable(λs[i],x) for x ∈ xs[i]]) for i=1:length(λs)]) where T<:Real

    δlogprob = find_δlogprob(logprob, λs, xs, essfs; WeightType=T)
    x = chain(xs...)

    F = typeof(logprob)
    T1 = typeof(δlogprob)
    T2 = typeof(x)
    return ReweightObj{F,T1,T2}(logprob, δlogprob, x)
end
