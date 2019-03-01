using NLsolve
#import IterTools: chain
import MCMCDiagnostics: ess_factor_estimate #Effective sample size estimate factor


# Given elements ln(aᵢ) in a vector. Calculates the sum ln( Σᵢaᵢ) by using the identity ln(a + b) = ln(a) + ln(1 + exp( ln(b) - ln(a) ))
# For this to be accurate it is wise to have a > b, we thus sort the list in decreasing order before calculating the sum.
function logSum(ln_elements::Vector; sorted=false)
    if !sorted
        sort!(ln_elements, rev=true)
    end
    ln_sum = ln_elements[1]
    for i = 2:length(ln_elements)
        ln_sum += log( 1 + exp( ln_elements[i] - ln_sum ) )
    end
    ln_sum
end


# Uses the energies at temperature βₘ₋₁ to estimate the difference Δᵐₘ₋₁lnZ between the partition functions at
# βₘ and βₘ₋₁. Then we use this to estimate all ΔᵐlnZ = Δᵐ₁lnZ.
function singleHistogramGuess(logprobλx::Vector{Vector{Vector{R}}}, Nₖ::Vector{I}, N₀::I) where {R<:Real, I<:Int}
    ΔlnZs = Array{R, 1}(undef, N₀-1) # We don't save for m = 1 since Δ¹₁lnZ = 0
    ln_series = logprobλx[2][1] .- logprobλx[1][1]
    ΔlnZs[2-1] = logSum(ln_series) - log(Nₖ[1])
    for m = 3:N₀
        ln_series = logprobλx[m][m-1] .- logprobλx[m-1][m-1]
        ΔlnZs[m-1] = ΔlnZs[m-2] + logSum(ln_series) - log(Nₖ[m-1])
    end
    ΔlnZs
end


function initialGuess(logprobλx::Vector{Vector{Vector{R}}}, Nₖ::Vector{I}, N₀::I) where {R<:Real, I<:Int}
    #[logprobλx[m][1][1] - logprobλx[1][1][1] for m = 2:N₀] # Bojesen guess
    #[1.0 for m = 2:N₀]
    singleHistogramGuess(logprobλx, Nₖ, N₀)
end

function initialGuess(logprob::Function, λs::AbstractVector, xs::AbstractVector{<:AbstractVector}; WeightType=Float64)
    N₀ = length(λs)
    Nₖ = length.(xs)
    logprobλx = [[[WeightType(logprob(λₘ, xs[k][i])) for i = 1:Nₖ[k]] for k = 1:N₀] for λₘ in λs]
    initialGuess(logprobλx, Nₖ, N₀)
end


# Denominator of Ferrenberg Swendsen equations. A constant factor χ is used to make sure exponents are < 0.
function FSdenom(logprobλx::Vector{Vector{Vector{R}}}, ΔlogZs::Vector{T}, eff_lengths::Vector{R}, χ::T2, k::I, 
                 i::I, N₀::I) where {T, T2, R<:Real, I<:Int} #All logZs are compared with logZ[λ[1]]

    d = eff_lengths[1]*exp(logprobλx[1][k][i] - χ) #Contribution from λ1
    for l = 2:N₀
        d += eff_lengths[l]*exp(logprobλx[l][k][i] - ΔlogZs[l-1] - χ)
    end
    return d
end


#In-place calculating for all λs at the same time...
function FSexpression(logprobλx::Vector{Vector{Vector{R}}}, ΔlogZs::Vector{T}, essfs::Vector{R}, eff_lengths::Vector{R},
                      Nₖ::Vector{I}, N₀::I) where {T, R<:Real, I<:Int}

    accumulated = similar(ΔlogZs)
    accumulated .= 0.0
    for k = 1:N₀
        for i = 1:Nₖ[k]
            # First we need to find the constant factor
            χ = max(logprobλx[1][k][i], maximum([logprobλx[m][k][i] - ΔlogZs[m-1] for m = 2:N₀]))
            # Then we use this in calculating the denominator
            d = FSdenom(logprobλx, ΔlogZs, eff_lengths, χ, k, i, N₀)
            # And finally add the contributions for the different ΔlogZs[m]s
            for m = 2:N₀
                accumulated[m-1] += essfs[m]*exp(logprobλx[m][k][i] - ΔlogZs[m-1]-χ)/d
            end
        end
    end
    accumulated
end


# Calculates the denominator of the FS equations.
function ΔlnZDenom(ΔlnZs::Vector{T}, lnprobλx::Vector{Vector{Vector{R}}}, ln_eff_lengths::Vector{R}, N₀::I, k::I, i::I) where {T, R<:Real, I<:Int}
    
    # First we find the logarithm of all elements
    ln_d = Array{T}(undef, N₀)
    ln_d[1] = ln_eff_lengths[1] + lnprobλx[1][k][i]
    for l = 2:N₀
        ln_d[l] = ln_eff_lengths[l] + lnprobλx[l][k][i] - ΔlnZs[l-1]
    end
    logSum(ln_d)
end


# Calculates the inner sum of the double sum in the FS equations for every ΔlnZ
function ΔlnZInnerSumSeries(ΔlnZs::Vector{T}, lnprobλx::Vector{Vector{Vector{R}}}, ln_eff_lengths::Vector{R}, N₀::I, 
                            Nₖ::Vector{I}, k::I) where {T, R<:Real, I<:Int}

    # Storage for the logarithms of the terms in the sum. There is one sum pr. m
    ln_sᵏₘᵢ = [Array{T, 1}(undef, Nₖ[k]) for m = 1:N₀-1]
    
    for i = 1:Nₖ[k]
        ln_denom = ΔlnZDenom(ΔlnZs, lnprobλx, ln_eff_lengths, N₀, k, i)
        
        # Now that we have the energy and the denominator we can store the logarithm for each m
        for m = 2:N₀
            ln_sᵏₘᵢ[m-1][i] = lnprobλx[m][k][i] - ln_denom
        end
    end
    
    # Now each vector ln_sᵏₘᵢ[m] is a series of (unsorted) logarithms. 
    ln_sᵏₘᵢ
end


# Calculates an estimate of ΔlnZs based on a previous guess of ΔlnZs. This should converge eventually.
function ΔlnZ(ΔlnZs::Vector{T}, lnprobλx::Vector{Vector{Vector{R}}}, ln_gₖ::Vector{R}, ln_eff_lengths::Vector{R},
              Nₖ::Vector{I}, N₀::I) where {T, R<:Real, I<:Int}

    new_ΔlnZs = Array{T}(undef, N₀-1)
    # Storage for outer sum-series
    ln_ssₘₖ = [Array{T, 1}(undef, N₀) for m = 1:N₀-1]
    for k = 1:N₀
        # Calculate inner sum-series for each m
        ln_sᵏₘᵢ = ΔlnZInnerSumSeries(ΔlnZs, lnprobλx, ln_eff_lengths, N₀, Nₖ, k)
        
        for m = 2:N₀
            ln_ssₘₖ[m-1][k] = ln_gₖ[k] + logSum(ln_sᵏₘᵢ[m-1])
        end
    end
    
    # Now the vector ln_ssₘₖ[m] is a series of logarithms of each term in the k-sum.
    for m = 2:N₀
        new_ΔlnZs[m-1] = logSum(ln_ssₘₖ[m-1])
    end
    new_ΔlnZs
end


# Sets up the Ferrenberg Swendsen equations for ΔᵐlnZᵦ = lnZᵦₘ - lnZᵦ₁ for lnZᵦ = -fᵦ, either on a fixed point form
# or on a kernel form depending on whether logarithm keywork is true or not. Then solves this equation using the
# Nlsolve package for solving non-linear equations.
function solveFSEquations(logprobλx::Vector{Vector{Vector{T}}}, essfs::Vector{R}, eff_lengths::Vector{R}, ln_eff_lengths::Vector{R},
                          Nₖ::Vector{I}, N₀::I; iterations=1000, verbose=true, logarithms=true, 
                          initial_guess=initialGuess(logprobλx, Nₖ, N₀)) where {T, R<:Real, I<:Int}
    ln_gₖ = log.(essfs)
    
    function f!(residual, ΔlogZs)
        residual .= FSexpression(logprobλx, ΔlogZs, essfs, eff_lengths, Nₖ, N₀) .- 1
    end
    function g!(residual, ΔlogZs)
        residual .= ΔlnZ(ΔlogZs, logprobλx, ln_gₖ, ln_eff_lengths, Nₖ, N₀) .- ΔlogZs
    end

    #Solving the equations, thus determining the relative ΔlogZs (compared to the first series)
    t = 0
    ΔlogZs = initial_guess
    if logarithms
        ΔlogZs, t, _, _, _ = @timed nlsolve(g!, initial_guess; autodiff=:forward, show_trace=verbose, iterations=iterations).zero
    else
        ΔlogZs, t, _, _, _ = @timed nlsolve(f!, initial_guess; autodiff=:forward, show_trace=verbose, iterations=iterations).zero
    end
    if verbose
        println("Calculated ΔlogZs:\n$(ΔlogZs)")
        println("In $(t/60) m")
    end
    ΔlogZs
end

function maxVal(A::Vector{Vector{Vector{T}}}) where T<:Number
    mx = A[1][1][1]
    for Aₘ ∈ A
        for Aₘₖ ∈ Aₘ
            for Aₘₖᵢ ∈ Aₘₖ
                mx = max(mx, Aₘₖᵢ)
            end
        end
    end
    mx
end

# This function mutates the logprobλx arrays to save memory.
# find_δlogprob returns a vector with the weight scalings necessary for each sample in xs.
function find_δlogprob(logprobλx::Vector{Vector{Vector{T}}}, ΔlogZs::Vector{T1}, essfs::Vector{R}, ln_eff_lengths::Vector{R},
                       N₀::I, Nₖ::Vector{I}; WeightType=Float64) where {T, T1, R<:Real, I<:Int}

    for k = 1:N₀
        for i = 1:Nₖ[k]
            logprobλx[1][k][i] += ln_eff_lengths[1]
            for m = 2:N₀
                logprobλx[m][k][i] += -ΔlogZs[m-1] + ln_eff_lengths[m]
            end
        end
    end
    mx = maxVal(logprobλx)
    for m = 1:N₀
        for k = 1:N₀
            logprobλx[m][k] .-= mx
        end
    end
#    @test maxVal(logprobλx) == 0.0

    # Finally calculate the shifts such that the "true weight" is logprob(λ, x) + δlogprob[ix]
    t1 = typeof(logprobλx[1][1][1])
    δlogprob = Array{Array{WeightType, 1}, 1}(undef, N₀)
    for k = 1:N₀
        δlogprob[k] = Array{WeightType, 1}(undef, Nₖ[k])
        for i = 1:Nₖ[k]
            ln_series = [logprobλx[m][k][i] for m = 1:N₀]
            δlogprob[k][i] = log(essfs[k]) - logSum(ln_series)
        end
    end

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
                    essfs::AbstractVector{<:Real} = [essf_estimate([autocorrelation_observable(λs[i],x) for x ∈ xs[i]]) for i=1:length(λs)],
                    iterations=1000, verbose=true, simple=false, logarithms=true, 
                    initial_guess=initialGuess(logprob, λs, xs; WeightType=T)) where T<:Real

    length(xs) == length(essfs) == (N₀ = length(λs)) || error("The number of input dataseries, effective sample size factors (essfs), and input parameters (λs) must match!")
    

    δlogprob = [Array{T, 1}(undef, length(xs[k])) for k = 1:N₀]
    ΔlogZs = Array{T, 1}(undef, N₀-1)
    if !simple

        #Setup for solving the multiple histogram equations:
        Nₖ = length.(xs)
        eff_lengths = Nₖ .* essfs #Effective Sample SizeS
        ln_eff_lengths = log.(eff_lengths)

        #logprobλx = WeightType[logprob(λi,xj) for λi ∈ λs, xj ∈ vcat(xs...)]
        logprobλx = [[[T(logprob(λₘ, xs[k][i])) for i = 1:Nₖ[k]] for k = 1:N₀] for λₘ in λs]

        ΔlogZs = solveFSEquations(logprobλx, essfs, eff_lengths, ln_eff_lengths, Nₖ, N₀;
                                  iterations=iterations, verbose=verbose, logarithms=logarithms, initial_guess=initial_guess)

        δlogprob .= find_δlogprob(logprobλx, ΔlogZs, essfs, ln_eff_lengths, N₀, Nₖ; WeightType=T)

        #δlogprob .= find_δlogprob(logprob, λs, xs, essfs; WeightType=T, iterations=iterations, verbose=verbose, logarithms=logarithms)
    end

    #x = chain(xs...)
    #x = Base.Iterators.flatten(xs)

    F = typeof(logprob)
    T1 = typeof(δlogprob)
    T2 = typeof(xs)
    return ReweightObj{F,T1,T2}(logprob, δlogprob, xs, ΔlogZs)
end

# Given a set of series of observables that corresponds to the set of series of energies used in constructing the
# ReweightObj rw, this function generates the multi-histogram reweighted observables at new couplings λ in rwt_λs
function reweight(O_by_λ::Vector{Vector{T}}, rw::ReweightObj, rwt_λs::Vector{R}) where {T, R<:Real}
    ot = typeof(O_by_λ[1][1])
    rwt_obs = fill(ot(0), length(rwt_λs))
    
    for (i, λ) = enumerate(rwt_λs)
        weights = evaluate(rw, λ)
        for k = 1:length(O_by_λ)
            rwt_obs[i] += sum(weights[k] .* O_by_λ[k])
        end
        rwt_obs[i] /= sum(sum.(weights))
    end
    
    rwt_obs
end
