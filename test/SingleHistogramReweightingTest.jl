using FerrenbergSwendsenReweighting
using Base.Test

#Basic Boltzmann weight reweighting:
@test begin
    λ0 = 0.0
    x0 = randn(100)
    λ = 1.0
    typeof(reweights(λ0, x0, λ)) <: SingleHistogramReweights
end


#Float32 output
@test begin
    λ0 = 0.0
    x0 = randn(100)
    λ = 1.0
    rw = reweights(λ0, x0, λ, WeightType=Float32)
    eltype(rw) == Float32
end


#Different weight function:
@test begin
    λ0 = 0.0
    x0 = randn(100)
    λ = 1.0
    logprob(λx, x) = λx*x^2
    reweights(logprob, λ0, x0, λ)
    true
end
