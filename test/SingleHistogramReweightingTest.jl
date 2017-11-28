using FerrenbergSwendsenReweighting
using Base.Test

#Basic Boltzmann weight reweighting:
@test begin
    λ0 = 0.0
    x0 = randn(100)
    rw = Reweights(λ0, x0)
    rw(λ0) == ones(x0) #Flat distribution
end


#Float32 output
@test begin
    λ0 = 0.0
    x0 = randn(100)
    rw = Reweights{Float32}(λ0, x0)
    eltype(rw(1.0)) == Float32
end


#Different weight function:
@test begin
    λ0 = 0.0
    x0 = randn(100)
    logprob(λ, x) = λ*x^2
    rw = Reweights(logprob, λ0, x0)
    true
end
