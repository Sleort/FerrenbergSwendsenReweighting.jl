using FerrenbergSwendsenReweighting
using Base.Test


#Basic Boltzmann weight reweighting:
@test begin
    λ0 = 0.0
    x0 = randn(100)
    rw = ReweightObj(λ0, x0)
    evaluate(rw,λ0) == ones(x0) #Flat distribution
end


#Float32 output
@test begin
    λ0 = 0.0
    x0 = randn(100)
    rw = ReweightObj{Float32}(λ0, x0)
    w = evaluate(rw, 1.0)
    eltype(w) == Float32
end


#Different weight function:
@test begin
    λ0 = 0.0
    x0 = randn(100)
    logprob(λ, x) = λ*x^2
    rw = ReweightObj(logprob, λ0, x0)
    true
end


#Check length function
@test begin
    λ0 = 0.0
    x0 = randn(100)
    rw = ReweightObj(λ0, x0)
    length(rw) == 100
end


#In-place evaluation:
@test begin
    λ0 = 0.0
    x0 = randn(100)
    rw = ReweightObj(λ0, x0)
    w = zeros(length(rw))
    evaluate!(rw, λ0, w)
    w == ones(x0) #Flat distribution
end
