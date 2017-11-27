"Single histogram reweighting object"
function Reweights{T}(logprob::Function, λ0, x::AbstractVector) where T<:Real
    δlogprob = T[-logprob(λ0, xi) for xi ∈ x]
    T1 = typeof(δlogprob)
    T2 = typeof(x)
    return Reweights{T1,T2}(logprob, δlogprob, x)
end
