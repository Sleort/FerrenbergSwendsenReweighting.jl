"Single histogram reweighting object"
function Reweights{T}(logprob::Function, λ0, x::AbstractVector) where T<:Real
    δlogprob = T[-logprob(λ0, xi) for xi ∈ x]
    w = similar(δlogprob)
    return Reweights{T}(logprob, δlogprob, x, w)
end
