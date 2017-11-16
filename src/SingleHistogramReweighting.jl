struct SingleHistogramReweights{S<:Real, T<:Real, V<:AbstractVector{T}} <: AbstractWeights{S, T, V}
    values::V
    sum::S
end

SingleHistogramReweights(vs::V, s::S=sum(vs)) where {S<:Real, V<:AbstractVector{<:Real}} = SingleHistogramReweights{S, eltype(vs), V}(vs, s)



"Find single histogram reweighting weights"
function reweights(logprob::Function, λ0, x0::AbstractVector, λ; WeightType=Float64)
    w = Vector{WeightType}(length(x0)) #Allocate weight storage.    AVOID EXPLICIT Float64 HERE?
    reweights!(logprob, w, x0, λ0, λ)
end

"In-place single histogram reweighting weights"
function reweights!(logprob::Function, w::Vector{<:AbstractFloat}, λ0, x0::AbstractVector, λ)
    Δlogprob(xi) = logprob(λ, xi) - logprob(λ0, xi)
    w .= Δlogprob.(x0)
    mw = maximum(w)
    w .= exp.(w .- mw) #Make sure the weights are maximum 1 (avoid overflow)
    w = SingleHistogramReweights(w) #NB: Is this the correct weight type to use here??
    return w
end
reweights!(logprob::Function, rw::SingleHistogramReweights, args...) = reweights!(logprob, rw.values, args...) #Re-use old weight storage
