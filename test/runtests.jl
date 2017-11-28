using FerrenbergSwendsenReweighting
using Base.Test

@testset "Single histogram reweighting tests" begin include("SingleHistogramReweightingTest.jl") end
@testset "Multiple histogram reweighting tests" begin include("MultipleHistogramReweightingTest.jl") end
