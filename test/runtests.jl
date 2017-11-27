using FerrenbergSwendsenReweighting
using Base.Test

@testset "Single histogram reweighting test" begin include("SingleHistogramReweightingTest.jl") end
@testset "Multiple histogram reweighting test" begin include("MultipleHistogramReweightingTest.jl") end
