using Test

@testset "ParaConv" begin
    include("ParaConvCat.jl")
end

@testset "Convex Model Predictive Control" begin
    include("CMPC.jl")
end

#=@testset "Optimization Algorithms" begin
    include("Optimizers.jl")
end=#

#=@testset "Linear Quadratic Regulators" begin
    include("LQR.jl")
end=#
