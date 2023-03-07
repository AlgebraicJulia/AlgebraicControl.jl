using Test

#=@testset "Optimization Algorithms" begin
    include("Optimizers.jl")
end=#

@testset "Linear Quadratic Regulators" begin
    include("LQR.jl")
end
