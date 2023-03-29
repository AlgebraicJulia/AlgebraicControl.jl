module CMPC

export one_step_bifunction, MPC_bifunction

using JuMP
using Ipopt
using LinearAlgebra
using Convex
using ..ParaConvCat
using ..Categories



function one_step_bifunction(f::Function,g::Function,A,B)
    n = size(A)[1]
    n_B = size(B)[1]
    n == n_B || error("Incosistent dimensions in dynamics.")
    # Check if B is nx1
    if size(B) == (size(B)[1],)
        m = 1
    else
        m = size(B)[2]
    end
    impl = (u1,x1,x2) -> ConvexBifunction(
        f(u1[1],x1),
        vcat(g(u1[1],x1), x2 == A*x1 + B*u1[1])
    )
    return OpenParaConvexBifunction(n,n,[m],impl)
end

function one_step_bifunction(x_dim, u_dim, f::Function,g::Function,h::Function)
    impl = (u1,x1,x2) -> ConvexBifunction(
        f(u1[1],x1),
        vcat(g(u1[1],x1), x2 == h(u1[1],x1))
    )
    return OpenParaConvexBifunction(x_dim,x_dim,[u_dim],impl)
end

function MPC_bifunction(one_step::OpenParaConvexBifunction, N::Int)
    res = one_step
    for i in 1:N-2
        res = compose(ParaConv(), res, one_step)
    end
    return res
end


end

