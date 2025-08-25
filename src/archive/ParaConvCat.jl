module ParaConvCat

export ConvexBifunction, to_problem, OpenParaConvexBifunction, ParaConv, to_cvx

using ..Categories
import ..Categories: dom, codom, id, compose
using Convex

struct ConvexBifunction
    obj::Convex.AbstractExpr
    cons::Vector{Constraint}
end

function to_problem(F::ConvexBifunction)::Problem
    if F.obj isa Constant
        return satisfy(F.cons)
    else
        return minimize(F.obj, F.cons)
    end
end


struct OpenParaConvexBifunction
    dom::Int
    codom::Int
    para::Vector{Int}
    impl::Function #param_vars × dom_var × codom_var ⇢ ConvexBifunction
end


(F::OpenParaConvexBifunction)(ps::Vector{Variable}, x::Vector, y::Vector) = begin
    x_var = Variable(length(x))
    fix!(x_var, x)
    F(ps, x_var, y)
end

(F::OpenParaConvexBifunction)(ps::Vector{Variable}, x::Vector, y::Variable) = begin
    x_var = Variable(length(x))
    fix!(x_var, x)
    F(ps, x_var, y)
end

(F::OpenParaConvexBifunction)(ps::Vector{Variable}, x::Variable, y::Vector) = begin
    y_var = Variable(length(y))
    fix!(y_var, y)
    F(ps, x, y_var)
end


(F::OpenParaConvexBifunction)(ps::Vector{Variable}, x::Variable, y::Variable) =
    F.impl(ps, x, y)


function to_cvx(F::OpenParaConvexBifunction, args...)
    bf = F(args...)
    return to_problem(bf)
end

struct ParaConv <: Category{Int, OpenParaConvexBifunction} end

dom(::ParaConv, F::OpenParaConvexBifunction) = F.dom
codom(::ParaConv, F::OpenParaConvexBifunction) = F.codom
id(::ParaConv, X::Int) = begin
    impl = (_, x1, x2) ->
        ConvexBifunction(Constant(0), [x1 == x2])
    return OpenParaConvexBifunction(X, X, Variable[], impl)
end
compose(::ParaConv, F::OpenParaConvexBifunction, G::OpenParaConvexBifunction) = begin
    @assert codom(ParaConv(), F) == dom(ParaConv(), G)
    y = Variable(codom(ParaConv(), F))
    F_nparams = length(F.para)
    G_nparams = length(G.para)
    impl = (ps, x, z) -> begin
        if F_nparams > 0
            FCB = F(ps[1:F_nparams], x, y)
        else
            FCB = F(Variable[], x, y)
        end
        GCB = G(ps[F_nparams+1:end], y, z)
        return ConvexBifunction(FCB.obj + GCB.obj, vcat(FCB.cons, GCB.cons))
    end
    return OpenParaConvexBifunction(F.dom, G.codom, vcat(F.para, G.para), impl)
end


end

