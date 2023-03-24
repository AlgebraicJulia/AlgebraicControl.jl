module ParaConvCat

export ConvexBifunction, to_problem, OpenParaConvexBifunction, ParaConv

using ..Categories
import ..Categories: dom, codom, id, compose
using Convex

#=struct ParaFunction
    dom::Int
    codom::Int
    para::Int
    impl::Function
end

struct ParaConvexBifunction
    dom::Int
    codom::Int
    para::Int
    f::ParaFunction # f(u,x1,x2) Objective function
    g::ParaFunction # g(u,x1,x2) <= 0 Inequality constraints
    h::ParaFunction # h(u,x1,x2) == 0 Equality constraints
end

struct ParaConv <: Category{Int, ConvexBifunction} end

dom(::ParaConv, F::ParaConvexBifunction) = F.dom
codom(::ParaConv, F::ParaConvexBifunction) = F.codom
compose(c::ParaConv, F::ParaConvexBifunction, G::ParaConvexBifunction) = begin
    codom(c, F) == dom(c, G) || error("Domain mismatch in bifunction composition")

end=#



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
(F::OpenParaConvexBifunction)(ps::Vector{Variable}, x::Variable, y::Variable) =
    F.impl(ps, x, y)

struct ParaConv <: Category{Int, OpenParaConvexBifunction} end

dom(::ParaConv, F::OpenParaConvexBifunction) = F.dom
codom(::ParaConv, F::OpenParaConvexBifunction) = F.codom
id(::ParaConv, X::Int) = begin
    impl = (_, x1, x2) ->
        ConvexBifunction(Constant(0), [x1 == x2])
    return OpenParaConvexBifunction(X, X, 0, impl)
end
compose(::ParaConv, F::OpenParaConvexBifunction, G::OpenParaConvexBifunction) = begin
    @assert codom(ParaConv(), F) == dom(ParaConv(), G)
    y = Variable(codom(ParaConv(), F))
    F_nparams = length(F.para)
    G_nparams = length(G.para)
    impl = (ps, x, z) -> begin
        FCB = F(ps[1:F_nparams], x, y)
        GCB = G(ps[F_nparams+1:end], y, z)
        return ConvexBifunction(FCB.obj + GCB.obj, vcat(FCB.cons, GCB.cons))
    end
    return OpenParaConvexBifunction(F.dom, G.codom, vcat(F.para, G.para), impl)
end


end

