using Convex

struct ConvexProgram
    objective::Convex.AbstractExpr
    constraints::Vector{Convex.Constraint}
end

struct OpenConvexProgram
    dom_var::Int # Dimension of the domain decision variable
    codom_var::Int # Dimension of the codomain decision variable
    internal_vars::Vector{Int} # Dimensions of any internal variables
    impl::Function # internal_vars × dom_var × codom_var ⇢ ConvexProgram
end

(F::OpenConvexProgram)(ps::Vector{Variable}, x::Variable, y::Variable) =
    F.impl(ps, x, y)

dom(F::OpenConvexProgram) = F.dom_var
codom(F::OpenConvexProgram) = F.codom_var

function close_program(F::OpenConvexProgram, input_val::Vector{Float64}, output_val::Vector{Float64}, ivs::Vector{Variable})::ConvexProgram
    x_var = Variable(length(input_val))
    fix!(x_var, input_val)
    y_var = Variable(length(output_val))
    fix!(y_var, output_val)
    return F(ivs, x_var, y_var)
end

function compose(F::OpenConvexProgram, G::OpenConvexProgram)
    @assert codom(F) == dom(G)
    y = Variable(codom(F))
    F_nparams = length(F.internal_vars)
    impl = (ps, x, z) -> begin
        if F_nparams > 0
            FCB = F(ps[1:F_nparams], x, y)
        else
            FCB = F(Variable[], x, y)
        end
        GCB = G(ps[F_nparams+1:end], y, z)
        return ConvexProgram(FCB.objective + GCB.objective, vcat(FCB.constraints, GCB.constraints))
    end
    return OpenConvexProgram(dom(F), codom(G), vcat(F.internal_vars, G.internal_vars), impl)
end

# Compose an endomorphism with itself n times
function compose(F::OpenConvexProgram, n::Int)
    @assert dom(F) == codom(F)
    res = F
    for i in 1:n-1
        res = compose(res, F)
    end
    return res
end
