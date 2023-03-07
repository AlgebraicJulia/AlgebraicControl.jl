module LQR

export generating_problem, build_LQR_opt, LQR_optimizer, init_prob

using ..Problems
using ..BasicCats
using ..Optimizers

function LQR_optimizer(A, B, Q, R, x₀, N)
    p0 = init_prob(A, B, Q, R, x₀)
    p = generating_problem(A, B, Q, R)
    o = build_LQR_opt(p, N-1)
    return compose(o, OpenOptimizer(p0))
end

function init_prob(A, B, Q, R, x₀)
    n_x = size(A)[1]
    n_u = size(R)[1]
    # This is a really stupid hack. Idk why Julia works like This
    if n_u == 1
        R = R[1]
    end
    f(u) = begin
        return u'*R*u
    end
    h(u) = begin
        if n_u == 1
            return A*x₀ + B*u[1]
        else
            return A*x₀ + B*u
        end
    end
    obj = SmoothFunction(n_u, 1, f)
    cons = SmoothFunction(n_u, n_x, h)
    return EqConstrainedProb(obj, cons)

end


function generating_problem(A, B, Q, R)
    n_x = size(A)[1]
    n_u = size(R)[1]
    # This is a really stupid hack. Idk why Julia works like This
    if n_u == 1
        R = R[1]
    end
    f(in) = begin
        x = in[1:n_x]
        u = in[n_x+1:end]
        return x'*Q*x + u'*R*u
    end
    h(in) = begin
        x = in[1:n_x]
        u = in[n_x+1:end]
        if n_u == 1
            u = u[1]
        end
        return A*x + B*u
    end
    obj = SmoothFunction(n_x + n_u, 1, f)
    cons = SmoothFunction(n_x + n_u, n_x, h)
    return EqConstrainedProb(obj, cons)
end

function build_LQR_opt(p::EqConstrainedProb, N::Int)
    n_u = p.h.dom - p.h.codom
    println(n_u)
    o = OpenOptimizer(p; para=n_u)
    return CompositeOptimizer(repeat([o], N))
end

end