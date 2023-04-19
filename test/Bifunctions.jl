using AlgebraicControl.BasicCats

#=@data RSpace begin
    Rn(Int)
    Oplus(RSpace,RSpace)
end=#

struct RSpace
    ns::Vector{Int}
end

struct SmoothFunction
    dom::RSpace
    codom::RSpace
    impl::Function
end
(f::SmoothFunction)(x::Vector) = begin
    inputs = []
    prev = 0
    for n in f.dom.ns
        if n==1
            push!(inputs, x[prev+n])
        else
            push!(inputs, x[prev+1:prev+n])
        end
        prev += n
    end
    return f.impl(inputs...)
end
#(f::SmoothFunction)(x::Float64) = f.impl()

struct AtomicBifunction
    f::SmoothFunction
    h::SmoothFunction
    g::SmoothFunction
end

function AtomicBifunction(f::Function)
end