module BasicCats
export SmoothFunction, Smooth

using ..Categories
import ..Categories: dom, codom, compose, id
using ..MonoidalCats
import ..MonoidalCats: otimes, munit, duplicate, delete, pair, proj1, proj2


struct SmoothFunction
    dom::Int
    codom::Int
    impl::Function
end
(f::SmoothFunction)(x) = f.impl(x)

struct Smooth <: CartesianCategory{Int, SmoothFunction} end

dom(::Smooth, f::SmoothFunction) = f.dom
codom(::Smooth, f::SmoothFunction) = f.codom
compose(c::Smooth, f::SmoothFunction, g::SmoothFunction) = begin
    codom(c, f) == dom(c, g) || error("Smooth composition: domain mismatch")
    SmoothFunction(dom(c, f), codom(c, g), x -> g(f(x)))
end
id(::Smooth, x::Int) = SmoothFunction(x, x, x->x)

otimes(::Smooth, x::Int, y::Int) = x + y
otimes(c::Smooth, f::SmoothFunction, g::SmoothFunction) = begin
    n = dom(c, f)
    m = dom(c, g)
    SmoothFunction(n+m, otimes(c, codom(c, f), codom(c, g)),
        x -> vcat(f(x[1:n]), g(x[n+1:end]))
    )
end

munit(::Smooth) = 0

duplicate(::Smooth, x) = SmoothFunction(x, x+x, x->vcat(x,x))
delete(::Smooth, x) = SmoothFunction(x, 0, x->Float64[])

proj1(c::Smooth, n::Int, m::Int) = SmoothFunction(otimes(c, n, m), n, x -> x[1:n])
proj2(c::Smooth, n::Int, m::Int) = SmoothFunction(otimes(c, n, m), m, x -> x[n+1:end])

pair(c::Smooth, f::SmoothFunction, g::SmoothFunction) = begin
    @assert dom(c, f) == dom(c, g)
    SmoothFunction(dom(c, f), otimes(c, codom(c,f), codom(c,g)), 
        x -> vcat(f(x), g(x))
    )
end

end