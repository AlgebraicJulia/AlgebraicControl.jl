module LensCats
export LensCat, Lens

using ..Categories
using ..MonoidalCats
import ..Categories: dom, codom, compose, id
import ..MonoidalCats: otimes, munit

⋅ = compose

struct Lens{Ob, Hom}
    dom::Pair{Ob, Ob}
    codom::Pair{Ob, Ob}
    get::Hom
    put::Hom
end

struct LensCat{Ob, Hom} <: MonoidalCategory{Pair{Ob, Ob}, Lens{Ob, Hom}}
    base::CartesianCategory{Ob, Hom}
end

dom(c::LensCat, f::Lens) = f.dom
codom(c::LensCat, f::Lens) = f.codom

id(c::LensCat, x::Pair) = Lens(x, x, id(c.base, x[1]), proj2(c.base, x[1], x[2]))

compose(c::LensCat, f::Lens, g::Lens) = begin
    @assert codom(c, f) == dom(c, g)
    C = c.base
    get = compose(C, f.get, g.get)
    π0 = proj1(C, dom(c, f)[1], codom(c, g)[2])
    π1 = proj2(C, dom(c, f)[1], codom(c, g)[2])
    put = ⋅(C, pair(C, π0, ⋅(C, pair(C, ⋅(C, π0, f.get), π1), g.put)), f.put)
    return Lens(dom(c, f), codom(c, g), get, put)
end

otimes(c::LensCat, x::Pair, y::Pair) = otimes(c.base, x[1], y[1]) => otimes(c.base, x[2], y[2])
otimes(c::LensCat, f::Lens, g::Lens) = Lens(
    otimes(c, dom(c, f), dom(c, g)),
    otimes(c, codom(c, f), codom(c, g)),
    otimes(c.base, f.get, g.get),
    otimes(c.base, f.put, g.put)
)

munit(c::LensCat) = munit(c.base) => munit(c.base)

end