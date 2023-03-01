module ParaCats
export ParaCat, ParaHom

using ..Categories
using ..MonoidalCats
import ..Categories: dom, codom, compose, id

struct ParaHom{Ob, Hom}
    P::Ob
    dom::Ob
    codom::Ob
    impl::Hom
end

struct ParaCat{Ob, Hom} <: Category{Ob, ParaHom{Ob,Hom}}
    base::MonoidalCategory{Ob, Hom}
end

dom(c::ParaCat, f::ParaHom) = f.dom
codom(c::ParaCat, f::ParaHom) = f.codom

compose(c::ParaCat, f::ParaHom, g::ParaHom) = begin
    codom(c, f) == dom(c, g) || error("ParaCat composition: domain mismatch")
    ParaHom(otimes(c.base, f.P, g.P),
        dom(c, f),
        codom(c, g),
        compose(c.base, otimes(c.base, id(c.base, f.P), f.impl), g.impl)
    )
end

id(c::ParaCat, x) = begin
    ParaHom(munit(c.base), x, x, id(c.base, x))
end

end