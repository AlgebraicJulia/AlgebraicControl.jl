# Cartesian Category of types and (pure) julia functions
module JulCat

using ..Categories
import ..Categories: dom, codom, compose, id
using ..MonoidalCats
import ..MonoidalCats: otimes, munit, duplicate, delete, pair, proj1, proj2

struct TypedFunction{dom, codom}
    impl::Function
end
(f::TypedFunction{τ1, τ2})(x::τ1) where {τ1, τ2} = f.impl(x)

struct Jul <: Category{Type, TypedFunction} end

dom(::Jul, f::TypedFunction{τ1, τ2}) where {τ1, τ2} = τ1
codom(::Jul, f::TypedFunction{τ1, τ2}) where {τ1, τ2} = τ2


end