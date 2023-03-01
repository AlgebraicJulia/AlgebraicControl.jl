module MonoidalCats
export MonoidalCategory, otimes, munit, CartesianCategory, duplicate, delete,
    proj1, proj2, pair

using ..Categories

abstract type MonoidalCategory{Ob, Hom} <: Category{Ob, Hom} end

function otimes(c::MonoidalCategory{Ob, Hom}, x::Ob, y::Ob)::Ob where {Ob, Hom}
    error("unimplemented")
end

function otimes(c::MonoidalCategory{Ob, Hom}, f::Hom, g::Hom)::Hom where {Ob, Hom}
    error("unimplemented")
end

function munit(c::MonoidalCategory{Ob, Hom})::Ob where {Ob, Hom}
    error("unimplemented")
end

abstract type CartesianCategory{Ob, Hom} <: MonoidalCategory{Ob, Hom} end

function duplicate(c::CartesianCategory{Ob, Hom}, x::Ob)::Hom where {Ob, Hom}
    error("unimplemented")
end

function delete(c::CartesianCategory{Ob, Hom}, x::Ob)::Hom where {Ob, Hom}
    error("unimplemented")
end

function pair(c::CartesianCategory{Ob, Hom}, f::Hom, g::Hom)::Hom where {Ob, Hom}
    error("unimplimented")
end

function proj1(c::CartesianCategory{Ob, Hom}, x::Ob, y::Ob)::Hom where {Ob, Hom}
    error("unimplemented")
end

function proj2(c::CartesianCategory{Ob, Hom}, x::Ob, y::Ob)::Hom where {Ob, Hom}
    error("unimplemented")
end

end