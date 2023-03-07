module ParaProbs

using ..Problems

abstract type ParaProb <: Problem end

struct ParaECProb <: ParaProb
    dom::Int
    codom::Int
    param::Int
    f::Function
    h::Function
end

function open_prob(f::SmoothFunction, h::SmoothFunction)
    

end