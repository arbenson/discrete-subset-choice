vec2ntuple(vec::Vector{Int64}) = NTuple{length(vec), Int64}(vec)

function index_points(sizes::Vector{Int64})
    inds = cumsum(sizes) + 1
    unshift!(inds, 1)
    return inds
end
