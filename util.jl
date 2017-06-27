# Represent a universal choice dataset as a vector of sizes (lengths of choices)
# and a vector of all items in all choices.
mutable struct UniversalChoiceDataset
    sizes::Vector{Int64}
    choices::Vector{Int64}
end

# Probabilistic universal subset choice model.
#
# -z is a vector of length max_size with the probability of choosing
#  a size-k subset being z[k]
# -probs are the item probabilities
# -gammas is a vector of length max_size of normalization parameters
# -H is a vector of length max_size where each element
#  is a dictionary that maps a choice to a probability.
mutable struct UniversalChoiceModel
    z::Vector{Float64}
    probs::Vector{Float64}
    gammas::Vector{Float64}
    H::Vector{Dict{NTuple,Float64}}
end

# Read text data
function read_data(dataset::AbstractString)
    f = open(dataset)
    choices = Int64[]
    sizes = Int64[]
    for line in eachline(f)
        choice = [parse(Int64, v) for v in split(line)]
        append!(choices, choice)
        push!(sizes, length(choice))
    end
    return UniversalChoiceDataset(sizes, choices)
end

# iterator over choices
#
# This is only syntactic sugar.  It pre-allocates the entire size of the output.
# Unfortunately, a "real" iterator with a Channel is too slow.
function iter_choices(data::UniversalChoiceDataset)
    curr_ind = 1
    choice_vec = Vector{Vector{Int64}}()
    for size in data.sizes
        choice = data.choices[curr_ind:(curr_ind + size - 1)]
        push!(choice_vec, choice)
        curr_ind += size
    end
    return zip(data.sizes, choice_vec)
end

#=
# iterator over choices (SLOW!!)
function iter_choices(sizes::Vector{Int64}, choices::Vector{Int64})
    curr_ind = 1
    function temp(c::Channel{Tuple{Int64,Vector{Int64}}})
        for size in sizes
            choice = choices[curr_ind:(curr_ind + size - 1)]
            push!(c, (size, choice))
            curr_ind += size
        end
    end
    return Channel(temp, ctype=Tuple{Int64,Vector{Int64}})
end
=#
