using Base.Threads

# Represent a variable choice dataset as a vector of slate sizes, a vector of
# all slates, a vector of choice set sizes, and a vector of all choices.
mutable struct VariableChoiceDataset
    slate_sizes::Vector{Int64}
    slates::Vector{Int64}
    choice_sizes::Vector{Int64}
    choices::Vector{Int64}
end

# Utility-based variable choice model.
#
# -z is a vector of length max_size with the probability of choosing
#  a size-k subset being z[k]
# -utilities are the item utilities
# -H is a vector of length max_size where each element
#  is a dictionary that maps a choice set to a utility.
mutable struct VariableChoiceModel
    z::Vector{Float64}
    utilities::Vector{Float64}
    H::Vector{Dict{NTuple,Float64}}
end

# Read text data
function read_data(dataset::AbstractString)
end

function in_hotset(model::VariableChoiceModel, choice::Vector{Int64})
    choice_size = length(choice)
    key = NTuple{choice_size, Int64}(choice)
    return haskey(model.H[choice_size], key)
end

function hotset_utility(model::VariableChoiceModel, choice::Vector{Int64})
    choice_size = length(choice)
    key = NTuple{choice_size, Int64}(choice)
    return model.H[choice_size][key]
end

function log_likelihood(model::VariableChoiceModel, data::VariableChoiceDataset)
end

function add_to_hotset(model::VariableChoiceModel, choice_to_add::Vector{Int64})
    if in_hotset(choice_to_add, model); error("Choice already in hot set."); end
    lc = length(choice_to_add)
    choice_tup = NTuple{lc, Int64}(choice_to_add)
    model.H[lc][choice_tup] = 0.0
end

function initialize_model(data::VariableChoiceDataset)
end
