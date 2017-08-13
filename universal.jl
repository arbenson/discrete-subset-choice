include("common.jl")
using StatsBase
using Base.Threads

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
    H::Dict{NTuple,Float64}
    # Some data to make updating computations faster
    # how often each item appears not in a hot set    
    item_counts::Vector{Int64}
    # how many times each subset appears
    subset_counts::Dict{NTuple, Int64}
    # number of times choice set of each size appears
    size_counts::Vector{Int64}
end

function get_subset_counts(data::UniversalChoiceDataset)
    counts = Dict{NTuple, Int64}()
    inds = index_points(data.sizes)
    for i in 1:length(data.sizes)
        size = data.sizes[i]        
        choice = data.choices[inds[i]:(inds[i + 1] - 1)]
        choice_tup = vec2ntuple(choice)
        if !haskey(counts, choice_tup); counts[choice_tup] = 0; end
        counts[choice_tup] += 1
    end
    return counts
end

function subsets_and_counts(data::UniversalChoiceDataset)
    subsets = Vector{NTuple}()
    counts = Int64[]
    for (subset, count) in get_subset_counts(data)
        push!(subsets, subset)
        push!(counts, count)
    end
    return subsets, counts
end

# Read text data
function read_data(dataset::AbstractString)
    f = open(dataset)
    choices = Int64[]
    sizes = Int64[]
    for line in eachline(f)
        choice = [parse(Int64, v) for v in split(line)]
        # Note: all choices are sorted!
        sort!(choice)
        append!(choices, choice)
        push!(sizes, length(choice))
    end
    return UniversalChoiceDataset(sizes, choices)
end

# iterator over choices
#
# This is only syntactic sugar.  It allocates the entire size of the output.
# Unfortunately, a "real" iterator with a Channel is too slow.
function iter_choices(data::UniversalChoiceDataset)
    curr_ind = 1
    choice_vec = Vector{Vector{Int64}}()
    for size in data.sizes
        choice = data.choices[curr_ind:(curr_ind + size - 1)]
        push!(choice_vec, choice)
        curr_ind += size
    end
    assert(curr_ind == length(data.choices) + 1)
    return zip(data.sizes, choice_vec)
end

in_hotset(model::UniversalChoiceModel, choice::Vector{Int64}) = haskey(model.H, vec2ntuple(choice))
in_hotset(model::UniversalChoiceModel, choice::NTuple) = haskey(model.H, choice)
hotset_prob(model::UniversalChoiceModel, choice::Vector{Int64}) = model.H[vec2ntuple(choice)]
hotset_prob(model::UniversalChoiceModel, choice::NTuple) = model.H[choice]

function log_likelihood(model::UniversalChoiceModel, subsets::Vector{NTuple}, counts::Vector{Int64})
    ns = length(subsets)
    all_lls = zeros(Float64, ns)
    Threads.@threads for i = 1:ns
        choice = subsets[i]
        size = length(choice)
        ll = log(model.z[size])
        if in_hotset(model, choice)
            ll += log(hotset_prob(model, choice))
        else
            ll += log(model.gammas[length(choice)])
            for item in choice; ll += log(model.probs[item]); end
        end
        all_lls[i] = counts[i] * ll
    end
    return sum(all_lls)
end


function log_likelihood(model::UniversalChoiceModel,
                        data::UniversalChoiceDataset)
    subsets, counts = subsets_and_counts(data)
    return log_likelihood(model, subsets, counts)
end

function normalization_values(max_size::Int64, H::Dict{NTuple,Float64},
                              item_probs::Vector{Float64})
    gammas = zeros(Float64, max_size)

    # No hotsets of size 1
    gammas[1] = 1
    
    base_probs = zeros(Float64, max_size)
    hotset_probs = zeros(Float64, max_size)
    for (subset, val) in H
        ns = length(subset)
        base_probs[ns] += prod([item_probs[item] for item in subset])
        hotset_probs[ns] += val
    end

    if max_size == 1; return gammas; end
    # normalization for size-2 choice probabilities
    sum_pi2   = sum([p^2 for p in item_probs])
    sum_pipj  = (1 - sum_pi2) / 2
    gammas[2] = (1 - hotset_probs[2]) / (sum_pipj + sum_pi2 - base_probs[2])

    if max_size == 2; return gammas; end    
    # normalization for size-3 choice probabilities
    sum_pi3    = sum([p^3 for p in item_probs])
    sum_pipjpj = sum_pi2 - sum_pi3
    sum_pipjpk = (1 - sum_pi3 - 3 * sum_pipjpj) / 6
    gammas[3]  = (1 - hotset_probs[3]) / (sum_pipjpk + sum_pipjpj + sum_pi3 - base_probs[3])

    if max_size == 3; return gammas; end
    # normalization for size-4 choice probabilities
    sum_pi4      = sum([p^4 for p in item_probs])
    sum_pi2pj2   = (sum_pi2 ^ 2 - sum_pi4) / 2
    sum_pipj3    = sum_pi3 - sum_pi4
    sum_pipjpk2  = (sum_pi2 - sum_pi4 - 2 * sum_pipj3 - 2 * sum_pi2pj2) / 2    
    sum_pipjpkpl = (1 - sum_pi4 - 4 * sum_pipj3 - 6 * sum_pi2pj2 - 12 * sum_pipjpk2) / 24    
    gammas[4]    = (1 - hotset_probs[4]) / (sum_pipjpkpl + sum_pi2pj2 + sum_pipj3 + sum_pipjpk2 + sum_pi4 - base_probs[4])
    
    if max_size == 4; return gammas; end
    # normalization for size-5 choice probabilities
    sum_pi5        = sum([p^5 for p in item_probs])
    sum_pipj4      = sum_pi4 - sum_pi5
    sum_pi2pj3     = sum_pi2 * sum_pi3 - sum_pi5
    sum_pipjpk3    = (sum_pi3 - sum_pi5 - 2 * sum_pipj4 - sum_pi2pj3) / 2    
    sum_pipj2pk2   = (sum_pi2^2 - sum_pipj4 - 2 * sum_pi2pj3 - sum_pi5) / 2    
    sum_pipjpkpl2  = sum_pi2 * sum_pipjpk - sum_pipjpk3
    sum_pipjpkplpm = (1.0 - sum_pi5 - 5 * sum_pipj4 - 20 * sum_pipjpk3 - 60 * sum_pipjpkpl2 - 30 * sum_pipj2pk2 - 10 * sum_pi2pj3) / 120
    gammas[5]      = (1 - hotset_probs[5]) / (sum_pipjpkplpm + sum_pi5 + sum_pipj4 + sum_pi2pj3 + sum_pipjpk3 + sum_pipjpkpl2 + sum_pipj2pk2 - base_probs[5])
    
    # TODO: support larger sizes.  We really shouldn't hard code the above.
    # There is a triangular system we can solve to automatically determine these
    # values, but it is annoying to set up.
    if max_size > 5; error("Support only for choices of size <= 5 items"); end

    return gammas
end

function add_to_hotset!(model::UniversalChoiceModel, choice_to_add::Vector{Int64})
    if in_hotset(model, choice_to_add); error("Choice already in hot set."); end

    # Update hotset probability
    choice_tup = vec2ntuple(choice_to_add)
    choice_count = model.subset_counts[choice_tup]
    model.H[choice_tup] = choice_count / model.size_counts[length(choice_to_add)]

    # Update item counts and probability
    for item in choice_to_add
        model.item_counts[item] -= choice_count
    end
    
    total = sum(model.item_counts)
    Threads.@threads for i = 1:length(model.item_counts)
        model.probs[i] = model.item_counts[i] / total
    end

    # Update normalization parameters
    model.gammas = normalization_values(length(model.gammas), model.H, model.probs)
end

function remove_from_hotset!(model::UniversalChoiceModel, choice_to_rm::Vector{Int64})
    if !in_hotset(model, choice_to_rm); error("Choice not in hot set."); end

    # Update hotset probability
    choice_tup = vec2ntuple(choice_to_rm)
    delete!(model.H, choice_tup)

    # Update item counts and probability
    choice_count = model.subset_counts[choice_tup]
    for item in choice_to_rm; model.item_counts[item] += choice_count; end
    total = sum(model.item_counts)
    for (i, count) in enumerate(model.item_counts); model.probs[i] = count / total; end

    # Update normalization parameters
    model.gammas = normalization_values(length(model.gammas), model.H, model.probs)
end


function initialize_model(data::UniversalChoiceDataset)
    max_size = maximum(data.sizes)

    # fixed-size probability are just the empirical fraction of selections
    z = zeros(Float64, max_size)
    size_counts = zeros(Int64, max_size)
    for (key, value) in countmap(data.sizes)
        z[key] = value / length(data.sizes)
        size_counts[key] = value
    end

    # item probabilities are initialized to the empirical fraction of
    # appearances in choice sets
    item_counts = zeros(Int64, maximum(data.choices))
    for (size, choice) in iter_choices(data)
        for item in choice; item_counts[item] += 1; end
    end
    item_probs = zeros(Float64, length(item_counts))
    total = sum(item_counts)
    for (i, count) in enumerate(item_counts); item_probs[i] = count / total; end

    # Initialize empty hotsets
    H = Dict{NTuple,Float64}()

    # Initialize normalization constants (gammas)
    gammas = normalization_values(max_size, H, item_probs)

    # Counts of how many times each choice set appears
    subset_counts = get_subset_counts(data::UniversalChoiceDataset)

    return UniversalChoiceModel(z, item_probs, gammas, H, item_counts,
                                subset_counts, size_counts)
end
