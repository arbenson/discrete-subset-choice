using StatsBase

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
        sort!(choice)
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

function in_hotset(choice::Vector{Int64}, model::UniversalChoiceModel)
    choice_size = length(choice)
    key = NTuple{choice_size, Int64}(choice)
    return haskey(model.H[choice_size], key)
end

function hotset_prob(choice::Vector{Int64}, model::UniversalChoiceModel)
    choice_size = length(choice)
    key = NTuple{choice_size, Int64}(choice)
    return model.H[choice_size][key]
end

function log_likelihood(data::UniversalChoiceDataset, model::UniversalChoiceModel)
    ll = 0.0
    for (size, choice) in iter_choices(data)
        ll += log(model.z[size])
        if in_hotset(choice, model)
            ll += log(hotset_prob(choice, model))
        else
            ll += log(model.gammas[length(choice)])
            for item in choice; ll += log(model.probs[item]); end
        end
    end
    return ll
end

function normalization_values(max_size::Int64, H::Vector{Dict{NTuple,Float64}},
                              item_probs::Vector{Float64})
    gammas = zeros(Float64, max_size)

    # No hotsets of size 1
    gammas[1] = 1

    if max_size == 1; return gammas; end
    # normalization for size-2 choice probabilities
    sum_pi2       = sum([p^2 for p in item_probs])
    sum_pipj      = (1 - sum_pi2) / 2
    hotset_probs2 = sum([val for (key, val) in H[2]])
    gammas[2]     = (1 - hotset_probs2) / (sum_pipj + sum_pi2)

    if max_size == 2; return gammas; end    
    # normalization for size-3 choice probabilities
    sum_pi3       = sum([p^3 for p in item_probs])
    sum_pipjpj    = sum_pi2 - sum_pi3
    sum_pipjpk    = (1 - sum_pi3 - 3 * sum_pipjpj) / 6
    hotset_probs3 = sum([val for (key, val) in H[3]])
    gammas[3]     = (1 - hotset_probs3) / (sum_pipjpk + sum_pipjpj + sum_pi3)

    if max_size == 3; return gammas; end
    # normalization for size-4 choice probabilities
    sum_pi4       = sum([p^4 for p in item_probs])
    sum_pi2pj2    = (sum_pi2 ^ 2 - sum_pi4) / 2
    sum_pipj3     = sum_pi3 - sum_pi4
    sum_pipjpk2   = (sum_pi2 - sum_pi4 - 2 * sum_pipj3 - 2 * sum_pi2pj2) / 2    
    sum_pipjpkpl  = (1 - sum_pi4 - 4 * sum_pipj3 - 6 * sum_pi2pj2 - 12 * sum_pipjpk2) / 24    
    hotset_probs4 = sum([val for (key, val) in H[4]])
    gammas[4]     = (1 - hotset_probs4) / (sum_pipjpkpl + sum_pi2pj2 + sum_pipj3 + sum_pipjpk2 + sum_pi4)
    
    if max_size == 4; return gammas; end
    # normalization for size-5 choice probabilities
    sum_pi5       = sum([p^5 for p in item_probs])
    sum_pipj4     = sum_pi4 - sum_pi5
    sum_pi2pj3    = sum_pi2 * sum_pi3 - sum_pi5
    sum_pipjpk3   = (sum_pi3 - sum_pi5 - 2 * sum_pipj4 - sum_pi2pj3) / 2    
    sum_pipj2pk2 = (sum_pi2^2 - sum_pipj4 - 2 * sum_pi2pj3 - sum_pi5) / 2    
    sum_pipjpkpl2 = sum_pi2 * sum_pipjpk - sum_pipjpk3
    sum_pipjpkplpm = (1.0 - sum_pi5 - 5 * sum_pipj4 - 20 * sum_pipjpk3 - 60 * sum_pipjpkpl2 - 30 * sum_pipj2pk2 - 10 * sum_pi2pj3) / 120
    hotset_probs5 = sum([val for (key, val) in H[5]])
    gammas[5] = (1 - hotset_probs5) / (sum_pipjpkplpm + sum_pi5 + sum_pipj4 + sum_pi2pj3 + sum_pipjpk3 + sum_pipjpkpl2 + sum_pipj2pk2)
    
    # TODO: support larger sizes.  We really shouldn't hard code the above.
    # There is a triangular system we can solve to automatically determine these
    # values, but it is annoying to set up.
    if max_size > 5; error("Support only for choices of size <= 5 items"); end

    return gammas
end

function update_hotset_and_model(data::UniversalChoiceDataset, model::UniversalChoiceModel,
                                 choice_to_add::Vector{Int64})
    if in_hotset(choice_to_add, model); error("Choice already in hot set."); end

    choice_to_add_count = 0
    item_counts = zeros(Int64, maximum(data.choices))
    for (size, choice) in iter_choices(data)
        if choice == choice_to_add
            choice_to_add_count += 1
        elseif !in_hotset(choice, model)
            for item in choice
                item_counts[item] += 1
            end
        end
    end

    # Update hot set probability for the new choice
    choices_same_size = countmap(data.sizes)[length(choice_to_add)]
    lc = length(choice_to_add)
    model.H[lc][NTuple{lc, Int64}(choice_to_add)] = choice_to_add_count / choices_same_size
    # Update item counts
    total = sum(item_counts)
    for (i, count) in enumerate(item_counts); model.probs[i] = count / total; end
    # Update normalization parameters
    model.gammas = normalization_values(maximum(data.sizes), model.H, model.probs)
end

function initialize_model(data::UniversalChoiceDataset)
    max_size = maximum(data.sizes)

    # fixed-size probability are just the empirical fraction of selections
    z = zeros(Float64, max_size)
    for (key, value) in countmap(data.sizes); z[key] = value / length(data.sizes); end

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
    H = Vector{Dict{NTuple,Float64}}(max_size)
    for i in 1:max_size; H[i] = Dict{NTuple{i, Int64}, Float64}(); end

    gammas = normalization_values(max_size, H, item_probs)

    return UniversalChoiceModel(z, item_probs, gammas, H)
end
