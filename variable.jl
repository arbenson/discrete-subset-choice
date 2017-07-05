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
    f = open(dataset)
    slate_sizes = Int64[]
    slates = Int64[]
    choice_sizes = Int64[]    
    choices = Int64[]
    for line in eachline(f)
        info = split(line, ";")
        slate = [parse(Int64, v) for v in split(info[1])]
        choice = [parse(Int64, v) for v in split(info[2])]
        push!(slate_sizes, length(slate))
        append!(slates, slate)
        push!(choice_sizes, length(choice))
        append!(choices, choice)
    end
    return VariableChoiceDataset(slate_sizes, slates,
                                 choice_sizes, choices)
end

# iterator over choices
function iter_slates_choices(data::VariableChoiceDataset)
    curr_slate_ind = 1
    curr_choice_ind = 1
    slate_vec = Vector{Vector{Int64}}()    
    choice_vec = Vector{Vector{Int64}}()
    for (slate_size, choice_size) in zip(data.slate_sizes, data.choice_sizes)
        slate = data.slates[curr_slate_ind:(curr_slate_ind + slate_size - 1)]
        choice = data.choices[curr_choice_ind:(curr_choice_ind + choice_size - 1)]
        push!(slate_vec, slate)
        push!(choice_vec, choice)
        curr_slate_ind += slate_size
        curr_choice_ind += choice_size
    end
    assert(curr_slate_ind == length(data.slates) + 1)
    assert(curr_choice_ind == length(data.choices) + 1 )          
    return zip(data.slate_sizes, slate_vec, data.choice_sizes, choice_vec)
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

function expsum_util1(model::VariableChoiceModel, slate::Vector{Int64})
    ns = length(slate)
    total = 0.0
    for i = 1:ns; total += exp(model.utilities[i]); end
    return total
end

function expsum_util2(model::VariableChoiceModel, slate::Vector{Int64})
    ns = length(slate)
    total = 0.0
    for i = 1:ns
        si = model.utilities[i]
        for j = i:ns; total += exp(si + model.utilities[j]); end
    end
    return total
end

function expsum_util3(model::VariableChoiceModel, slate::Vector{Int64})
    ns = length(slate)
    total = 0.0
    for i = 1:ns
        si = model.utilities[i]
        for j = i:ns
            sj = si + model.utilities[j]
            for k = j:ns; total += exp(sj + model.utilities[k]); end
        end
    end
    return total
end

function expsum_util4(model::VariableChoiceModel, slate::Vector{Int64})
    ns = length(slate)
    total = 0.0
    for i = 1:ns
        si = model.utilities[i]
        for j = i:ns
            sj = si + model.utilities[j]
            for k = j:ns
                sk = sj + model.utilities[k]
                for l = k:ns; total += exp(sk + model.utilities[l]); end
            end
        end
    end
    return total
end

function expsum_util5(model::VariableChoiceModel, slate::Vector{Int64})
    ns = length(slate)
    total = 0.0
    for i = 1:ns
        si = model.utilities[i]
        for j = i:ns
            sj = si + model.utilities[j]
            for k = j:ns
                sk = sj + model.utilities[k]
                for l = k:ns
                    sl = sk + model.utilities[l]
                    for m = l:ns; total += exp(sl + model.utilties[m]); end
                end
            end
        end
    end
    return total
end

# Given a slate, takes the sum of the exponential of the set utilities for all
# size-k subsets of the slate.  There is one function for each of k = 1, 2, 3, 4, 5.
function expsum_util(model::VariableChoiceModel, slate::Vector{Int64}, size::Int64)
    if     size == 1; return expsum_util1(model, slate)
    elseif size == 2; return expsum_util2(model, slate)
    elseif size == 3; return expsum_util3(model, slate)
    elseif size == 4; return expsum_util4(model, slate)
    elseif size == 5; return expsum_util5(model, slate)
    else error(@sprintf("Cannot handle size %d", size))
    end
end
        
function log_likelihood(model::VariableChoiceModel, data::VariableChoiceDataset)
    ns = length(data.slate_sizes)
    ll = zeros(Float64, ns)
    slate_inds = cumsum(data.slate_sizes) + 1
    unshift!(slate_inds, 1)
    choice_inds = cumsum(data.choice_sizes) + 1
    unshift!(choice_inds, 1)        
    Threads.@threads for i = 1:ns
        slate = data.slates[slate_inds[i]:(slate_inds[i + 1] - 1)]        
        choice = data.choices[choice_inds[i]:(choice_inds[i + 1] - 1)]
        size = length(choice)
        ll[i] += log(model.z[size])        
        for item in choice; ll[i] += model.utilities[item]; end
        ll[i] -= log(expsum_util(slate, size))
    end
    return sum(ll)
end

function add_to_hotset(model::VariableChoiceModel, choice_to_add::Vector{Int64})
    if in_hotset(choice_to_add, model); error("Choice already in hot set."); end
    lc = length(choice_to_add)
    choice_tup = NTuple{lc, Int64}(choice_to_add)
    model.H[lc][choice_tup] = 0.0
end

function initialize_model(data::VariableChoiceDataset)
end
