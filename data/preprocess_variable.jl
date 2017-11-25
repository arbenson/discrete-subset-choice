include("../variable.jl")

using DataStructures

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

function max_size_cutoff(data::VariableChoiceDataset, choice_size_cutoff::Int64, slate_size_cutoff::Int64)
    cutoff_slate_sizes = Int64[]
    cutoff_slates = Int64[]
    cutoff_choice_sizes = Int64[]    
    cutoff_choices = Int64[]
    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data)
        if choice_size <= choice_size_cutoff && slate_size <= slate_size_cutoff
            push!(cutoff_slate_sizes, slate_size)
            append!(cutoff_slates, slate)
            push!(cutoff_choice_sizes, choice_size)
            append!(cutoff_choices, choice)
        end
    end
    return VariableChoiceDataset(cutoff_slate_sizes, cutoff_slates,
                                 cutoff_choice_sizes, cutoff_choices)
end

function item_frequencies(data::VariableChoiceDataset)
    choice_counts = counter(Int64)
    slate_counts = counter(Int64)
    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data)    
        for item in choice
            push!(choice_counts, item)
            push!(slate_counts, item)
        end
    end
    assert(sum(choice_counts) == length(data.choices))
    return choice_counts, slate_counts
end

function item_frequency_cutoff(data::VariableChoiceDataset, selected_cutoff::Int64, total_cutoff::Int64)
    choice_counts, slate_counts = item_frequencies(data)
    cutoff_slate_sizes = Int64[]
    cutoff_slates = Int64[]
    cutoff_choice_sizes = Int64[]    
    cutoff_choices = Int64[]

    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data) 
        keep_choice = true
        for item in slate
            if choice_counts[item] < selected_cutoff || slate_counts[item] < total_cutoff
                keep_choice = false
                break
            end
        end
        if keep_choice
            push!(cutoff_slate_sizes, slate_size)
            append!(cutoff_slates, slate)
            push!(cutoff_choice_sizes, choice_size)
            append!(cutoff_choices, choice)
        end
    end
    return VariableChoiceDataset(cutoff_slate_sizes, cutoff_slates,
                                 cutoff_choice_sizes, cutoff_choices)
end

function apply_cutoffs(dataset::AbstractString, choice_size_cutoff::Int64,
                       slate_size_cutoff::Int64, selected_cutoff::Int64,
                       total_cutoff::Int64)
    @show dataset
    orig_data = read_data(dataset)
    @show length(orig_data.slate_sizes)
    data = max_size_cutoff(orig_data, choice_size_cutoff, slate_size_cutoff)

    while true
        @show length(data.slate_sizes)
        cutoff_data = item_frequency_cutoff(data, selected_cutoff, total_cutoff)
        # Stop once we fail to remove anything
        if length(cutoff_data.slate_sizes) == length(data.slate_sizes); break; end
        data = cutoff_data
    end

    # write out results
    basename = split(dataset, ".")[1]
    output = open("$basename-$(choice_size_cutoff)-$(slate_size_cutoff)-$(selected_cutoff)-$(total_cutoff).txt", "w")

    item_map = Dict{Int64, Int64}()
    function get_id(key::Int64)
        if !haskey(item_map, key)
            n = length(item_map) + 1
            item_map[key] = n
        end
        return item_map[key]
    end
    
    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data)
        slate = [get_id(k) for k in slate]
        sort!(slate)
        choice = [get_id(k) for k in choice]
        sort!(choice)
        write(output, @sprintf("%s;%s\n", join(slate, " "), join(choice, " ")))
    end
end

function main()
    apply_cutoffs("yc-cats.txt", 5, 10, 4, 8)
    apply_cutoffs("yc-items.txt", 5, 10, 4, 8)    
end

main()
