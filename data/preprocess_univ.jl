include("../universal.jl")

using DataStructures

function max_size_cutoff(data::UniversalChoiceDataset, cutoff::Int64)
    cutoff_choices = Int64[]
    cutoff_sizes = Int64[]
    for (size, choice) in iter_choices(data)
        if size <= cutoff
            append!(cutoff_choices, choice)
            push!(cutoff_sizes, size)
        end
    end
    return UniversalChoiceDataset(cutoff_sizes, cutoff_choices)
end

function item_frequencies(data::UniversalChoiceDataset)
    counts = counter(Int64)
    for (size, choice) in iter_choices(data)
        for item in choice
            push!(counts, item)
        end
    end
    assert(sum(counts) == length(data.choices))
    return counts
end

function item_frequency_cutoff(data::UniversalChoiceDataset, freq_cutoff::Int64)
    item_freqs = item_frequencies(data)    
    cutoff_choices = Int64[]
    cutoff_sizes = Int64[]
    for (size, choice) in iter_choices(data)
        keep_choice = true
        for item in choice
            if item_freqs[item] < freq_cutoff
                keep_choice = false
                break
            end
        end
        if keep_choice
            append!(cutoff_choices, choice)
            push!(cutoff_sizes, size)
        end
    end
    return UniversalChoiceDataset(cutoff_sizes, cutoff_choices)
end

function apply_cutoffs(dataset::AbstractString, size_cutoff::Int64, freq_cutoff::Int64)
    @show dataset
    orig_data = read_data(dataset)
    @show length(orig_data.sizes)
    data = max_size_cutoff(orig_data, size_cutoff)

    while true
        @show length(data.sizes)
        cutoff_data = item_frequency_cutoff(data, freq_cutoff)
        # Stop once we fail to remove anything
        if length(cutoff_data.sizes) == length(data.sizes); break; end
        data = cutoff_data
    end

    # write out results
    basename = split(dataset, ".")[1]

    item_map = Dict{Int64, Int64}()
    function get_id(key::Int64)
        if !haskey(item_map, key)
            n = length(item_map) + 1
            item_map[key] = n
        end
        return item_map[key]
    end

    # random permutation
    out = [(size, choice) for (size, choice) in iter_choices(data)]
    shuffle!(out)

    output = open("$basename-$size_cutoff-$freq_cutoff-clean.txt", "w")
    for (size, choice) in out
        choice = [get_id(c) for c in choice]
        sort!(choice)
        write(output, string(join(choice, ' '), "\n"))
    end
    close(output)
end

function main()
    #apply_cutoffs("bakery.txt", 5, 25)
    apply_cutoffs("walmart-items.txt", 5, 25)
    apply_cutoffs("walmart-depts.txt", 5, 25)
    apply_cutoffs("kosarak.txt", 5, 25)
    #apply_cutoffs("instacart.txt", 5, 25)
    #apply_cutoffs("lastfm-genres.txt", 5, 25)
end

main()
