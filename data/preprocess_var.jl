include("../variable.jl")

using DataStructures

function max_size_cutoff(data::VariableChoiceDataset, cutoff::Int64)
    cutoff_slate_sizes = Int64[]
    cutoff_slates = Int64[]
    cutoff_choice_sizes = Int64[]    
    cutoff_choices = Int64[]
    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data)
        if choice_size <= cutoff
            push!(cutoff_slate_sizes, slate_size)
            append!(cutoff_slates, slate)
            push!(cutoff_choice_sizes, choice_size)
            append!(cutoff_choices, choice)
        end
    end
    return VariableChoiceDataset(cutoff_slate_sizes, cutoff_slates,
                                 cutoff_choice_sizes, cutoff_choices)
end

function choice_item_frequencies(data::VariableChoiceDataset)
    counts = counter(Int64)
    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data)    
        for item in choice
            push!(counts, item)
        end
    end
    assert(sum(counts) == length(data.choices))
    return counts
end

function item_frequency_cutoff(data::VariableChoiceDataset, freq_cutoff::Int64)
    item_freqs = choice_item_frequencies(data) 
    cutoff_slate_sizes = Int64[]
    cutoff_slates = Int64[]
    cutoff_choice_sizes = Int64[]    
    cutoff_choices = Int64[]

    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data)     
        keep_choice = true
        for item in slate
            if item_freqs[item] < freq_cutoff
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

function apply_cutoffs(dataset::AbstractString, size_cutoff::Int64, freq_cutoff::Int64)
    @show dataset
    orig_data = read_data(dataset)
    @show length(orig_data.slate_sizes)
    data = max_size_cutoff(orig_data, size_cutoff)

    while true
        @show length(data.slate_sizes)
        cutoff_data = item_frequency_cutoff(data, freq_cutoff)
        # Stop once we fail to remove anything
        if length(cutoff_data.slate_sizes) == length(data.slate_sizes); break; end
        data = cutoff_data
    end

    # write out results
    basename = split(dataset, ".")[1]
    output = open("$basename-$size_cutoff-$freq_cutoff.txt", "w")
    for (slate_size, slate, choice_size, choice) in iter_slates_choices(data)
        sort!(slate)
        sort!(choice)
        write(output, @sprintf("%s;%s\n", join(slate, " "), join(choice, " ")))
    end
end

function main()
    apply_cutoffs("yc-items.txt", 5, 5)
    apply_cutoffs("yc-cats.txt", 5, 5)    
end

main()
