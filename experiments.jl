include("universal.jl")

function subset_counts(data::UniversalChoiceDataset)
    # Get the counts
    counts = Dict{NTuple, Int64}()
    for (size, choice) in iter_choices(data)
        choice_tup = NTuple{length(choice), Int64}(choice)
        if length(choice) > 1
            if !haskey(counts, choice_tup); counts[choice_tup] = 0; end
            counts[choice_tup] += 1
        end
    end
    return counts
end

function update_by_frequency(data::UniversalChoiceDataset, num_updates::Int64,
                             output_filename::AbstractString)
    counts = [(count, choice_tup) for (choice_tup, count) in subset_counts(data)]
    sort!(counts, rev=true)
    choices_to_add = [collect(choice_tup) for (count, choice_tup) in counts[1:num_updates]]

    model = initialize_model(data)
    log_likelihoods = Float64[]
    push!(log_likelihoods, log_likelihood(data, model))

    for (i, choice) in enumerate(choices_to_add)
        println(@sprintf("iteration %d of %d", i, num_updates))
        update_hotset_and_model(data, model, choice)
        ll = log_likelihood(data, model)
        push!(log_likelihoods, ll)
    end

    output = open(output_filename, "w")
    for (i, ll) in enumerate(log_likelihoods)
        write(output, @sprintf("%d %f\n", i - 1, ll))
    end
end

function update_by_lift(data::UniversalChoiceDataset, num_updates::Int64,
                        output_filename::AbstractString)
    # individual item counts
    item_counts = zeros(Int64, maximum(data.choices))
    for (size, choice) in iter_choices(data)
        for item in choice; item_counts[item] += 1; end
    end
    lifts = Vector{Tuple{Float64,NTuple}}()
    for (choice_tup, subset_count) in subset_counts(data)
        subset_item_counts = [item_counts[item] for item in choice_tup]
        push!(lifts, (subset_count^2 / prod(subset_item_counts), choice_tup))
    end
    sort!(lifts, rev=true)
    choices_to_add = [collect(choice_tup) for (_, choice_tup) in lifts[1:num_updates]]

    model = initialize_model(data)    
    log_likelihoods = Float64[]
    push!(log_likelihoods, log_likelihood(data, model))
    
    for (i, choice) in enumerate(choices_to_add)
        println(@sprintf("iteration %d of %d", i, num_updates))
        update_hotset_and_model(data, model, choice)
        ll = log_likelihood(data, model)
        push!(log_likelihoods, ll)
    end

    output = open(output_filename, "w")
    for (i, ll) in enumerate(log_likelihoods)
        write(output, @sprintf("%d %f\n", i - 1, ll))
    end
end


function frequency_experiments()
    function run_experiment(dataset_file::AbstractString)
        data = read_data(dataset_file)
        basename = split(split(dataset_file, "/")[end], ".")[1]
        num_items = length(unique(data.choices))
        num_updates = min(num_items, 1000)
        update_by_frequency(data, num_updates, "output/$basename-freq.txt")
        #update_by_lift(data, num_updates, "output/$basename-lift.txt")              
    end

    run_experiment("data/bakery-5-10.txt")
    #run_experiment("data/walmart-depts-5-10.txt")
    #run_experiment("data/walmart-items-5-10.txt")
    #run_experiment("data/kosarak-5-25.txt")
    #run_experiment("data/lastfm-genres-5-25.txt")
    #run_experiment("data/instacart-5-25.txt")
end

frequency_experiments()
