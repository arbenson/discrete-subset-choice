include("universal.jl")

# Get the num most frequent choices.
function most_freq_choice_tups(data::UniversalChoiceDataset, num::Int64)
    counts = [(count, tup) for (tup, count) in get_subset_counts(data) if length(tup) > 1]
    sort!(counts, rev=true)
    return [collect(tup) for (count, tup) in counts[1:num]]
end

maximum_H_size(data::UniversalChoiceDataset) = min(length(unique(data.choices)), 1000)

function universal_improvements(data::UniversalChoiceDataset, num_updates::Int64,
                                basename::AbstractString, update_type::AbstractString,
                                timing::Bool=false)
    # Vector for randomly splitting into training / test data
    n = length(data.sizes)
    log_likelihoods = zeros(Float64, num_updates + 1)

    # Split into training / test data
    training_end = convert(Int64, floor(0.8 * n))
    training_sizes   = Int64[]
    training_choices = Int64[]
    test_sizes       = Int64[]
    test_choices     = Int64[]
    for (ind, (size, choice)) in enumerate(iter_choices(data))
        if ind <= training_end
            push!(training_sizes, size)
            append!(training_choices, choice)
        else
            push!(test_sizes, size)
            append!(test_choices, choice)
        end
    end
    training_data = UniversalChoiceDataset(training_sizes, training_choices)
    test_data = UniversalChoiceDataset(test_sizes, test_choices)
    model = initialize_model(training_data)
    log_likelihoods[1] = log_likelihood(model, test_data)
    
    item_counts = zeros(Int64, maximum(data.choices))
    for (size, choice) in iter_choices(training_data)
        for item in choice; item_counts[item] += 1; end
    end

    choices_to_add = most_freq_choice_tups(training_data, num_updates)
    if     update_type == "f"
        # Keep frequency-based updates
    elseif update_type == "nl"
        # normalized lift-based updates
        lifts = Vector{Tuple{Float64,NTuple}}()
        for (choice_tup, subset_count) in get_subset_counts(training_data)
            if length(choice_tup) > 1
                subset_item_counts = [item_counts[item] for item in choice_tup]
                push!(lifts, (subset_count^2 / prod(subset_item_counts), choice_tup))
            end
        end
        sort!(lifts, rev=true)
        choices_to_add = [collect(choice_tup) for (_, choice_tup) in lifts[1:num_updates]]
    elseif update_type == "l"
        # Lift-based updates
        lifts = Vector{Tuple{Float64,NTuple}}()
        for (choice_tup, subset_count) in get_subset_counts(training_data)
            if length(choice_tup) > 1
                subset_item_counts = [item_counts[item] for item in choice_tup]
                push!(lifts, (subset_count / prod(subset_item_counts), choice_tup))
            end
        end
        sort!(lifts, rev=true)
        choices_to_add = [collect(choice_tup) for (_, choice_tup) in lifts[1:num_updates]]
    else
        error("Unknown update type")
    end

    for (i, choice) in enumerate(choices_to_add)
        add_to_hotset!(model, choice)
        if !timing;
            println(@sprintf("iteration %d of %d", i, num_updates))
            log_likelihoods[i + 1] = log_likelihood(model, test_data)
        end
    end

    # If we are doing a timing experiment, we are done.
    if timing
        log_likelihoods[end] = log_likelihood(model, test_data)
        return
    end

    if     update_type == "f";  output = open("output/$basename-freq.txt", "w")
    elseif update_type == "l";  output = open("output/$basename-lift.txt", "w")
    elseif update_type == "nl"; output = open("output/$basename-nlift.txt", "w")
    else   error("Unknown update type")    end
    for i = 0:num_updates
        ll = log_likelihoods[i + 1]
        write(output, "$i $ll\n")
    end
    close(output)
end

function negative_corrections_experiment(basename::AbstractString)
    println("$basename...")
    data = read_data("data/$basename.txt")
    num_updates = maximum_H_size(data)
    choices_to_add = most_freq_choice_tups(data, num_updates)
    model = initialize_model(data)
    num_negative_corrections = Int64[]
    for (i, choice) in enumerate(choices_to_add)
        add_to_hotset!(model, choice)
        # Get negative corrections
        count = 0
        for (subset, val) in model.H
            separable_prob = prod([model.probs[item] for item in subset])
            gamma = model.gammas[length(subset)]
            correction = val - gamma * separable_prob
            if correction < 0; count += 1; end
        end
        println("iteration $i of $num_updates: $count corretions < 0")
        push!(num_negative_corrections, count)
    end
    
    open("output/$basename-freq-neg-corrections.txt", "w") do output
        for (i, num) in enumerate(num_negative_corrections)
            write(output, @sprintf("%d %d\n", i, num))
        end
    end
end

function biggest_corrections_experiment(basename::AbstractString)
    data = read_data("data/$basename.txt")
    num_updates = maximum_H_size(data)    
    choices_to_add = most_freq_choice_tups(data, num_updates)
    model = initialize_model(data)
    num_negative_corrections = Int64[]
    for choice in choices_to_add; add_to_hotset!(model, choice); end

    results = []
    for (subset, val) in model.H
        separable_prob = prod([model.probs[item] for item in subset])
        gamma = model.gammas[length(subset)]
        correction = val - gamma * separable_prob
        push!(results, (correction, subset))
    end

    sort!(results)
    open("output/$basename-biggest-corrections.txt", "w") do output
        for i = 1:5
            line = @sprintf("%d:\t%s\t%f\n", i, join(results[i][2], " "), results[i][1])
            println(line)
            write(output, line)
        end
        for i = 1:5
            ind = length(results) - i + 1
            line = @sprintf("%d:\t%s\t%f\n", i, join(results[ind][2], " "), results[ind][1])
            println(line)
            write(output, line)
        end
    end
end

function timing_experiment(basename::AbstractString)
    data = read_data("data/$basename.txt")
    num_updates = maximum_H_size(data)    
    # warmup
    universal_improvements(data, num_updates, basename, "f", true)
    # run
    tic()
    @time universal_improvements(data, num_updates, basename, "f", true)
    tf = toc()
    open("output/$basename-times.txt", "w") do f
        write(f, @sprintf("freq: %f\n", tf))
    end
end

# Use frequency, lift, and normalized lift heuristics to generate models.
function universal_likelihood_experiments(basename::AbstractString)
    println("$basename...")
    data = read_data("data/$basename.txt")
    num_updates = min(length(unique(data.choices)), 1000)
    println("frequency heuristic...")
    universal_improvements(data, num_updates, basename, "f")
    println("lift heuristic...")    
    universal_improvements(data, num_updates, basename, "l")
    println("normalized lift heuristic...")        
    universal_improvements(data, num_updates, basename, "nl")
end


# timing_experiment("kosarak-5-25")
# biggest_corrections_experiment("lastfm-genres-5-25")

