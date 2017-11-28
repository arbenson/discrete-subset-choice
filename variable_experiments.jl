include("variable.jl")

function variable_model_freq_improvements(data::VariableChoiceDataset,
                                          num_updates::Int64,
                                          basename::AbstractString)
    # Vector for randomly splitting into training / test data
    n = length(data.choice_sizes)
    training_end = convert(Int64, floor(0.8 * n))    
    log_likelihoods = zeros(Float64, num_updates + 1)    

    # Split into training / test data
    training_slate_sizes  = Int64[]
    training_slates       = Int64[]        
    training_choice_sizes = Int64[]
    training_choices      = Int64[]
    test_slate_sizes      = Int64[]
    test_slates           = Int64[]        
    test_choice_sizes     = Int64[]
    test_choices          = Int64[]
    
    slate_inds = index_points(data.slate_sizes)
    choice_inds = index_points(data.choice_sizes)
    for ind = 1:length(data.choice_sizes)
        slate = data.slates[slate_inds[ind]:(slate_inds[ind + 1] - 1)]        
        choice = data.choices[choice_inds[ind]:(choice_inds[ind + 1] - 1)]
        if ind <= training_end
            push!(training_slate_sizes, length(slate))
            append!(training_slates, slate)
            push!(training_choice_sizes, length(choice))
            append!(training_choices, choice)
        else
            push!(test_slate_sizes, length(slate))
            append!(test_slates, slate)
            push!(test_choice_sizes, length(choice))
            append!(test_choices, choice)
        end
    end
    training_data = VariableChoiceDataset(training_slate_sizes, training_slates,
                                          training_choice_sizes, training_choices)
    test_data = VariableChoiceDataset(test_slate_sizes, test_slates,
                                      test_choice_sizes, test_choices)

    model = initialize_model(training_data)
    learn_model!(model, training_data)
    log_likelihoods[1] = log_likelihood(model, test_data)
    @show model.z
    
    # Get subsets by frequency
    choice_inds = index_points(training_data.choice_sizes)
    subset_counts = Dict{NTuple, Int64}()
    for i = 1:length(training_data.choice_sizes)
        choice = training_data.choices[choice_inds[i]:(choice_inds[i + 1] - 1)]
        choice_tup = vec2ntuple(choice)
        if length(choice_tup) > 1
            if !haskey(subset_counts, choice_tup); subset_counts[choice_tup] = 0; end
            subset_counts[choice_tup] += 1
        end
    end
    counts = [(count, choice_tup) for (choice_tup, count) in subset_counts]
    sort!(counts, rev=true)
    choices_to_add = [collect(choice_tup) for (count, choice_tup) in counts[1:num_updates]]

    for (i, choice) in enumerate(choices_to_add)
        println(@sprintf("iteration %d of %d", i, num_updates))
        add_to_H!(model, choice)
        learn_utilities!(model, training_data)
        log_likelihoods[i + 1] = log_likelihood(model, test_data)
        ntest = length(test_data.choice_sizes)
        ll_new  = log_likelihoods[i + 1]
        ll_orig = log_likelihoods[1]
        improve = exp((ll_new - ll_orig) / ntest)
        @show ll_new, ll_orig, ntest, improve
    end

    output = open("output/$basename-freq.txt", "w")
    for i = 0:num_updates
        write(output, @sprintf("%d %s\n", i, join(log_likelihoods[i + 1, :], " ")))
    end
end

function variable_likelihood_experiment(dataset_file::AbstractString)
    data = read_data("data/$dataset_file")
    basename = split(split(dataset_file, "/")[end], ".")[1]
    num_items = length(unique(data.choices))
    variable_model_freq_improvements(data, 20, basename)
end
