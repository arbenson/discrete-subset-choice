include("variable.jl")

function variable_improvements(data::VariableChoiceDataset, num_updates::Int64,
                               basename::AbstractString, update_type::AbstractString)
    # Vector for randomly splitting into training / test data
    n = length(data.choice_sizes)
    inds = collect(1:n)
    shuffle!(inds)
    num_folds = 10
    log_likelihoods = zeros(Float64, num_updates + 1, num_folds)
    
    for fold in 1:num_folds
        println(@sprintf("fold %d of %d", fold, num_folds))

        # Split into training / test data
        training = ones(Int64, n)
        fold_size = convert(Int64, floor(n / num_folds))
        test_start_ind = (fold - 1) * fold_size + 1
        test_end_ind = fold * fold_size
        if fold == num_folds; test_end_ind = n; end
        training[test_start_ind:test_end_ind] = 0
        training = training[inds]
        training_slate_sizes = Int64[]
        training_slates = Int64[]        
        training_choice_sizes = Int64[]
        training_choices = Int64[]
        test_slate_sizes = Int64[]
        test_slates = Int64[]        
        test_choice_sizes = Int64[]
        test_choices = Int64[]
        for (ind, (slate_size, choice_size)) in zip(data.slate_sizes, data.choice_sizes)
            slate = data.slates[curr_slate_ind:(curr_slate_ind + slate_size - 1)]
            choice = data.choices[curr_choice_ind:(curr_choice_ind + choice_size - 1)]
            if training[ind] == 1
                push!(training_slate_sizes, slate_size)
                append!(training_slates, slate)
                push!(training_choice_sizes, choice_size)
                append!(training_choices, choice)
            else
                push!(test_slate_sizes, slate_size)
                append!(test_slates, slate)
                push!(test_choice_sizes, choice_size)
                append!(test_choices, choice)
            end
        end
        training_data =
            VariableChoiceDataset(training_slate_sizes, training_slates,
                                  training_choice_sizes, training_choices)
        training_data =
            VariableChoiceDataset(test_slate_sizes, test_slates,
                                  test_choice_sizes, test_choices)
        log_likelihoods[1, fold] = log_likelihood(model, test_data)

        # update here
        item_counts = zeros(Int64, maximum(data.choices))
        for (size, choice) in iter_choices(training_data)
            for item in choice; item_counts[item] += 1; end
        end
        counts = [(count, choice_tup) for (choice_tup, count) in get_subset_counts(training_data)]
        sort!(counts, rev=true)

        choices_to_add = [collect(choice_tup) for (count, choice_tup) in counts[1:num_updates]]
        for (i, choice) in enumerate(choices_to_add)
            println(@sprintf("iteration %d of %d", i, num_updates))
            add_to_hotset!(model, choice)
            log_likelihoods[i + 1, fold] = log_likelihood(model, test_data)
        end
    end

    output = open("output/$basename-freq.txt", "w")
    for i = 0:num_updates
        write(output, @sprintf("%d %s\n", i, join(log_likelihoods[i + 1, :], " ")))
    end
end


function universal_improvement_experiments()
    function run_universal_improvement_experiment(dataset_file::AbstractString)
        data = read_data(dataset_file)
        basename = split(split(dataset_file, "/")[end], ".")[1]
        num_items = length(unique(data.choices))
        num_updates = min(num_items, 1000)
        #universal_improvements(data, num_updates, basename, "f")
        #universal_improvements(data, num_updates, basename, "nl")
        #universal_improvements(data, num_updates, basename, "l")
        #universal_improvements(data, num_updates, basename, "lev")
        #universal_improvements(data, num_updates, basename, "g")
        #negative_corrections(data, num_updates, basename)
        biggest_corrections(data, num_updates, basename)
    end

    #run_universal_improvement_experiment("data/bakery-5-25.txt")
    #run_universal_improvement_experiment("data/walmart-depts-5-25.txt")
    #run_universal_improvement_experiment("data/walmart-items-5-25.txt")
    #run_universal_improvement_experiment("data/lastfm-genres-5-25.txt")
    #run_universal_improvement_experiment("data/kosarak-5-25.txt")
    run_universal_improvement_experiment("data/instacart-5-25.txt")
end

universal_improvement_experiments()
