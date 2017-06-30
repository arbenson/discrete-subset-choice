include("universal.jl")

using Base.Test

function test()
    function test_normalization(data::UniversalChoiceDataset, model::UniversalChoiceModel)
        p = model.probs
        n = length(unique(data.choices))
        @show n
        
        p1 = sum([val for (key, val) in model.H[1]])
        g = model.gammas[1]    
        for i = 1:n; p1 += g * p[i]; end
        @test p1 ≈ 1

        p2 = sum([val for (key, val) in model.H[2]])
        g = model.gammas[2]
        subset = [0, 0]
        for i = 1:n
            subset[1] = i
            for j = i:n
                subset[2] = j
                if !in_hotset(subset, model); p2 += g * p[i] * p[j]; end
            end
        end
        @test p2 ≈ 1    
        
        p3 = sum([val for (key, val) in model.H[3]])
        g = model.gammas[3]
        subset = [0, 0, 0]        
        for i = 1:n
            subset[1] = i            
            for j = i:n
                subset[2] = j
                for k = j:n
                    subset[3] = k
                    if !in_hotset(subset, model); p3 += g * p[i] * p[j] * p[k]; end
                end
            end
        end
        @test p3 ≈ 1    
        
        p4 = sum([val for (key, val) in model.H[4]])
        g = model.gammas[4]
        subset = [0, 0, 0, 0]
        for i = 1:n
            subset[1] = i
            for j = i:n
                subset[2] = j
                for k = j:n
                    subset[3] = k
                    for l = k:n
                        subset[4] = l
                        if !in_hotset(subset, model); p4 += g * p[i] * p[j] * p[k] * p[l]; end
                    end
                end
            end
        end
        @test p4 ≈ 1    
        
        p5 = sum([val for (key, val) in model.H[5]])
        g = model.gammas[5]
        subset = [0, 0, 0, 0, 0]
        for i = 1:n
            subset[1] = i
            for j = i:n
                subset[2] = j
                for k = j:n
                    subset[3] = k
                    for l = k:n
                        subset[4] = l
                        for m = l:n
                            subset[5] = m
                            if !in_hotset(subset, model); p5 += g * p[i] * p[j] * p[k] * p[l] * p[m]; end
                        end
                    end
                end
            end
        end
        @test p5 ≈ 1
    end

    function update_hotset_and_test(data::UniversalChoiceDataset, model::UniversalChoiceModel, choice::Vector{Int64})
        update_hotset_and_model(data, model, choice)
        test_normalization(data, model)        
    end

    data = read_data("data/bakery-5-10.txt")
    model = initialize_model(data)
    test_normalization(data, model)
    update_hotset_and_test(data, model, [1, 2])
    update_hotset_and_test(data, model, [5, 26])
    update_hotset_and_test(data, model, [8, 13, 14])
    update_hotset_and_test(data, model, [18, 19, 20, 21])
    update_hotset_and_test(data, model, [14, 20, 23, 24, 25])
    
    data = read_data("data/walmart-depts-5-10.txt")
    model = initialize_model(data)
    test_normalization(data, model)
    update_hotset_and_test(data, model, [17, 24, 24, 32])
    update_hotset_and_test(data, model, [13, 31])
    update_hotset_and_test(data, model, [5, 9, 12, 12])
    update_hotset_and_test(data, model, [5, 9, 23])
    update_hotset_and_test(data, model, [5, 5, 5, 5, 5])
    update_hotset_and_test(data, model, [36, 36])
    update_hotset_and_test(data, model, [5, 9, 20, 30])
end

test()
