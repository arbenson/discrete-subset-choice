include("universal.jl")

using Base.Test

function test()
    function test_normalization(data::UniversalChoiceDataset, model::UniversalChoiceModel)
        p = model.probs
        n = length(model.probs)
        
        p1 = 0.0
        g = model.gammas[1]
        for i = 1:n; p1 += g * p[i]; end
        @test p1 ≈ 1

        p2 = sum([val for (key, val) in model.H if length(key) == 2])
        g = model.gammas[2]
        subset = [0, 0]
        for i = 1:n
            subset[1] = i
            for j = i:n
                subset[2] = j
                if !in_H(model, subset); p2 += g * p[i] * p[j]; end
            end
        end
        @test p2 ≈ 1    

        p3 = sum([val for (key, val) in model.H if length(key) == 3])
        g = model.gammas[3]
        subset = [0, 0, 0]        
        for i = 1:n
            subset[1] = i            
            for j = i:n
                subset[2] = j
                for k = j:n
                    subset[3] = k
                    if !in_H(model, subset); p3 += g * p[i] * p[j] * p[k]; end
                end
            end
        end
        @test p3 ≈ 1    

        p4 = sum([val for (key, val) in model.H if length(key) == 4])
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
                        if !in_H(model, subset); p4 += g * p[i] * p[j] * p[k] * p[l]; end
                    end
                end
            end
        end
        @test p4 ≈ 1    

        p5 = sum([val for (key, val) in model.H if length(key) == 5])        
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
                            if !in_H(model, subset); p5 += g * p[i] * p[j] * p[k] * p[l] * p[m]; end
                        end
                    end
                end
            end
        end
        @test p5 ≈ 1
    end

    function update_H_and_test(data::UniversalChoiceDataset, model::UniversalChoiceModel, choice::Vector{Int64})
        println("adding $choice to H...")
        add_to_H!(model, choice)
        test_normalization(data, model)        
    end

    data = read_data("data/bakery-5-25.txt")
    model = initialize_model(data)
    test_normalization(data, model)
    update_H_and_test(data, model, [23, 40])
    update_H_and_test(data, model, [32, 34])
    update_H_and_test(data, model, [16, 22, 29, 39])
    update_H_and_test(data, model, [40, 46])
    update_H_and_test(data, model, [6, 14, 24])
    update_H_and_test(data, model, [10, 25, 26, 39, 48])
    update_H_and_test(data, model, [16, 39])    

    data = read_data("data/walmart-depts-5-25.txt")
    model = initialize_model(data)
    test_normalization(data, model)
    update_H_and_test(data, model, [10, 32])
    update_H_and_test(data, model, [10, 10, 30])
    update_H_and_test(data, model, [7, 7])
    update_H_and_test(data, model, [10, 32, 32])
    update_H_and_test(data, model, [3, 3, 3, 10])
    update_H_and_test(data, model, [5, 10, 17])
end

test()
