include("universal.jl")

using PyPlot

function universal_likelihood_gains_plot(basename::AbstractString, titlename::AbstractString)
    function read_output(filename::AbstractString)
        input_file = open(filename)
        hotset_sizes = []
        model_lls = []
        f = open(filename)
        for line in eachline(f)
            exp_data = split(line)
            push!(hotset_sizes, parse(Int64, exp_data[1]))
            fold_lls = [parse(Float64, v) for v in exp_data[2:end]]
            push!(model_lls, fold_lls)
        end
        close(f)
        # Get mean and std of ll improvements

        data = read_data("data/$basename.txt")
        num_choices = length(data.sizes)

        num_folds = length(model_lls[1])
        num_fold_choices = convert(Int64, floor(num_choices / num_folds)) * ones(Int64, num_folds)
        num_extra = num_choices - sum(num_fold_choices)
        num_fold_choices[end] += num_extra

        mean_improvements = []
        std_improvements = []
        for lls in model_lls
            improvements = exp.((lls - model_lls[1]) ./ num_fold_choices)
            push!(mean_improvements, mean(improvements))
            push!(std_improvements, std(improvements))            
        end
        
        return (hotset_sizes, mean_improvements, std_improvements)
    end

    PyPlot.pygui(true)
    (hotset_sizes, mean_improvements, std_improvements) = read_output("output/$basename-freq.txt")
    plot(hotset_sizes, mean_improvements, label="Frequency")
    (hotset_sizes, mean_improvements, std_improvements) = read_output("output/$basename-lift.txt")
    plot(hotset_sizes, mean_improvements, label="Lift")
    (hotset_sizes, mean_improvements, std_improvements) = read_output("output/$basename-nlift.txt")
    plot(hotset_sizes, mean_improvements, label="Norm. Lift")        

    legend()
    xlabel("Number of corrections")
    ylabel("Mean relative likelihood gain")
    title(titlename)
    show()
    #close()
end

function negative_corrections_plot()
    function read_output(basename::AbstractString)
        f = open("output/$basename-freq-neg-corrections.txt")
        hotset_sizes = Int64[]
        num_neg_corrections = Int64[]
        for line in eachline(f)
            sz, num = [parse(Int64, v) for v in split(line)]
            push!(hotset_sizes, sz)
            push!(num_neg_corrections, num)
        end
        return (hotset_sizes, num_neg_corrections)        
    end

    PyPlot.pygui(true)    

    (hotset_sizes, num_neg_corrections) = read_output("bakery-5-10")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", lw=3, label="Bakery")
    (hotset_sizes, num_neg_corrections) = read_output("walmart-items-5-10")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="--", label="WalmartItems")    
    (hotset_sizes, num_neg_corrections) = read_output("walmart-depts-5-10") 
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls=":", lw=3, label="WalmartDepts")
    (hotset_sizes, num_neg_corrections) = read_output("lastfm-genres-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", lw=1, label="LastfmGenres")
    (hotset_sizes, num_neg_corrections) = read_output("kosarak-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-.", label="Kosarak")
    (hotset_sizes, num_neg_corrections) = read_output("instacart-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", label="Instacart")    

    legend()
    xlabel("Number of corrections")
    ylabel("Fraction negative corrections")
    show()
    #title(titlename)
end

function main()
    universal_likelihood_gains_plot("bakery-5-25", "Bakery")
    universal_likelihood_gains_plot("walmart-items-5-25", "WalmartItems")
    universal_likelihood_gains_plot("walmart-depts-5-25", "WalmartDepts")
    universal_likelihood_gains_plot("lastfm-genres-5-25", "LastfmGenres")
    universal_likelihood_gains_plot("kosarak-5-25", "Kosarak")
    universal_likelihood_gains_plot("instacart-5-25", "Instacart")    
end

#main()
#negative_corrections_plot()
