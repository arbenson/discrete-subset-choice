include("variable.jl")

using PyPlot

function read_output(basename::AbstractString)
    hotset_sizes = []
    model_lls = []
    f = open("output/$basename-freq.txt")
    for line in eachline(f)
        exp_data = split(line)
        push!(hotset_sizes, parse(Int64, exp_data[1]))
        fold_lls = [parse(Float64, v) for v in exp_data[2:end]]
        push!(model_lls, fold_lls)
    end
    close(f)
    
    data = read_data("data/$basename.txt")
    num_choices = length(data.choice_sizes)
    
    num_folds = length(model_lls[1])
    num_fold_choices = convert(Int64, floor(num_choices / num_folds)) * ones(Int64, num_folds)
    num_extra = num_choices - sum(num_fold_choices)
    num_fold_choices[end] += num_extra
    
    mean_improvements = []
    std_improvements = []
    for lls in model_lls
        improvements = exp.((lls - model_lls[1]) ./ num_fold_choices)
        #improvements = lls ./ model_lls[1]
        push!(mean_improvements, mean(improvements))
        push!(std_improvements, std(improvements))            
    end
    
    return (hotset_sizes, mean_improvements, std_improvements)
end

function variable_likelihood_gains_plot()
    PyPlot.pygui(true)
    #(hotset_sizes, means, stds) = read_output("yc-cats-5-8-5")
    #(hotset_sizes, means, stds) = read_output("yc-cats-3-10-4-8")
    #(hotset_sizes, means, stds) = read_output("yc-cats-5-15-4-8")
    #(hotset_sizes, means, stds) = read_output("yc-cats-5-5")
    (hotset_sizes, means, stds) = read_output("yc-cats-5-10-4-8")

    @show means
    plot(hotset_sizes, means, ls="-", lw=2.5, label="YcCats")

    fsz = 20  # font size
    legend(frameon=false, fontsize=fsz-2)
    ax = axes()
    ax[:tick_params]("both", labelsize=fsz-4, length=8, width=2) 
    xlabel("Number of corrections", fontsize=fsz)
    ylabel("Mean relative likelihood gain", fontsize=fsz)
    tight_layout()
    savefig("variable-gains.eps")
    close()
end

variable_likelihood_gains_plot()
