include("universal.jl")

using PyPlot

function read_output(basename::AbstractString, filename::AbstractString)
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


function universal_likelihood_gains_plot(basename::AbstractString, titlename::AbstractString)
    PyPlot.pygui(true)
    (hotset_sizes, means, stds) = read_output(basename, "output/$basename-freq.txt")
    plot(hotset_sizes, means, ls="-", lw=2, label="Frequency", color="#e41a1c")
    fill_between(hotset_sizes, means + stds, means - stds, color="#fbb4ae")
    (hotset_sizes, means, stds) = read_output(basename, "output/$basename-lift.txt")
    plot(hotset_sizes, means, ls="--", lw=2, label="Lift", color="#377eb8")
    fill_between(hotset_sizes, means + stds, means - stds, color="#b3cde3")    
    (hotset_sizes, means, stds) = read_output(basename, "output/$basename-nlift.txt")
    plot(hotset_sizes, means, ls=":", lw=3, label="Norm. Lift", color="#4daf4a")
    fill_between(hotset_sizes, means + stds, means - stds, color="#ccebc5")

    fsz = 20  # font size
    if titlename == "WalmartDepts"
        legend(frameon=false, fontsize=fsz-2)
    end
    ax = axes()
    ax[:tick_params]("both", labelsize=fsz-4, length=8, width=2) 
    xlabel("Number of corrections", fontsize=fsz)
    ylabel("Mean relative likelihood gain", fontsize=fsz)
    title(titlename, fontsize=fsz)
    tight_layout()
    savefig("universal-gains-$basename.eps")
    close()
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

    (hotset_sizes, num_neg_corrections) = read_output("bakery-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="--", lw=7, label="Bakery")
    (hotset_sizes, num_neg_corrections) = read_output("walmart-items-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls=":", lw=5, label="WalmartItems")    
    (hotset_sizes, num_neg_corrections) = read_output("walmart-depts-5-25") 
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls=":", lw=3, label="WalmartDepts")
    (hotset_sizes, num_neg_corrections) = read_output("lastfm-genres-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", lw=2.5, label="LastfmGenres")
    (hotset_sizes, num_neg_corrections) = read_output("kosarak-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", lw=1.75, label="Kosarak")
    (hotset_sizes, num_neg_corrections) = read_output("instacart-5-25")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", lw=1, label="Instacart")    

    fsz = 18  # font size
    legend(frameon=false, fontsize=fsz-4, loc="upper right")    
    xlabel("Number of corrections", fontsize=fsz)
    ylabel("Fraction negative corrections", fontsize=fsz)
    ax = axes()    
    ax[:tick_params]("both", labelsize=fsz-4, length=8, width=2)
    ax[:tick_params]("both", which="minor", length=4, width=1)
    tight_layout()
    savefig("negative_corrections.eps")
    close()
end

function main()
    universal_likelihood_gains_plot("bakery-5-25", "Bakery")
    universal_likelihood_gains_plot("walmart-depts-5-25", "WalmartDepts")
    universal_likelihood_gains_plot("walmart-items-5-25", "WalmartItems")
    universal_likelihood_gains_plot("lastfm-genres-5-25", "LastfmGenres")
    universal_likelihood_gains_plot("kosarak-5-25", "Kosarak")
    universal_likelihood_gains_plot("instacart-5-25", "Instacart")
end
main()
#negative_corrections_plot()

