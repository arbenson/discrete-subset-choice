include("variable.jl")

using PyPlot

function read_output(basename::AbstractString)
    H_sizes = Int64[]
    model_lls = Float64[]
    open("output/$basename-freq.txt") do f
        for line in eachline(f)
            size, ll = split(line)
            push!(H_sizes, parse(Int64, size))
            push!(model_lls, parse(Float64, ll))
        end
    end
    
    data = read_data("data/$basename.txt")
    n = length(data.choice_sizes)
    num_test = n - convert(Int64, floor(0.8 * n))        
    improvements = exp.((model_lls - model_lls[1]) / num_test)
    return (H_sizes, improvements)
end

function variable_likelihood_gains_plot()
    PyPlot.pygui(true)
    (H_sizes, improvements) = read_output("yc-cats-5-10-4-8")
    plot(H_sizes, improvements, ls="-",  lw=2.5,  color="#1b9e77", label="YcCats")
    (H_sizes, improvements) = read_output("yc-items-5-10-4-8")
    plot(H_sizes, improvements, ls="--", lw=1.75, color="#d95f02", label="YcItems")

    fsz = 18  # font size
    legend(frameon=false, fontsize=fsz-4, loc="lower right")
    ax = axes()
    ax[:tick_params]("both", labelsize=fsz-4, length=8, width=2) 
    xlabel("Number of corrections", fontsize=fsz)
    ylabel("Mean per-choice likelihood gain", fontsize=fsz)
    tight_layout()
    savefig("variable-gains.eps")
    show()
end
