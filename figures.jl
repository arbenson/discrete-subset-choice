include("universal.jl")

using PyPlot

function make_plot(basename::AbstractString, titlename::AbstractString)
    freq_output = open("output/$basename-freq.txt")
    hotset_sizes = []
    model_lls = []
    for line in eachline(freq_output)
        exp_data = split(line)
        push!(hotset_sizes, parse(Int64, exp_data[1]))
        push!(model_lls, parse(Float64, exp_data[2]))
    end

    data = read_data("data/$basename.txt")
    num_choices = length(data.sizes)
    
    ll0 = model_lls[1]
    improvements = [exp((ll - ll0) / num_choices) for ll in model_lls]
    
    PyPlot.pygui(true)
    plot(hotset_sizes, improvements)
    xlabel("Number of corrections")
    ylabel("Mean relative likelihood gain")
    title(titlename)
    show()
    #close()
end

function main()
    #make_plot("bakery-5-10", "Bakery")
    #make_plot("walmart-items-5-10", "WalmartItems")
    make_plot("walmart-depts-5-10", "WalmartDepts")
    #make_plot("kosarak-5-25", "Kosarak")
end

main()
