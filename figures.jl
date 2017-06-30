include("universal.jl")

using PyPlot

function read_output(filename::AbstractString)
    input_file = open(filename)
    hotset_sizes = []
    model_lls = []
    f = open(filename)
    for line in eachline(f)
        exp_data = split(line)
        push!(hotset_sizes, parse(Int64, exp_data[1]))
        push!(model_lls, parse(Float64, exp_data[2]))
    end
    close(f)
    return (hotset_sizes, model_lls)
end

function make_plot(basename::AbstractString, titlename::AbstractString)
    data = read_data("data/$basename.txt")
    num_choices = length(data.sizes)

    PyPlot.pygui(true)
    
    (hotset_sizes, model_lls) = read_output("output/$basename-freq.txt")
    ll0 = model_lls[1]
    freq_improvements = [exp((ll - ll0) / num_choices) for ll in model_lls]
    plot(hotset_sizes, freq_improvements, label="Frequency")

    (hotset_sizes, model_lls) = read_output("output/$basename-lift.txt")
    ll0 = model_lls[1]
    lift_improvements = [exp((ll - ll0) / num_choices) for ll in model_lls]
    plot(hotset_sizes, lift_improvements, label="Lift")    

    legend()
    xlabel("Number of corrections")
    ylabel("Mean relative likelihood gain")
    title(titlename)
    show()
    #close()
end

function main()
    make_plot("bakery-5-10", "Bakery")
    make_plot("walmart-items-5-10", "WalmartItems")
    make_plot("walmart-depts-5-10", "WalmartDepts")
    make_plot("kosarak-5-25", "Kosarak")
    make_plot("lastfm-genres-5-25", "LastfmGenres")
    make_plot("instacart-5-25", "Instacart")    
end

#main()
