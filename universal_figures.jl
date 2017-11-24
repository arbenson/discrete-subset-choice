include("universal.jl")

using PyPlot

function read_output(basename::AbstractString, filename::AbstractString)
    hotset_sizes = Int64[]
    model_lls = Float64[]
    f = open(filename)
    for line in eachline(f)
        exp_data = split(line)
        push!(hotset_sizes, parse(Int64, exp_data[1]))
        push!(model_lls, parse(Float64, exp_data[2]))
    end
    close(f)

    data = read_data("data/$basename.txt")
    n = length(data.sizes)
    ntest = n - convert(Int64, floor(0.8 * n))
    
    return (hotset_sizes, model_lls, ntest)
end

function universal_likelihood_gains_plot(basename::AbstractString, titlename::AbstractString)
    (freq_sizes,  freq_lls,  freq_ntest)  = read_output(basename, "output/$basename-single-freq.txt")
    (lift_sizes,  lift_lls,  lift_ntest)  = read_output(basename, "output/$basename-single-lift.txt")
    (nlift_sizes, nlift_lls, nlift_ntest) = read_output(basename, "output/$basename-single-nlift.txt")

    base_ll = freq_lls[1]
    ntest = freq_ntest
    norm_ll(lls::Vector{Float64}) = exp.((lls - base_ll) / ntest)

    PyPlot.pygui(true)
    plot(freq_sizes,  norm_ll(freq_lls),  ls="-",  lw=4.5, label="Frequency")
    plot(lift_sizes,  norm_ll(lift_lls),  ls="--", lw=4.5, label="Lift")
    plot(nlift_sizes, norm_ll(nlift_lls), ls=":",  lw=4.5, label="Norm. Lift")

    dpp_file = "output/$basename-dpp.txt"
    if isfile(dpp_file)
        f = open(dpp_file)
        dpp_ll = parse(Float64, readline(f))
        close(f)
        plot(freq_sizes, norm_ll(dpp_ll * ones(length(freq_sizes))), ls="-", lw=1.5, label="DPP")
    end

    fsz = 18  # font size
    if titlename == "Bakery"
        legend(frameon=false, fontsize=fsz, loc=(0.57, 0.42))
    end
    ax = axes()
    ax[:tick_params]("both", labelsize=fsz-4, length=8, width=2) 
    xlabel("Number of corrections", fontsize=fsz)
    ylabel("Mean per-choice likelihood gain", fontsize=fsz)
    title(titlename, fontsize=fsz)
    tight_layout()
    #show()
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

    (hotset_sizes, num_neg_corrections) = read_output("bakery-5-25-clean")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="--", lw=7, label="Bakery")
    (hotset_sizes, num_neg_corrections) = read_output("walmart-items-5-25-clean")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls=":", lw=5, label="WalmartItems")    
    (hotset_sizes, num_neg_corrections) = read_output("walmart-depts-5-25-clean") 
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls=":", lw=3, label="WalmartDepts")
    (hotset_sizes, num_neg_corrections) = read_output("lastfm-genres-5-25-clean")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", lw=2.5, label="LastfmGenres")
    (hotset_sizes, num_neg_corrections) = read_output("kosarak-5-25-clean")
    semilogx(hotset_sizes, num_neg_corrections ./ hotset_sizes, ls="-", lw=1.75, label="Kosarak")
    (hotset_sizes, num_neg_corrections) = read_output("instacart-5-25-clean")
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
    universal_likelihood_gains_plot("bakery-5-25-clean", "Bakery")
    universal_likelihood_gains_plot("walmart-depts-5-25-clean", "WalmartDepts")
    universal_likelihood_gains_plot("walmart-items-5-25-clean", "WalmartItems")
    universal_likelihood_gains_plot("lastfm-genres-5-25-clean", "LastfmGenres")
    universal_likelihood_gains_plot("instacart-5-25-clean", "Instacart")
    universal_likelihood_gains_plot("kosarak-5-25-clean", "Kosarak")
end
main()
#negative_corrections_plot()
