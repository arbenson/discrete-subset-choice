include("universal.jl")
using StatsBase
using Printf

function statistics(dataset_str::AbstractString)
    println("$dataset_str...")
    choice_dataset = read_data(dataset_str)
    sizes = choice_dataset.sizes
    println(@sprintf("\t%d items", length(unique(choice_dataset.choices))))
    println(@sprintf("\t%d choices", length(sizes)))
    size_dist = [(k, v) for (k, v) in countmap(sizes)]
    sort!(size_dist)
    for (k, v) in size_dist
        println(@sprintf("\tz_%d = %0.2f", k, v / length(sizes)))
    end
end

function main()
    statistics("data/bakery-5-25.txt")
    statistics("data/walmart-items-5-25.txt")
    statistics("data/walmart-depts-5-25.txt")
    statistics("data/kosarak-5-25.txt")
    statistics("data/instacart-5-25.txt")
    statistics("data/lastfm-genres-5-25.txt")
end

main()
