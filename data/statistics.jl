using StatsBase

function read_data(dataset::AbstractString)
    f = open(dataset)
    choices = Int64[]
    sizes = Int64[]
    for line in eachline(f)
        choice = [parse(Int64, v) for v in split(line)]
        append!(choices, choice)
        push!(sizes, length(choice))
    end
    return (choices, sizes)
end

function statistics(dataset::AbstractString)
    @show dataset
    choices, sizes = read_data(dataset)
    println(@sprintf("%d items", length(unique(choices))))
    println(@sprintf("%d choices", length(sizes)))
    size_dist = [(k, v) for (k, v) in countmap(sizes)]
    sort!(size_dist)
    for (k, v) in size_dist
        println(@sprintf("z_%d = %0.2f", k, v / length(sizes)))
    end
end

function main()
    statistics("bakery-5-25-clean.txt")
    statistics("walmart-items-5-25-clean.txt")
    statistics("walmart-depts-5-25-clean.txt")
    statistics("kosarak-5-25-clean.txt")
    statistics("lastfm-genres-5-25-clean.txt")
    statistics("instacart-5-25-clean.txt")
end

main()
