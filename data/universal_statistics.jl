include("../universal.jl")

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
    statistics("bakery-5-25.txt")
    statistics("walmart-items-5-25.txt")
    statistics("walmart-depts-5-25.txt")
    statistics("kosarak-5-25.txt")
    statistics("instacart-5-25.txt")
    statistics("lastfm-genres-5-25.txt")
end

main()
