include("variable.jl")

function statistics(dataset_str::AbstractString)
    println("$dataset_str...")
    choice_dataset = read_data(dataset_str)
    slate_sizes = choice_dataset.slate_sizes
    println(@sprintf("\t%d items", length(unique(choice_dataset.slates))))
    println(@sprintf("\t%d choices", length(choice_dataset.choice_sizes)))
    println(@sprintf("\tsizes %d--%d", minimum(slate_sizes), maximum(slate_sizes)))
end

function main()
    statistics("data/yc-cats-5-10-4-8.txt")
    statistics("data/yc-items-5-10-4-8.txt")
end

main()
