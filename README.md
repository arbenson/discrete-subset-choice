# discrete-subset-choice

This repository accompanies the paper

- Austin R. Benson, Ravi Kumar, and Andrew Tomkins. A Discrete Choice Model for Subset Selection. Proceedings of WSDM, 2018.


## Universal choice datasets

Compute summary statistics of datasets (Table 1 in paper).
```
julia universal_statistics.jl
```

Reproduce figures.
```
julia> include("universal_figures.jl")
julia> universal_likelihood_gains_plots()  # Figure 1
julia> negative_corrections_plot()  # Figure 2
```

Re-run experiments.
```
julia> negative_corrections_experiments()
```

## Variable choice datasets

Compute summary statistics of datasets (Table 4 in paper):
```
julia variable_statistics.jl
```

Reproduce figures.
```
julia> include("variable_figures.jl")
julia> variable_likelihood_gains_plot()
```

Re-run experiments.
```
julia> include("variable_experiments.jl")
julia> run_variable_model_experiments("yc-cats-5-10-4-8.txt")
julia> run_variable_model_experiments("yc-items-5-10-4-8.txt")
```

