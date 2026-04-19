import Pkg

# Install all required Julia packages for the MSTP engine.
# Run once after installing Julia: julia setup.jl
pkgs = [
    "BlockDiagonals",
    "DataStructures",
    "Distributions",
    "HiGHS",
    "JuMP",
    "LinearAlgebra",
    "Random",
    "SDDP",
    "Statistics",
]

Pkg.add(pkgs)
Pkg.precompile()
println("MSTP Julia engine setup complete.")
