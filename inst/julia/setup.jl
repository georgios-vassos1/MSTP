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

# Gurobi requires a valid licence; install but do not error if unavailable
try
    Pkg.add("Gurobi")
catch e
    @warn "Gurobi could not be installed — falling back to HiGHS will still work." e
end

Pkg.precompile()
println("MSTP Julia engine setup complete.")
