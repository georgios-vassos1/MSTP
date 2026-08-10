import Pkg

# Instantiate the pinned MSTP Julia environment (Project.toml + Manifest.toml in
# this directory). Reproducibility relies on this lock, not on whatever versions
# happen to be in the user's global environment. Run once after installing
# Julia: `julia inst/julia/setup.jl`.
Pkg.activate(@__DIR__)
Pkg.instantiate()   # resolves against Manifest.toml; downloads only if missing
Pkg.precompile()
println("MSTP Julia engine setup complete (pinned env: ", @__DIR__, ").")
