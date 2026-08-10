# MSTP Julia engine — single entry point loaded by R via JuliaCall.
# Sources all modules in dependency order.

const MSTP_JULIA_DIR = @__DIR__

# Activate the pinned engine environment before any `using`, so the modules
# resolve SDDP/JuMP/HiGHS against Manifest.toml rather than the global env.
import Pkg
if Base.active_project() != joinpath(MSTP_JULIA_DIR, "Project.toml")
    Pkg.activate(MSTP_JULIA_DIR; io = devnull)
end

include(joinpath(MSTP_JULIA_DIR, "utils.jl"))
include(joinpath(MSTP_JULIA_DIR, "types.jl"))
include(joinpath(MSTP_JULIA_DIR, "model.jl"))
include(joinpath(MSTP_JULIA_DIR, "api.jl"))
