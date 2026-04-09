# MSTP Julia engine — single entry point loaded by R via JuliaCall.
# Sources all modules in dependency order.

const MSTP_JULIA_DIR = @__DIR__

include(joinpath(MSTP_JULIA_DIR, "utils.jl"))
include(joinpath(MSTP_JULIA_DIR, "types.jl"))
include(joinpath(MSTP_JULIA_DIR, "model.jl"))
include(joinpath(MSTP_JULIA_DIR, "api.jl"))
