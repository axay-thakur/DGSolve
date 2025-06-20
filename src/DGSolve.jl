module DGSolve
using Reexport
@reexport using LinearAlgebra # For inv, transpose, eigen, diagm, I
@reexport using Statistics # For mean, etc. (not directly used in current code, but good to keep if needed)
@reexport using Plots # For plotting functions (plot, plot!, savefig, display)
@reexport using Polynomials # For polyfit and polyval equivalents (fit, Polynomial, .)
@reexport using SpecialFunctions # For special functions like exp, erf, etc.

# Include all source files in dependency order
include("mesh_utils.jl")
include("quadrature.jl")
include("basis_functions.jl")
include("mass_matrix.jl")
include("claw.jl")
include("dg_operators.jl")
include("time_stepping.jl") # Depends on dg_operators for advance_ode_nonlin_rk4's odefun argument
include("plotting.jl") # Depends on Polynomials
include("../examples.jl") # Main script for running simulation

function __init__()
println("""
+++++++++++++++++++++++
=======DGSolve========
+++++++++++++++++++++++
""")
return nothing
end

end # module DGSolver