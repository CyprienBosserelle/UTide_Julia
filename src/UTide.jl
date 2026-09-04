module UTide
include("utilities.jl")
include("time_conversion.jl")
include("astronomy.jl")
include("ut_constants.jl")
include("harmonics.jl")
include("ellipse_params.jl")
include("robustfit.jl")
include("periodogram.jl")
include("constituent_selection.jl")
include("diagnostics.jl")
include("confidence.jl")
include("solve.jl")
include("reconstruct.jl")
using .Solve, .Reconstruct, .UtConstants
using .Utilities: Opt, ScalarCoef, VectorCoef, ReconstructOutput
export solve, reconstruct, ut_constants, Opt, ScalarCoef, VectorCoef, ReconstructOutput
end # module
