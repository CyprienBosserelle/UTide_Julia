module UtConstants
include("ut_constants_data.jl")
using .UtConstantsData

export ut_constants, constit_index_dict, cycles_per_hour, hours_per_cycle

const ut_constants = UtConstantsData.UT_CONSTANTS
const constit_names = ut_constants.const.name
const constit_index_dict = Dict(name => i for (i, name) in enumerate(vec(constit_names)))
const _uc = ut_constants.const
const cycles_per_hour = Dict(Symbol(name) => freq for (name, freq) in zip(vec(_uc.name), vec(_uc.freq)))
const hours_per_cycle = Dict(Symbol(name) => 1.0/freq for (name, freq) in zip(vec(_uc.name), vec(_uc.freq)))

end # module
