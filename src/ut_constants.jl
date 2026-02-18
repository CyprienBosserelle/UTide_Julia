module UtConstants
using ..Utilities
export ut_constants, constit_index_dict, cycles_per_hour, hours_per_cycle
_base_dir = joinpath(@__DIR__, "..", "data")
_ut_constants_fname = joinpath(_base_dir, "ut_constants.mat")
ut_constants = loadbunch(_ut_constants_fname, masked=false)
ut_constants = convert_unicode_arrays(ut_constants)
constit_names = ut_constants.const.name
constit_index_dict = Dict(name => i for (i, name) in enumerate(vec(constit_names)))
_uc = ut_constants.const
cycles_per_hour = Bunch(Dict(Symbol(name) => freq for (name, freq) in zip(vec(_uc.name), vec(_uc.freq))))
hours_per_cycle = Bunch(Dict(Symbol(name) => 1.0/freq for (name, freq) in zip(vec(_uc.name), vec(_uc.freq))))
end # module
