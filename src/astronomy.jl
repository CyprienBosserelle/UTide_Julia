module Astronomy
using LinearAlgebra
export ut_astron
const _sc = [270.434164, 13.1763965268, -0.0000850, 0.000000039]
const _hc = [279.696678, 0.9856473354, 0.00002267, 0.000000000]
const _pc = [334.329556, 0.1114040803, -0.0007739, -0.00000026]
const _npc = [-259.183275, 0.0529539222, -0.0001557, -0.000000050]
const _ppc = [281.220844, 0.0000470684, 0.0000339, 0.000000070]
const _coefs = vcat(_sc', _hc', _pc', _npc', _ppc')
function ut_astron(jd)
    jd_vec = collect(jd)[:]
    daten = 693595.5
    d = jd_vec .- daten
    D = d ./ 10000
    args = vcat(ones(length(jd_vec))', d', (D .* D)', (D.^3)')
    astro = mod.((_coefs * args) ./ 360, 1)
    tau = (jd_vec .% 1) .+ astro[2, :] .- astro[1, :]
    astro = vcat(tau', astro)
    dargs = vcat(zeros(length(jd_vec))', ones(length(jd_vec))', (2.0e-4 .* D)', (3.0e-4 .* D .* D)')
    ader = (_coefs * dargs) ./ 360.0
    dtau = 1.0 .+ ader[2, :] .- ader[1, :]
    ader = vcat(dtau', ader)
    return astro, ader
end
end # module
