module Diagnostics
export _PE, _SNR, ut_diagn
function _PE(coef)
    if haskey(coef, :Lsmaj) E = coef.Lsmaj .^ 2 .+ coef.Lsmin .^ 2; PE = (100.0 / sum(E)) .* E
    else PE = 100.0 .* coef.A .^ 2 ./ sum(coef.A .^ 2) end
    return PE
end
function _SNR(coef)
    if haskey(coef, :Lsmaj) SNR = (coef.Lsmaj .^ 2 .+ coef.Lsmin .^ 2) ./ ((coef.Lsmaj_ci ./ 1.96) .^ 2 .+ (coef.Lsmin_ci ./ 1.96) .^ 2)
    else SNR = (coef.A .^ 2) ./ (coef.A_ci ./ 1.96) .^ 2 end
    return SNR
end
function ut_diagn(coef)
    PE, SNR = _PE(coef), _SNR(coef)
    indPE = sortperm(PE, rev=true)
    coef.diagn = Dict(:name => coef.name[indPE], :PE => PE[indPE], :SNR => SNR[indPE])
    return coef
end
end # module
