module Reconstruct
using ..Utilities, ..TimeConversion, ..Harmonics, LinearAlgebra
export reconstruct
function reconstruct(t, coef; epoch=nothing, verbose=true, constit=nothing, min_SNR=2.0, min_PE=0.0)
    t_norm = normalize_time(t, epoch)
    t_mpl = t_norm
    goodmask = .!isnan.(t_norm)
    t_comp = t_norm[goodmask]
    u, v = _reconstruct_inner(t_comp, goodmask, coef, verbose, constit, min_SNR, min_PE)

    return ReconstructOutput(
        t_in = t,
        epoch = epoch,
        constit = constit,
        min_SNR = min_SNR,
        min_PE = min_PE,
        t_mpl = t_mpl,
        h = isnothing(v) ? u : nothing,
        u = !isnothing(v) ? u : nothing,
        v = v
    )
end
function _reconstruct_inner(t, goodmask, coef, verbose, constit, min_SNR, min_PE)
    twodim = coef.aux.opt.twodim
    if !isnothing(constit) ind = filter(!isnothing, [findfirst(==(c), coef.name) for c in constit])
    elseif (min_SNR == 0 && min_PE == 0) || (coef.aux.opt.nodiagn != 0) ind = 1:length(coef.name)
    else
        E = twodim ? coef.Lsmaj .^ 2 .+ coef.Lsmin .^ 2 : coef.A .^ 2
        N = twodim ? (coef.Lsmaj_ci ./ 1.96) .^ 2 .+ (coef.Lsmin_ci ./ 1.96) .^ 2 : (coef.A_ci ./ 1.96) .^ 2
        SNR, PE = E ./ N, 100.0 .* E ./ sum(E)
        ind = findall(SNR .>= min_SNR .&& PE .>= min_PE)
    end
    rpd = pi / 180.0
    if twodim
        ap = 0.5 .* ((coef.Lsmaj[ind] .+ coef.Lsmin[ind]) .* exp.(1im .* (coef.theta[ind] .- coef.g[ind]) .* rpd))
        am = 0.5 .* ((coef.Lsmaj[ind] .- coef.Lsmin[ind]) .* exp.(1im .* (coef.theta[ind] .+ coef.g[ind]) .* rpd))
    else ap = 0.5 .* coef.A[ind] .* exp.(-1im .* coef.g[ind] .* rpd); am = conj.(ap) end
    ngflgs = [coef.aux.opt.nodsatlint, coef.aux.opt.nodsatnone, coef.aux.opt.gwchlint, coef.aux.opt.gwchnone]
    if verbose print("prep/calcs ... ") end
    E_mat = ut_E(t, coef.aux.reftime, coef.aux.frq[ind], coef.aux.lind[ind], coef.aux.lat, ngflgs, nothing)
    fit = E_mat * ap .+ conj.(E_mat) * am
    u, trend = fill(NaN, size(goodmask)), !coef.aux.opt.notrend
    local v = nothing
    if twodim
        v = fill(NaN, size(goodmask)); u[goodmask], v[goodmask] = real(fit) .+ coef.umean, imag(fit) .+ coef.vmean
        if trend u[goodmask] .+= coef.uslope .* (t .- coef.aux.reftime); v[goodmask] .+= coef.vslope .* (t .- coef.aux.reftime) end
    else u[goodmask] = real(fit) .+ coef.mean; if trend u[goodmask] .+= coef.slope .* (t .- coef.aux.reftime) end end
    if verbose println("done.") end
    return u, v
end
end # module
