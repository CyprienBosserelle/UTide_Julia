module Solve
using ..Utilities, ..TimeConversion, ..Confidence, ..ConstituentSelection, ..Diagnostics, ..EllipseParams, ..Harmonics, ..RobustFit, LinearAlgebra, Statistics
export solve
const default_opts = Dict(:constit => "auto", :order_constit => nothing, :conf_int => "linear", :method => "ols", :trend => true, :phase => "Greenwich", :nodal => true, :infer => nothing, :MC_n => 200, :Rayleigh_min => 1.0, :robust_kw => Dict(:weight_function => "cauchy"), :white => false, :verbose => true, :epoch => nothing)
function _process_opts(opts, is_2D)
    merged_opts = Dict{Symbol, Any}(default_opts)
    for (k, v) in opts merged_opts[Symbol(k)] = v end
    merged_opts[:order_constit] = isnothing(merged_opts[:order_constit]) ? "PE" : merged_opts[:order_constit]
    return _translate_opts(merged_opts)
end
function _translate_opts(opts)
    notrend = !get(opts, :trend, true)
    cnstit = get(opts, :constit, "auto")
    ordercnstit = get(opts, :order_constit, "PE")
    infer = get(opts, :infer, nothing)
    conf_int = true
    linci = true
    nodiagn = 0

    ci_val = get(opts, :conf_int, "linear")
    if ci_val == "MC"
        linci = false
    elseif ci_val == "none"
        conf_int = false
        nodiagn = 1
    end

    nodsatlint = false
    nodsatnone = false
    nodal_val = get(opts, :nodal, true)
    if nodal_val == "linear_time"
        nodsatlint = true
    elseif nodal_val == false
        nodsatnone = true
    end

    gwchlint = false
    gwchnone = false
    phase_val = get(opts, :phase, "Greenwich")
    if phase_val == "linear_time"
        gwchlint = true
    elseif phase_val == "raw"
        gwchnone = true
    end

    rmin = get(opts, :Rayleigh_min, 1.0)
    white = get(opts, :white, false)
    RunTimeDisp = get(opts, :verbose, true)
    epoch = get(opts, :epoch, nothing)

    return Opt(
        constfolder="",
        twodim=false,
        minsnr=2.0,
        minpe=0.0,
        nodsatlint=nodsatlint,
        nodsatnone=nodsatnone,
        gwchlint=gwchlint,
        gwchnone=gwchnone,
        notrend=notrend,
        cnstit=cnstit,
        ordercnstit=ordercnstit,
        infer=infer,
        conf_int=conf_int,
        linci=linci,
        rmin=rmin,
        white=white,
        RunTimeDisp=RunTimeDisp,
        epoch=epoch,
        nodiagn=nodiagn,
        equi=false,
        lsfrqosmp=1,
        nrlzn=200,
        method=get(opts, :method, "ols"),
        robust_kw=get(opts, :robust_kw, Dict{Symbol, Any}(:weight_function => "cauchy"))
    )
end
function solve(t, u, v=nothing; lat=nothing, opts...)
    if isnothing(lat) error("Latitude must be supplied") end
    compat_opts = _process_opts(opts, !isnothing(v))
    return _solv1(t, u, v, lat, compat_opts)
end
function _solv1(tin, uin, vin, lat, opt)
    packed = _slvinit(tin, uin, vin, lat, opt)
    tin, t, u, v, tref, lor, elor, opt = packed
    nt = length(t)
    if opt.RunTimeDisp print("solve: ") end
    cnstit, temp_coef = ut_cnstitsel(tref, opt.rmin / (24.0 * lor), opt.cnstit, opt.infer)
    if isnothing(vin)
        coef = ScalarCoef(
            name=temp_coef.name,
            freq=temp_coef.aux.frq,
            lind=temp_coef.aux.lind,
            lat=lat,
            reftime=tref,
            opt=opt,
            nR=temp_coef.nR,
            nNR=temp_coef.nNR,
            nI=temp_coef.nI
        )
    else
        coef = VectorCoef(
            name=temp_coef.name,
            freq=temp_coef.aux.frq,
            lind=temp_coef.aux.lind,
            lat=lat,
            reftime=tref,
            opt=opt,
            nR=temp_coef.nR,
            nNR=temp_coef.nNR,
            nI=temp_coef.nI
        )
    end
    if opt.RunTimeDisp print("matrix prep ... ") end
    ngflgs = [opt.nodsatlint, opt.nodsatnone, opt.gwchlint, opt.gwchnone]
    E = ut_E(t, tref, cnstit.NR.frq, cnstit.NR.lind, lat, ngflgs, nothing)
    B = [E conj.(E)]

    if !isnothing(opt.infer)
        Etilp = zeros(ComplexF64, nt, coef.nR)
        Etilm = zeros(ComplexF64, nt, coef.nR)
        for (k, ref) in enumerate(cnstit.R)
            E_ref = ut_E(t, tref, [ref.frq], [ref.lind], lat, ngflgs, nothing)
            Q = ut_E(t, tref, ref.I.frq, ref.I.lind, lat, ngflgs, nothing) ./ E_ref
            Qsum_p = sum(Q .* ref.I.Rp', dims=2)
            Etilp[:, k] = E_ref[:, 1] .* (1.0 .+ vec(Qsum_p))
            Qsum_m = sum(Q .* conj.(ref.I.Rm)', dims=2)
            Etilm[:, k] = E_ref[:, 1] .* (1.0 .+ vec(Qsum_m))
        end
        B = [B Etilp conj.(Etilm)]
    end

    B = [B ones(nt, 1)]
    if !opt.notrend B = [B ((t .- tref) ./ lor)] end
    if opt.RunTimeDisp print("solution ... ") end
    xraw = !isnothing(v) ? u .+ 1im .* v : u
    if opt.method == "ols" m, W = B \ xraw, ones(nt)
    else rf = robustfit(B, xraw; opt.robust_kw...); m, W = rf.b, rf.w; coef.rf = rf end
    coef.weights, xmod = W, B * m
    if isnothing(v) xmod = real(xmod) end
    e, nI, nR, nNR = W .* (xraw .- xmod), coef.nI, coef.nR, coef.nNR
    ap, am = vcat(m[1:nNR], m[2*nNR+1 : 2*nNR+nR]), vcat(m[nNR+1 : 2*nNR], m[2*nNR+nR+1 : 2*nNR+2*nR])
    Xu, Yu = real(ap .+ am), -imag(ap .- am)
    if isnothing(v) coef.A, _, _, coef.g = ut_cs2cep(Xu, Yu); Xv, Yv = Float64[], Float64[]
    else Xv, Yv = imag(ap .+ am), real(ap .- am); coef.Lsmaj, coef.Lsmin, coef.theta, coef.g = ut_cs2cep(Xu, Yu, Xv, Yv) end

    if !isnothing(opt.infer)
        apI, amI = zeros(ComplexF64, nI), zeros(ComplexF64, nI)
        curr_ind = 1
        for k in 1:nR
            ref = cnstit.R[k]
            apI[curr_ind : curr_ind + ref.nI - 1] = ref.I.Rp .* ap[nNR + k]
            amI[curr_ind : curr_ind + ref.nI - 1] = ref.I.Rm .* am[nNR + k]
            curr_ind += ref.nI
        end
        XuI, YuI = real(apI .+ amI), -imag(apI .- amI)
        if isnothing(v)
            A_I, _, _, g_I = ut_cs2cep(XuI, YuI)
            coef.A, coef.g = vcat(coef.A, A_I), vcat(coef.g, g_I)
        else
            XvI, YvI = imag(apI .+ amI), real(apI .- amI)
            Lsmaj_I, Lsmin_I, theta_I, g_I = ut_cs2cep(XuI, YuI, XvI, YvI)
            coef.Lsmaj, coef.Lsmin, coef.theta, coef.g = vcat(coef.Lsmaj, Lsmaj_I), vcat(coef.Lsmin, Lsmin_I), vcat(coef.theta, theta_I), vcat(coef.g, g_I)
        end
    end
    if !isnothing(v)
        if opt.notrend coef.umean, coef.vmean = real(m[end]), imag(m[end])
        else coef.umean, coef.vmean, coef.uslope, coef.vslope = real(m[end-1]), imag(m[end-1]), real(m[end])/lor, imag(m[end])/lor end
    else
        if opt.notrend coef.mean = real(m[end])
        else coef.mean, coef.slope = real(m[end-1]), real(m[end])/lor end
    end
    if opt.conf_int coef = _confidence(coef, cnstit, opt, t, e, tin, elor, xraw, xmod, W, m, B, Xu, Yu, Xv, Yv) end
    if opt.nodiagn == 0 coef = ut_diagn(coef); coef.PE, coef.SNR = _PE(coef), _SNR(coef) end
    coef = _reorder(coef, opt)
    if opt.RunTimeDisp println("done.") end
    return coef
end
function _reorder(coef::ScalarCoef, opt)
    if opt.ordercnstit == "PE" ind = sortperm(coef.PE, rev=true)
    elseif opt.ordercnstit == "frequency" ind = sortperm(coef.freq)
    elseif opt.ordercnstit == "SNR" ind = sortperm(coef.SNR, rev=true)
    else ind = 1:length(coef.name) end

    coef.name = coef.name[ind]
    coef.PE = coef.PE[ind]
    coef.SNR = coef.SNR[ind]
    coef.A = coef.A[ind]
    coef.A_ci = coef.A_ci[ind]
    coef.g = coef.g[ind]
    coef.g_ci = coef.g_ci[ind]
    coef.freq = coef.freq[ind]
    coef.lind = coef.lind[ind]
    return coef
end

function _reorder(coef::VectorCoef, opt)
    if opt.ordercnstit == "PE" ind = sortperm(coef.PE, rev=true)
    elseif opt.ordercnstit == "frequency" ind = sortperm(coef.freq)
    elseif opt.ordercnstit == "SNR" ind = sortperm(coef.SNR, rev=true)
    else ind = 1:length(coef.name) end

    coef.name = coef.name[ind]
    coef.PE = coef.PE[ind]
    coef.SNR = coef.SNR[ind]
    coef.Lsmaj = coef.Lsmaj[ind]
    coef.Lsmaj_ci = coef.Lsmaj_ci[ind]
    coef.Lsmin = coef.Lsmin[ind]
    coef.Lsmin_ci = coef.Lsmin_ci[ind]
    coef.theta = coef.theta[ind]
    coef.theta_ci = coef.theta_ci[ind]
    coef.g = coef.g[ind]
    coef.g_ci = coef.g_ci[ind]
    coef.freq = coef.freq[ind]
    coef.lind = coef.lind[ind]
    return coef
end
function _slvinit(tin, uin, vin, lat, opts)
    if ndims(tin) != 1 || ndims(uin) != 1
        error("t and u must be 1-D arrays")
    end
    if !isnothing(vin) && ndims(vin) != 1
        error("v must be a 1-D array")
    end

    tin = normalize_time(tin, opts.epoch)

    # Step 1: remove invalid times from tin, uin, vin
    # Using isnan for Float64
    valid_t_mask = .!isnan.(tin)
    if !isnothing(uin) valid_t_mask .&= .!isnan.(uin) end
    if !isnothing(vin) valid_t_mask .&= .!isnan.(vin) end

    tin_valid = tin[valid_t_mask]
    uin_valid = isnothing(uin) ? nothing : uin[valid_t_mask]
    vin_valid = isnothing(vin) ? nothing : vin[valid_t_mask]

    t, u, v = tin_valid, uin_valid, vin_valid

    # Are the times equally spaced?
    if length(tin_valid) > 1 && var(diff(tin_valid)) < eps(Float64)
        opt_equi = true
        lor = tin_valid[end] - tin_valid[1]
        ntgood = length(tin_valid)
        elor = lor * ntgood / (ntgood - 1)
        tref = 0.5 * (tin_valid[1] + tin_valid[end])
    else
        opt_equi = false
        lor = maximum(t) - minimum(t)
        nt = length(t)
        elor = lor * nt / (nt - 1)
        tref = 0.5 * (minimum(t) + maximum(t))
    end

    opt = opts
    opt.twodim, opt.equi, opt.lsfrqosmp, opt.nrlzn = !isnothing(vin), opt_equi, 1, 200
    return tin_valid, t, u, v, tref, lor, elor, opt
end
end # module
