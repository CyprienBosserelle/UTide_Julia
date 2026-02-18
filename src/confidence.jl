module Confidence
using ..Utilities, ..EllipseParams, ..Periodogram, Statistics, LinearAlgebra, Distributions
export _confidence, ut_linci
function band_averaged_psd_by_constit(tin, t, e, elor, coef, opt)
    constits = coef.aux.frq
    e_ = (opt.equi && length(tin) > length(t)) ? complex_interp(tin, t, e) : e
    ba = band_psd(tin, e_, constits, equi=opt.equi, frqosamp=opt.lsfrqosmp)
    df = 1.0 / (elor * 24.0)
    ba.Puu .*= df
    Puu = zeros(length(constits))
    Pvv = opt.twodim ? zeros(length(constits)) : nothing
    Puv = opt.twodim ? zeros(length(constits)) : nothing
    if opt.twodim ba.Pvv .*= df; ba.Puv .*= df end
    for i in 1:size(ba.fbnd, 1)
        lo, hi = ba.fbnd[i, 1], ba.fbnd[i, 2]
        inside = findall(f -> f >= lo && f <= hi, constits)
        Puu[inside] .= ba.Puu[i]
        if opt.twodim Pvv[inside] .= ba.Pvv[i]; Puv[inside] .= ba.Puv[i] end
    end
    return Puu, Pvv, Puv
end
function cluster(x, ang=360.0) y = copy(x); ofs = -x[1] + ang / 2.0; return mod.(x .+ ofs, ang) .- ofs end
function _is_PD(A) try cholesky(A); return true catch; return false end end
function nearestSPD(A)
    B = (A + A') / 2.0
    s = svd(B); H = s.V * diagm(abs.(s.S)) * s.Vt
    Ahat = (B + H) / 2.0; Ahat = (Ahat + Ahat') / 2.0
    n, k = size(A, 1), 0
    while k == 0 || !_is_PD(Ahat)
        k += 1; maxeig = maximum(abs.(eigvals(Ahat)))
        for i in 1:n Ahat[i, i] += eps(maxeig) end
        if k > 100 @warn "nearestSPD did not converge"; return diagm(diag(A)) end
    end
    return Ahat
end
function _confidence(coef, cnstit, opt, t, e, tin, elor, xraw, xmod, W, m, B, Xu, Yu, Xv, Yv)
    Puu, Pvv, Puv = opt.white ? (nothing, nothing, nothing) : band_averaged_psd_by_constit(tin, t, e, elor, coef, opt)
    _Wx, _WB = W .* xraw, W .* B
    nt, nm = length(xraw), size(B, 2)
    varMSM = real(dot(xraw, _Wx) - dot(xmod, _Wx)) / (nt - nm)
    gamC = inv(B' * _WB) * varMSM
    gamP = inv(transpose(B) * _WB) * ((dot(transpose(xraw), _Wx) - dot(transpose(xmod), _Wx)) / (nt - nm))
    Gall, Hall = gamC + gamP, gamC - gamP
    nc = length(Xu)
    coef.g_ci = fill(NaN, size(coef.g))
    if opt.twodim coef.Lsmaj_ci, coef.Lsmin_ci, coef.theta_ci, varcov_mCw = copy(coef.g_ci), copy(coef.g_ci), copy(coef.g_ci), fill(NaN, nc, 4, 4)
    else coef.A_ci, varcov_mCw = copy(coef.g_ci), fill(NaN, nc, 2, 2) end
    if !opt.white varcov_mCc = copy(varcov_mCw) end
    for c in 1:nc
        G = [Gall[c, c] Gall[c, c+nc]; Gall[c+nc, c] Gall[c+nc, c+nc]]
        H = [Hall[c, c] Hall[c, c+nc]; Hall[c+nc, c] Hall[c+nc, c+nc]]
        varXu, varYu = real(G[1, 1] + G[2, 2] + 2 * G[1, 2]) / 2.0, real(H[1, 1] + H[2, 2] - 2 * H[1, 2]) / 2.0
        if opt.twodim varXv, varYv = real(H[1, 1] + H[2, 2] + 2 * H[1, 2]) / 2.0, real(G[1, 1] + G[2, 2] - 2 * G[1, 2]) / 2.0 end
        if opt.linci
            if !opt.twodim
                varcov_mCw[c, :, :] = diagm([varXu, varYu])
                if !opt.white den = varXu + varYu; varXu, varYu = Puu[c] * varXu / den, Puu[c] * varYu / den; varcov_mCc[c, :, :] = diagm([varXu, varYu]) end
                sig1, sig2 = ut_linci(Xu[c], Yu[c], sqrt(max(0, varXu)), sqrt(max(0, varYu)))
                coef.A_ci[c], coef.g_ci[c] = 1.96 * sig1, 1.96 * sig2
            else
                varcov_mCw[c, :, :] = diagm([varXu, varYu, varXv, varYv])
                if !opt.white den = varXv + varYv; varXv, varYv = Pvv[c] * varXv / den, Pvv[c] * varYv / den; varcov_mCc[c, :, :] = diagm([varXu, varYu, varXv, varYv]) end
                sig1, sig2 = ut_linci(Xu[c] + 1im * Xv[c], Yu[c] + 1im * Yv[c], sqrt(max(0, varXu)) + 1im * sqrt(max(0, varXv)), sqrt(max(0, varYu)) + 1im * sqrt(max(0, varYv)))
                coef.Lsmaj_ci[c], coef.Lsmin_ci[c], coef.g_ci[c], coef.theta_ci[c] = 1.96 * real(sig1), 1.96 * imag(sig1), 1.96 * real(sig2), 1.96 * imag(sig2)
            end
        else
            # MC
            covXuYu = imag(H[1, 1] - H[1, 2] + H[2, 1] - H[2, 2]) / 2.0
            Duu = [varXu covXuYu; covXuYu varYu]
            if !opt.twodim
                Duu = !opt.white ? nearestSPD(Puu[c] * Duu / tr(Duu)) : nearestSPD(Duu)
                mCall = rand(MvNormal([Xu[c], Yu[c]], Duu), opt.nrlzn)'
                A_mc, _, _, g_mc = ut_cs2cep(mCall)
                coef.A_ci[c], coef.g_ci[c] = 1.96 * median(abs.(A_mc .- median(A_mc))) / 0.6745, 1.96 * median(abs.(cluster(g_mc, 360.0) .- median(cluster(g_mc, 360.0)))) / 0.6745
            end
        end
    end
    if !isnothing(opt.infer)
        nNR = coef.nNR
        ind = nNR + coef.nR + 1
        for k in 1:coef.nR
            ref = cnstit.R[k]
            varcov = opt.white ? varcov_mCw : varcov_mCc
            idx = nNR + k
            varReap = 0.25 * varcov[idx, 1, 1]
            varImap = 0.25 * varcov[idx, 2, 2]
            if opt.twodim
                varReap += 0.25 * varcov[idx, 4, 4]
                varImap += 0.25 * varcov[idx, 3, 3]
            end
            for i in 1:length(ref.I.Rp)
                rp, rm = ref.I.Rp[i], ref.I.Rm[i]
                varX = (real(rp)^2 + real(rm)^2) * varReap + (imag(rp)^2 + imag(rm)^2) * varImap
                varY = (real(rp)^2 + real(rm)^2) * varImap + (imag(rp)^2 + imag(rm)^2) * varReap
                if !opt.twodim
                    sig1, sig2 = ut_linci(Xu[idx], Yu[idx], sqrt(max(0, varX)), sqrt(max(0, varY)))
                    coef.A_ci[ind], coef.g_ci[ind] = 1.96 * sig1, 1.96 * sig2
                else
                    sig1, sig2 = ut_linci(Xu[idx] + 1im * Xv[idx], Yu[idx] + 1im * Yv[idx],
                                         sqrt(max(0, varX)) + 1im * sqrt(max(0, varY)),
                                         sqrt(max(0, varY)) + 1im * sqrt(max(0, varX)))
                    coef.Lsmaj_ci[ind], coef.Lsmin_ci[ind], coef.g_ci[ind], coef.theta_ci[ind] = 1.96 * real(sig1), 1.96 * imag(sig1), 1.96 * real(sig2), 1.96 * imag(sig2)
                end
                ind += 1
            end
        end
    end
    return coef
end
function ut_linci(X, Y, sigX, sigY)
    Xu, sigXu, Yu, sigYu, Xv, sigXv, Yv, sigYv = real(X), real(sigX), real(Y), real(sigY), imag(X), imag(sigX), imag(Y), imag(sigY)
    rp, rm = 0.5 * sqrt((Xu + Yv)^2 + (Xv - Yu)^2), 0.5 * sqrt((Xu - Yv)^2 + (Xv + Yu)^2)
    sigXu2, sigYu2, sigXv2, sigYv2 = sigXu^2, sigYu^2, sigXv^2, sigYv^2
    ex, fx, gx, hx = (Xu + Yv) / rp, (Xu - Yv) / rm, (Yu - Xv) / rp, (Yu + Xv) / rm
    sig1 = sqrt((0.25 * (ex + fx))^2 * sigXu2 + (0.25 * (gx + hx))^2 * sigYu2 + (0.25 * (hx - gx))^2 * sigXv2 + (0.25 * (ex - fx))^2 * sigYv2)
    rn, rd = 2.0 * (Xu * Yu + Xv * Yv), Xu^2 - Yu^2 + Xv^2 - Yv^2
    den = rn^2 + rd^2
    sig2 = (180.0 / pi) * sqrt(((rd * Yu - rn * Xu) / den)^2 * sigXu2 + ((rd * Xu + rn * Yu) / den)^2 * sigYu2 + ((rd * Yv - rn * Xv) / den)^2 * sigXv2 + ((rd * Xv + rn * Yv) / den)^2 * sigYv2)
    if X isa Complex
        sig1 += 1im * sqrt((0.25 * (ex - fx))^2 * sigXu2 + (0.25 * (gx - hx))^2 * sigYu2 + (0.25 * (hx + gx))^2 * sigXv2 + (0.25 * (ex + fx))^2 * sigYv2)
        rn, rd = 2.0 * (Xu * Xv + Yu * Yv), Xu^2 + Yu^2 - (Xv^2 + Yv^2)
        den = rn^2 + rd^2
        sig2 += 1im * (180.0 / pi) * sqrt(((rd * Xv - rn * Xu) / den)^2 * sigXu2 + ((rd * Yv - rn * Yu) / den)^2 * sigYu2 + ((rd * Xu + rn * Xv) / den)^2 * sigXv2 + ((rd * Yu + rn * Yv) / den)^2 * sigYv2)
    end
    return sig1, sig2
end
end # module
