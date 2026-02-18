module Periodogram
using ..Utilities, Statistics, FFTW, DSP, LombScargle, Interpolations
export freq_bands, band_psd
const freq_bands = [0.00010 0.00417; 0.03192 0.04859; 0.07218 0.08884; 0.11243 0.12910; 0.15269 0.16936; 0.19295 0.20961; 0.23320 0.25100; 0.26000 0.29000; 0.30000 0.50000]
function fbndavg(P, freq, cfreq=nothing, fbands=nothing)
    if isnothing(fbands) fbands = freq_bands end
    P_work = copy(P)
    if !isnothing(cfreq)
        itp = linear_interpolation(freq, 1:length(freq), extrapolation_bc=Flat())
        i_constit = clamp.(round.(Int, itp(cfreq)), 1, length(freq))
        P_work[i_constit] .= NaN
    end
    nbands = size(fbands, 1)
    avP = fill(convert(eltype(P), NaN), nbands)
    for k in 1:nbands
        indices = findall(f -> f >= fbands[k, 1] && f <= fbands[k, 2], freq)
        if !isempty(indices) avP[k] = nanmean(P_work[indices]) end
    end
    return avP
end
nanmean(x) = (v = filter(!isnan, x); isempty(v) ? NaN : mean(v))
function _lomb_freqs(t, fbands=nothing; ofac=1, max_per_band=500)
    if isnothing(fbands) fbands = freq_bands end
    n = length(t)
    reclen = n * (t[end] - t[1]) / (n - 1)
    freq = collect(1:(floor(Int, n*ofac/2)-1)) ./ reclen
    freqs = Float64[]
    for k in 1:size(fbands, 1)
        indices = findall(f -> f >= fbands[k, 1] && f <= fbands[k, 2], freq)
        if length(indices) > max_per_band append!(freqs, range(freq[indices[1]], freq[indices[end]], length=max_per_band))
        elseif !isempty(indices) append!(freqs, freq[indices]) end
    end
    return freqs
end
function _psd_lomb(t, x, window=nothing, freq=nothing; ofac=1)
    out, xdm, n = Bunch(), x .- mean(x), length(x)
    if isnothing(window) w = ones(n)
    else itp = linear_interpolation(range(minimum(t), maximum(t), length=n), window, extrapolation_bc=Flat()); w = itp(t); xdm .*= w end
    delta_t = (t[end] - t[1]) / (n - 1)
    if isnothing(freq) freq = collect(1:(floor(Int, n*ofac/2)-1)) ./ (n * delta_t) end
    out.F, psdnorm = freq, 2.0 * delta_t * n / sum(w.^2)
    out.Pxx = psdnorm .* _lombscargle_match(t, real(xdm), freq)
    if eltype(x) <: Complex out.Pyy = psdnorm .* _lombscargle_match(t, imag(xdm), freq); out.Pxy = psdnorm .* _ls_cross(t, xdm, freq) end
    return out
end
function _lombscargle_match(t, x, freq)
    P = zeros(length(freq))
    for (i, f) in enumerate(freq)
        arg = (2.0 * pi * f) .* t
        tau = 0.5 * atan(sum(sin.(2.0 .* arg)), sum(cos.(2.0 .* arg)))
        c, s = cos.(arg .- tau), sin.(arg .- tau)
        P[i] = 0.5 * (sum(x .* c)^2 / sum(c.^2) + sum(x .* s)^2 / sum(s.^2))
    end
    return P
end
function _ls_cross(t, x, freq)
    nc, xr, xi = length(freq), real(x), imag(x)
    pxy = zeros(ComplexF64, nc)
    for (i, f) in enumerate(freq)
        arg = (2.0 * pi * f) .* t
        tau = 0.5 * atan(sum(sin.(2.0 .* arg)), sum(cos.(2.0 .* arg)))
        c, s = cos.(arg .- tau), sin.(arg .- tau)
        a, b = 1.0 / sqrt(sum(c.^2)), 1.0 / sqrt(sum(s.^2))
        tmpx, tmpy = (a * sum(xr .* c)) + 1im * (b * sum(xr .* s)), (a * sum(xi .* c)) + 1im * (-b * sum(xi .* s))
        f0 = exp(1im * tau); pxy[i] = 0.5 * (tmpx * f0 * tmpy * conj(f0))
    end
    return pxy
end
function _psd(e, window, fs)
    e_dm = (e .- mean(e)) .* window
    cs = eltype(e) <: Complex ? conj.(rfft(real(e_dm))) .* rfft(imag(e_dm)) : abs2.(rfft(real(e_dm)))
    if length(cs) > 1 cs[2:end-1] .*= 2.0 end
    return cs .* ((1.0 / fs) * (1.0 / sum(window.^2)))
end
function band_psd(t, e, cfrq; equi=true, frqosamp=1)
    P = Bunch(fbnd=freq_bands)
    if length(e) % 2 != 0 e, t = e[1:end-1], t[1:end-1] end
    nt, hn = length(e), hanning(length(e))
    local allfrq, Puu1s, ls_spec
    if equi
        dt = 24.0 * (t[2] - t[1])
        fs = 1.0 / dt
        Puu1s, allfrq = _psd(real(e), hn, fs), collect(0:(nt ÷ 2)) ./ (nt * dt)
    else
        ls_spec = _psd_lomb(t .* 24.0, e, hn, _lomb_freqs(t .* 24.0, ofac=frqosamp))
        Puu1s, allfrq = ls_spec.Pxx, ls_spec.F
    end
    P.Puu = fbndavg(Puu1s, allfrq, cfrq)
    if eltype(e) <: Complex
        if equi Pvv1s, Puv1s = _psd(imag(e), hn, fs), _psd(e, hn, fs)
        else Pvv1s, Puv1s = ls_spec.Pyy, ls_spec.Pxy end
        P.Pvv, P.Puv = fbndavg(Pvv1s, allfrq, cfrq), real.(fbndavg(Puv1s, allfrq, cfrq))
    end
    return P
end
end # module
