module Harmonics
using ..UtConstants, ..Astronomy, LinearAlgebra
export linearized_freqs, ut_E, FUV
const sat = UtConstants.ut_constants.sat
const const_data = UtConstants.ut_constants.const
const shallow = UtConstants.ut_constants.shallow
const _ishallow_raw = vec(const_data.ishallow)
const _nshallow_raw = vec(const_data.nshallow)
const _not_shallow = isnan.(_ishallow_raw)
const _nshallow = Int.(filter(!isnan, _nshallow_raw))
const _ishallow = Int.(filter(!isnan, _ishallow_raw)) .- 1
const _kshallow = findall(!, _not_shallow)
function linearized_freqs(tref)
    astro, ader = ut_astron(tref)
    freq = copy(const_data.freq)
    selected_indices = findall(_not_shallow)
    selected = (const_data.doodson[selected_indices, :] * ader) ./ 24.0
    freq[selected_indices] .= vec(selected)
    for (i, (i0, nshal, k)) in enumerate(zip(_ishallow, _nshallow, _kshallow))
        ik = i0 .+ (1:nshal)
        indices = round.(Int, shallow.iname[ik])
        freq[k] = sum(freq[indices] .* shallow.coef[ik])
    end
    return freq
end
function ut_E(t, tref, frq, lind, lat, ngflgs, prefilt)
    t, frq, lind = collect(t), collect(frq), collect(lind)
    nt, nc = length(t), length(frq)
    if ngflgs[2] && ngflgs[4]
        F, U = ones(nt, nc), zeros(nt, nc)
        V = (24.0 .* (t .- tref)) * frq'
    else
        F, U, V = FUV(t, tref, lind, lat, ngflgs)
    end
    return F .* exp.(1im .* (U .+ V) .* 2.0 .* pi)
end
function FUV(t, tref, lind, lat, ngflgs)
    t = collect(t)
    nt, nc = length(t), length(lind)
    local F, U
    if ngflgs[2]
        F, U = ones(nt, nc), zeros(nt, nc)
    else
        tt = ngflgs[1] ? [tref] : t
        ntt = length(tt)
        astro, ader = ut_astron(tt)
        actual_lat = abs(lat) < 5 ? sign(lat) * 5 : lat
        slat = sin(deg2rad(actual_lat))
        rr = copy(sat.amprat)
        j1 = sat.ilatfac .== 1
        rr[j1] .*= 0.36309 * (1.0 - 5.0 * slat^2) / slat
        j2 = sat.ilatfac .== 2
        rr[j2] .*= 2.59808 * slat
        uu = sat.deldood * astro[4:6, :] .+ sat.phcorr
        uu = mod.(uu, 1.0)
        mat = rr .* exp.(1im .* 2.0 .* pi .* uu)
        nfreq = length(const_data.isat)
        F_all = ones(ComplexF64, nfreq, ntt)
        iconst_1based = Int.(sat.iconst)
        # Group indices by iconst
        group_indices = Dict{Int, Vector{Int}}()
        for (idx, val) in enumerate(iconst_1based)
            push!(get!(group_indices, val, Int[]), idx)
        end
        for (ii, indices) in group_indices
            if ii <= nfreq
                for t_idx in 1:ntt
                    sum_val = 0.0 + 0.0im
                    for idx in indices
                        sum_val += mat[idx, t_idx]
                    end
                    F_all[ii, t_idx] = 1.0 + sum_val
                end
            end
        end
        U_all = angle.(F_all) ./ (2.0 .* pi)
        F_all = abs.(F_all)
        for (i0, nshal, k) in zip(_ishallow, _nshallow, _kshallow)
            ik = i0 .+ (1:nshal)
            j = round.(Int, shallow.iname[ik])
            exp1 = shallow.coef[ik]
            for t_idx in 1:ntt
                f_val = 1.0
                u_val = 0.0
                for (s_idx, j_idx) in enumerate(j)
                    f_val *= F_all[j_idx, t_idx] ^ abs(exp1[s_idx])
                    u_val += U_all[j_idx, t_idx] * exp1[s_idx]
                end
                F_all[k, t_idx] = f_val
                U_all[k, t_idx] = u_val
            end
        end
        F, U = F_all[lind, :]', U_all[lind, :]'
    end
    local V
    if ngflgs[4]
        freq = linearized_freqs(tref)
        V = (24.0 .* (t .- tref)) .* freq[lind]'
    else
        tt = ngflgs[3] ? [tref] : t
        astro, ader = ut_astron(tt)
        V_all = const_data.doodson * astro .+ const_data.semi
        V_all = mod.(V_all, 1.0)
        for (i0, nshal, k) in zip(_ishallow, _nshallow, _kshallow)
            ik = i0 .+ (1:nshal)
            j = round.(Int, shallow.iname[ik])
            exp1 = shallow.coef[ik]
            for t_idx in 1:ntt
                v_val = 0.0
                for (s_idx, j_idx) in enumerate(j)
                    v_val += V_all[j_idx, t_idx] * exp1[s_idx]
                end
                V_all[k, t_idx] = v_val
            end
        end
        V = V_all[lind, :]'
        if ngflgs[3]
            freq = linearized_freqs(tref)
            V = V .+ (24.0 .* (t .- tref)) .* freq[lind]'
        end
    end
    return F, U, V
end
end # module
