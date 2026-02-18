module ConstituentSelection
using ..UtConstants, ..Harmonics, ..Utilities, OrderedCollections
export ut_cnstitsel
function ut_cnstitsel(tref, minres, incnstit, infer)
    const_data = UtConstants.ut_constants.const
    cnstit, coef = Bunch(), Bunch()
    freqs = linearized_freqs(tref)
    cnstit.NR = Bunch()
    if incnstit == "auto"
        cnstit.NR.lind = findall(vec(const_data.df) .>= minres)
    else
        cnstit.NR.lind = [UtConstants.constit_index_dict[n] for n in incnstit]
    end
    if !isnothing(infer)
        RIset = Set(infer.inferred_names) ∪ Set(infer.reference_names)
        RI_index_set = Set([UtConstants.constit_index_dict[n] for n in RIset])
        cnstit.NR.lind = filter(ind -> !(ind in RI_index_set), cnstit.NR.lind)
    end
    cnstit.NR.frq = freqs[cnstit.NR.lind]
    cnstit.NR.name = const_data.name[cnstit.NR.lind]
    nNR = length(cnstit.NR.frq)
    nR, nI = 0, 0
    cnstit.R = []
    local allrefs
    if !isnothing(infer)
        nI = length(infer.inferred_names)
        allrefs = unique(infer.reference_names)
        nR = length(allrefs)
        for name in allrefs
            refstruct = Bunch(name=name)
            refstruct.lind = UtConstants.constit_index_dict[name]
            refstruct.frq = freqs[refstruct.lind]
            ind = findall(rname -> rname == name, infer.reference_names)
            refstruct.nI = length(ind)
            refstruct.I = Bunch(Rp=ComplexF64[], Rm=ComplexF64[], name=String[], lind=Int[], frq=Float64[])
            for (lk, ilk) in enumerate(ind)
                push!(refstruct.I.Rp, infer.amp_ratios[ilk] * exp(1im * infer.phase_offsets[ilk] * pi / 180.0))
                if length(infer.amp_ratios) > nI
                    push!(refstruct.I.Rm, infer.amp_ratios[ilk + nI] * exp(-1im * infer.phase_offsets[ilk + nI] * pi / 180.0))
                else push!(refstruct.I.Rm, conj(refstruct.I.Rp[lk])) end
                iname = infer.inferred_names[ilk]
                push!(refstruct.I.name, iname)
                lind = UtConstants.constit_index_dict[iname]
                push!(refstruct.I.lind, lind)
                push!(refstruct.I.frq, freqs[lind])
            end
            push!(cnstit.R, refstruct)
        end
    end
    coef_name = collect(cnstit.NR.name)
    coef_aux = Bunch(frq = collect(cnstit.NR.frq), lind = collect(cnstit.NR.lind), reftime = tref)
    if !isnothing(infer)
        append!(coef_name, allrefs)
        append!(coef_aux.frq, [ref.frq for ref in cnstit.R])
        append!(coef_aux.lind, [ref.lind for ref in cnstit.R])
        for ref in cnstit.R
            append!(coef_name, ref.I.name)
            append!(coef_aux.frq, ref.I.frq)
            append!(coef_aux.lind, ref.I.lind)
        end
    end
    coef.name, coef.aux = coef_name, coef_aux
    coef.nR, coef.nNR, coef.nI = nR, nNR, nI
    return cnstit, coef
end
end # module
