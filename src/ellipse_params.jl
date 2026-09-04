module EllipseParams
export ut_cs2cep
function ut_cs2cep(Xu, Yu=nothing, Xv=nothing, Yv=nothing)
    if isnothing(Yu)
        if size(Xu, 2) == 2
            Yu, Xu = Xu[:, 2], Xu[:, 1]
        elseif size(Xu, 2) == 4
            Yu, Xv, Yv, Xu = Xu[:, 2], Xu[:, 3], Xu[:, 4], Xu[:, 1]
        else error("invalid arguments") end
    end
    if isnothing(Xv)
        ap = Xu .- 1im .* Yu
        Lsmaj, Lsmin, theta = abs.(ap), zeros(eltype(Xu), size(Xu)), zeros(eltype(Xu), size(Xu))
        g = mod.(-rad2deg.(angle.(ap)), 360.0)
        return Lsmaj, Lsmin, theta, g
    end
    ap, am = ((Xu .+ Yv) .+ 1im .* (Xv .- Yu)) ./ 2.0, ((Xu .- Yv) .+ 1im .* (Xv .+ Yu)) ./ 2.0
    Ap, Am = abs.(ap), abs.(am)
    Lsmaj, Lsmin = Ap .+ Am, Ap .- Am
    epsp, epsm = rad2deg.(angle.(ap)), rad2deg.(angle.(am))
    theta = mod.((epsp .+ epsm) ./ 2.0, 180.0)
    g = mod.(-epsp .+ theta, 360.0)
    # Wait, theta was not defined in the g calculation line in my previous version?
    # Ah, g = mod.(-epsp .+ theta, 360.0). Correct.
    return Lsmaj, Lsmin, theta, g
end
end # module
