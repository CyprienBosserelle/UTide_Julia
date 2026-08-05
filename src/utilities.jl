module Utilities
using LinearAlgebra, Statistics, OrderedCollections, Interpolations

export Bunch, complex_interp, loadbunch, convert_unicode_arrays, Opt, ScalarCoef, VectorCoef, ReconstructOutput

mutable struct Opt
    constfolder::String
    twodim::Bool
    minsnr::Float64
    minpe::Float64
    nodsatlint::Bool
    nodsatnone::Bool
    gwchlint::Bool
    gwchnone::Bool
    notrend::Bool
    cnstit::Any
    ordercnstit::String
    infer::Any
    conf_int::Bool
    linci::Bool
    rmin::Float64
    white::Bool
    RunTimeDisp::Bool
    epoch::Any
    nodiagn::Int
    equi::Bool
    lsfrqosmp::Int
    nrlzn::Int
    method::String
    robust_kw::Dict{Symbol, Any}
end

function Opt(;
    constfolder="",
    twodim=false,
    minsnr=2.0,
    minpe=0.0,
    nodsatlint=false,
    nodsatnone=false,
    gwchlint=false,
    gwchnone=false,
    notrend=false,
    cnstit="auto",
    ordercnstit="PE",
    infer=nothing,
    conf_int=true,
    linci=true,
    rmin=1.0,
    white=false,
    RunTimeDisp=true,
    epoch=nothing,
    nodiagn=0,
    equi=false,
    lsfrqosmp=1,
    nrlzn=200,
    method="ols",
    robust_kw=Dict{Symbol, Any}(:weight_function => "cauchy")
)
    return Opt(
        constfolder, twodim, minsnr, minpe, nodsatlint, nodsatnone, gwchlint, gwchnone, notrend,
        cnstit, ordercnstit, infer, conf_int, linci, rmin, white, RunTimeDisp, epoch, nodiagn,
        equi, lsfrqosmp, nrlzn, method, robust_kw
    )
end

function Opt(d::AbstractDict)
    constfolder = get(d, :constfolder, "")
    twodim = get(d, :twodim, false)
    minsnr = get(d, :minsnr, 2.0)
    minpe = get(d, :minpe, 0.0)
    nodsatlint = get(d, :nodsatlint, false)
    nodsatnone = get(d, :nodsatnone, false)
    gwchlint = get(d, :gwchlint, false)
    gwchnone = get(d, :gwchnone, false)
    notrend = get(d, :notrend, false)
    cnstit = get(d, :cnstit, "auto")
    ordercnstit = get(d, :ordercnstit, "PE")
    infer = get(d, :infer, nothing)
    conf_int = get(d, :conf_int, true)
    linci = get(d, :linci, true)
    rmin = get(d, :rmin, 1.0)
    white = get(d, :white, false)
    RunTimeDisp = get(d, :RunTimeDisp, true)
    epoch = get(d, :epoch, nothing)
    nodiagn = get(d, :nodiagn, 0)
    equi = get(d, :equi, false)
    lsfrqosmp = get(d, :lsfrqosmp, 1)
    nrlzn = get(d, :nrlzn, 200)

    method = "ols"
    robust_kw = Dict{Symbol, Any}(:weight_function => "cauchy")
    if haskey(d, :newopts)
        no = d[:newopts]
        method = get(no, :method, "ols")
        robust_kw = get(no, :robust_kw, robust_kw)
    elseif haskey(d, :method)
        method = d[:method]
        robust_kw = get(d, :robust_kw, robust_kw)
    end

    return Opt(
        constfolder, twodim, minsnr, minpe, nodsatlint, nodsatnone, gwchlint, gwchnone, notrend,
        cnstit, ordercnstit, infer, conf_int, linci, rmin, white, RunTimeDisp, epoch, nodiagn,
        equi, lsfrqosmp, nrlzn, method, robust_kw
    )
end

mutable struct ScalarCoef
    name::Vector{String}
    A::Vector{Float64}
    A_ci::Vector{Float64}
    g::Vector{Float64}
    g_ci::Vector{Float64}
    mean::Float64
    slope::Float64
    freq::Vector{Float64}
    lind::Vector{Int}
    lat::Float64
    reftime::Float64
    opt::Opt
    PE::Vector{Float64}
    SNR::Vector{Float64}
    diagn::Dict{Symbol, Any}
    weights::Vector{Float64}
    rf::Any
    nR::Int
    nNR::Int
    nI::Int

    function ScalarCoef(;
        name=String[], A=Float64[], A_ci=Float64[], g=Float64[], g_ci=Float64[],
        mean=0.0, slope=0.0, freq=Float64[], lind=Int[], lat=0.0, reftime=0.0,
        opt=Opt(), PE=Float64[], SNR=Float64[], diagn=Dict{Symbol, Any}(),
        weights=Float64[], rf=nothing, nR=0, nNR=0, nI=0
    )
        new(name, A, A_ci, g, g_ci, mean, slope, freq, lind, lat, reftime, opt, PE, SNR, diagn, weights, rf, nR, nNR, nI)
    end
end

mutable struct VectorCoef
    name::Vector{String}
    Lsmaj::Vector{Float64}
    Lsmaj_ci::Vector{Float64}
    Lsmin::Vector{Float64}
    Lsmin_ci::Vector{Float64}
    theta::Vector{Float64}
    theta_ci::Vector{Float64}
    g::Vector{Float64}
    g_ci::Vector{Float64}
    umean::Float64
    vmean::Float64
    uslope::Float64
    vslope::Float64
    freq::Vector{Float64}
    lind::Vector{Int}
    lat::Float64
    reftime::Float64
    opt::Opt
    PE::Vector{Float64}
    SNR::Vector{Float64}
    diagn::Dict{Symbol, Any}
    weights::Vector{Float64}
    rf::Any
    nR::Int
    nNR::Int
    nI::Int

    function VectorCoef(;
        name=String[], Lsmaj=Float64[], Lsmaj_ci=Float64[], Lsmin=Float64[], Lsmin_ci=Float64[],
        theta=Float64[], theta_ci=Float64[], g=Float64[], g_ci=Float64[],
        umean=0.0, vmean=0.0, uslope=0.0, vslope=0.0, freq=Float64[], lind=Int[], lat=0.0, reftime=0.0,
        opt=Opt(), PE=Float64[], SNR=Float64[], diagn=Dict{Symbol, Any}(),
        weights=Float64[], rf=nothing, nR=0, nNR=0, nI=0
    )
        new(name, Lsmaj, Lsmaj_ci, Lsmin, Lsmin_ci, theta, theta_ci, g, g_ci, umean, vmean, uslope, vslope,
            freq, lind, lat, reftime, opt, PE, SNR, diagn, weights, rf, nR, nNR, nI)
    end
end

function Base.getproperty(c::ScalarCoef, s::Symbol)
    if s === :aux
        return (opt = getfield(c, :opt), lat = getfield(c, :lat), reftime = getfield(c, :reftime), frq = getfield(c, :freq), lind = getfield(c, :lind))
    else
        return getfield(c, s)
    end
end

function Base.getproperty(c::VectorCoef, s::Symbol)
    if s === :aux
        return (opt = getfield(c, :opt), lat = getfield(c, :lat), reftime = getfield(c, :reftime), frq = getfield(c, :freq), lind = getfield(c, :lind))
    else
        return getfield(c, s)
    end
end

Base.haskey(c::ScalarCoef, s::Symbol) = hasfield(ScalarCoef, s)
Base.haskey(c::VectorCoef, s::Symbol) = hasfield(VectorCoef, s)
Base.getindex(c::Union{ScalarCoef, VectorCoef}, s::Symbol) = getproperty(c, s)
Base.setindex!(c::Union{ScalarCoef, VectorCoef}, v, s::Symbol) = setproperty!(c, s, v)

mutable struct ReconstructOutput
    t_in::Any
    epoch::Any
    constit::Any
    min_SNR::Float64
    min_PE::Float64
    t_mpl::Vector{Float64}
    h::Union{Nothing, Vector{Float64}}
    u::Union{Nothing, Vector{Float64}}
    v::Union{Nothing, Vector{Float64}}

    function ReconstructOutput(;
        t_in=nothing, epoch=nothing, constit=nothing, min_SNR=2.0, min_PE=0.0,
        t_mpl=Float64[], h=nothing, u=nothing, v=nothing
    )
        new(t_in, epoch, constit, min_SNR, min_PE, t_mpl, h, u, v)
    end
end

function Base.getproperty(r::ReconstructOutput, s::Symbol)
    return getfield(r, s)
end
function Base.setproperty!(r::ReconstructOutput, s::Symbol, v)
    setfield!(r, s, v)
end
Base.haskey(r::ReconstructOutput, s::Symbol) = hasfield(ReconstructOutput, s)
Base.getindex(r::ReconstructOutput, s::Symbol) = getproperty(r, s)
Base.setindex!(r::ReconstructOutput, v, s::Symbol) = setproperty!(r, s, v)

mutable struct Bunch <: AbstractDict{Symbol, Any}
    data::OrderedDict{Symbol, Any}
end
Bunch(; kwargs...) = Bunch(OrderedDict{Symbol, Any}(kwargs))
Bunch(d::AbstractDict) = Bunch(OrderedDict{Symbol, Any}(Symbol(k) => v for (k, v) in d))
Base.getproperty(b::Bunch, s::Symbol) = s === :data ? getfield(b, :data) : b.data[s]
Base.setproperty!(b::Bunch, s::Symbol, v) = s === :data ? setfield!(b, :data, v) : (b.data[s] = v)
Base.getindex(b::Bunch, s::Symbol) = b.data[s]
Base.setindex!(b::Bunch, v, s::Symbol) = (b.data[s] = v)
Base.keys(b::Bunch) = keys(b.data)
Base.iterate(b::Bunch) = iterate(b.data)
Base.iterate(b::Bunch, state) = iterate(b.data, state)
Base.length(b::Bunch) = length(b.data)
Base.haskey(b::Bunch, s::Symbol) = haskey(b.data, s)
Base.get(b::Bunch, s::Symbol, default) = get(b.data, s, default)

function complex_interp(x, xp, fp)
    if eltype(fp) <: Complex
        re = linear_interpolation(xp, real(fp), extrapolation_bc=Flat())(x)
        img = linear_interpolation(xp, imag(fp), extrapolation_bc=Flat())(x)
        return re .+ 1im .* img
    else
        return linear_interpolation(xp, fp, extrapolation_bc=Flat())(x)
    end
end

function _crunch(arr; masked=true)
    if arr isa AbstractArray
        if length(arr) == 1 && ndims(arr) == 0
            return arr[]
        end
        if sum(size(arr) .> 1) <= 1
            return vec(arr)
        end
    end
    return arr
end

function _structured_to_bunch(arr; masked=true)
    if arr isa String return arr end
    if arr isa AbstractDict
        b = Bunch()
        for (k, v) in arr
            b[Symbol(k)] = _structured_to_bunch(v, masked=masked)
        end
        return b
    elseif arr isa AbstractArray && !isempty(arr) && all(x -> x isa AbstractDict, arr)
        if length(arr) == 1
            return _structured_to_bunch(arr[1], masked=masked)
        else
            return [_structured_to_bunch(x, masked=masked) for x in arr]
        end
    end
    return _crunch(arr, masked=masked)
end

function loadbunch(fname; masked=true)
    if endswith(fname, ".mat")
        xx = matread(fname)
    elseif endswith(fname, ".npz")
        #xx = npzread(fname)
    else
        error("Unrecognized file $fname")
    end
    out = Bunch()
    for (k, v) in xx
        if !startswith(k, "__")
            out[Symbol(k)] = _structured_to_bunch(v, masked=masked)
        end
    end
    return out
end

function convert_unicode_arrays(b::Bunch)
    out = Bunch()
    for (key, val) in b
        if val isa AbstractArray && eltype(val) <: AbstractString
            newval = rstrip.(val)
            if length(newval) == 1 newval = newval[1] end
            out[key] = newval
        elseif val isa Bunch
            out[key] = convert_unicode_arrays(val)
        elseif val isa AbstractDict
            out[key] = convert_unicode_arrays(Bunch(val))
        else
            out[key] = val
        end
    end
    return out
end
end # module
