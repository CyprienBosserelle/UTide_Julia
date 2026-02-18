module Utilities
using LinearAlgebra, Statistics, NPZ, OrderedCollections, Interpolations, MAT

export Bunch, complex_interp, loadbunch, convert_unicode_arrays

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
        xx = npzread(fname)
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
