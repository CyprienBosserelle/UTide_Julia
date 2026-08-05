module TimeConversion
using Dates
export normalize_time
const DAY_TO_GREGORIAN_EPOCH = 719163
function date2num(date::DateTime)
    return (Dates.value(date) - Dates.value(DateTime(1970, 1, 1))) / (24 * 60 * 60 * 1000)
end
function date2num(date::Date)
    return Float64(Dates.value(date) - Dates.value(Date(1970, 1, 1)))
end
function python_gregorian_datenum(date::Union{Date, DateTime})
    return date2num(date) + DAY_TO_GREGORIAN_EPOCH
end
function python_gregorian_datenum(t::AbstractArray{<:TimeType})
    return [python_gregorian_datenum(x) for x in t]
end
function normalize_time(t, epoch=nothing)
    if isnothing(epoch)
        if eltype(t) <: TimeType || t isa TimeType
            return python_gregorian_datenum(t)
        elseif eltype(t) <: Number
            return Float64.(t)
        else
            return python_gregorian_datenum(t)
        end
    end
    if eltype(t) <: Number
        if epoch == "python" return Float64.(t)
        elseif epoch == "matlab" return Float64.(t .- 366)
        else
            ofs = python_gregorian_datenum(epoch)
            return Float64.(t .+ ofs)
        end
    else
        error("Can not process time array as timestamp or datenum.")
    end
end
end # module
