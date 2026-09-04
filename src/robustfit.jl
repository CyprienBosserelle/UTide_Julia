module RobustFit
using ..Utilities, LinearAlgebra, Statistics
export robustfit
function andrews(r) r_abs = max.(sqrt(eps(Float64)), abs.(r)); return (r_abs .< pi) .* sin.(r_abs) ./ r_abs end
function bisquare(r) r_abs = abs.(r); return (r_abs .< 1) .* (1 .- r_abs.^2).^2 end
function cauchy(r) return 1.0 ./ (1.0 .+ abs.(r).^2) end
function fair(r) return 1.0 ./ (1.0 .+ abs.(r)) end
function huber(r) return 1.0 ./ max.(1.0, abs.(r)) end
function logistic(r) r_abs = max.(sqrt(eps(Float32)), abs.(r)); return tanh.(r_abs) ./ r_abs end
function ols_weight(r) return ones(length(r)) end
function talwar(r) return Float64.(abs.(r) .< 1) end
function welsch(r) return exp.(-(abs.(r).^2)) end
const wfuncdict = Dict("andrews"=>andrews,"bisquare"=>bisquare,"cauchy"=>cauchy,"fair"=>fair,"huber"=>huber,"logistic"=>logistic,"ols"=>ols_weight,"talwar"=>talwar,"welsch"=>welsch)
const tune_defaults = Dict("andrews"=>1.339,"bisquare"=>4.685,"cauchy"=>2.385,"fair"=>1.400,"huber"=>1.345,"logistic"=>1.205,"ols"=>1.0,"talwar"=>2.795,"welsch"=>2.985)
function complex_median(x)
    if eltype(x) <: Complex
        # Match numpy's behavior: lexicographical sort
        sorted_x = sort(x, by=v -> (real(v), imag(v)))
        n = length(sorted_x)
        if n % 2 == 1
            return sorted_x[(n+1)÷2]
        else
            return (sorted_x[n÷2] + sorted_x[n÷2+1]) / 2.0
        end
    else
        return median(x)
    end
end
function sigma_hat(x) m = complex_median(x); return median(abs.(x .- m)) / 0.6745 end
function leverage(X) pX = pinv(X); return vec(abs.(sum(X' .* pX, dims=1))) end
function r_normed(R, rfac) return rfac .* R ./ sigma_hat(R) end
function robustfit(X, y; weight_function="bisquare", tune=nothing, rcond=1e-15, tol=0.001, maxit=50)
    if isnothing(tune) tune = tune_defaults[weight_function] end
    _wfunc = wfuncdict[weight_function]
    if ndims(X) == 1 X = reshape(X, :, 1) end
    n, p = size(X)
    lev = leverage(X)
    out = Bunch(weight_function=weight_function, tune=tune, rcond=rcond, tol=tol, maxit=maxit, leverage=lev)
    rfac = 1.0 ./ (tune .* sqrt.(1.0 .- lev))
    oldrmeansq, oldlstsq_b, oldw = nothing, nothing, nothing
    iterations, w = 0, ones(length(y))
    b = nothing
    for i in 1:maxit
        wX, wy = w .* X, w .* y
        b = wX \ wy
        rsumsq = sum(abs2, wy .- wX * b)
        if i == 1 out.ols_b, out.ols_rms_resid = b, sqrt(rsumsq / n) end
        rmeansq = rsumsq / sum(w)
        if !isnothing(oldrmeansq)
            improvement = (oldrmeansq - rmeansq) / oldrmeansq
            if improvement < 0 b, w, iterations = oldlstsq_b, oldw, i - 1; break end
            if improvement < tol iterations = i; break end
        end
        oldlstsq_b, oldw, oldrmeansq = b, w, rmeansq
        resid = y .- X * b
        w = _wfunc(r_normed(resid, rfac))
        iterations = i
    end
    resid = y .- X * b
    out.iterations, out.b, out.w, out.rms_resid = iterations, b, w, sqrt(mean(abs2, resid))
    return out
end
end # module
