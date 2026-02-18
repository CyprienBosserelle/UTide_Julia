using UTide
using Test
using Statistics

@testset "UTide.jl" begin
    @testset "Scalar Solve/Reconstruct" begin
        t = collect(0:1:1000) ./ 24.0
        freq_M2 = 1.0/12.42
        u = 2.0 .* cos.(2.0 .* pi .* (t .* 24.0) .* freq_M2)
        lat = 45.0
        coef = solve(t, u, lat=lat, phase="raw", nodal=false)
        tide = reconstruct(t, coef)
        @test sqrt(mean(abs2, tide.h .- u)) < 0.05
        @test any(coef.name .== "M2")
    end

    @testset "Vector Solve/Reconstruct" begin
        t = collect(0:1:1000) ./ 24.0
        freq_M2 = 1.0/12.42
        wt = 2.0 .* pi .* (t .* 24.0) .* freq_M2
        u = 2.0 .* cos(30*pi/180) .* cos.(wt) .- 0.5 .* sin(30*pi/180) .* sin.(wt)
        v = 2.0 .* sin(30*pi/180) .* cos.(wt) .+ 0.5 .* cos(30*pi/180) .* sin.(wt)
        lat = 45.0
        coef = solve(t, u, v, lat=lat, phase="raw", nodal=false)
        tide = reconstruct(t, coef)
        @test sqrt(mean(abs2, tide.u .- u)) < 0.05
        @test sqrt(mean(abs2, tide.v .- v)) < 0.05
    end
end
