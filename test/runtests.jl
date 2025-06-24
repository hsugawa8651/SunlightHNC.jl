
# This is a Julia port of SunlightHNC 
# Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

using SunlightHNC
using Test

@testset "SunlightHNC.jl" begin

    @testset "Utitilies" begin
        @test_throws ArgumentError upperTriangularIndices(1, 2)
    end

    @testset "Constructot" begin
        @test_throws ArgumentError Wizard(; ncomp=0)
        @test_throws ArgumentError Wizard(; ng=1)
        @test_throws ArgumentError Wizard(; deltar=-1.0e-1)
        @test_throws ArgumentError Wizard(; maxsteps=-1)
        @test_throws ArgumentError Wizard(; tol=-1.0e-1)
        @test_throws ArgumentError Wizard(; nps=0)
        @test_throws ArgumentError Wizard(; maxsteps=0)
        @test_throws ArgumentError Wizard(; start_type=0)
    end

    @testset "fft" begin
        ng = 2^8
        w = Wizard(; ng=ng)
        @test length(w.fftw) == ng - 1
        a = similar(w.fftw)
        @test length(a) == ng - 1
        a .= 1.0e0 ./ sqrt.(1:ng-1)
        b = w.planFFT * a
        b ./= sqrt(2.0e0 * ng)
        #
        c = w.planFFT * b
        c ./= sqrt(2.0e0 * ng)
        #
        d = similar(a)
        for k in 1:ng-1
            d[k] = sum(a .* sinpi.((1:ng-1) .* k / ng))
        end
        d .*= sqrt(2.0e0 / ng)
        @test all(isapprox.(a, c))
        @test all(isapprox.(b, d))
    end

    @testset "dpd_potential" begin
        w = Wizard()
        @test_throws ArgumentError dpd_potential!(w, 11)
    end

    @testset "urpm_potential" begin
        w = Wizard()
        @test_throws ArgumentError urpm_potential!(w)
    end

    @testset "rpm_potential" begin
        w = Wizard()
        @test_throws ArgumentError rpm_potential!(w)
    end
end
