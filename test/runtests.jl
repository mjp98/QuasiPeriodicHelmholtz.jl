using QuasiPeriodicHelmholtz
using SpecialFunctions
using StaticArrays
using Test

import QuasiPeriodicHelmholtz: NaiveGreensFunction2D
import QuasiPeriodicHelmholtz: TailGreensFunction2D
import QuasiPeriodicHelmholtz: FreeSpaceGreensFunction2D
import QuasiPeriodicHelmholtz: hankelh1_1
import QuasiPeriodicHelmholtz: hankelh1_2
import QuasiPeriodicHelmholtz: hinc

import SpecialFunctions:γ, hankelh1
logabs = log∘abs

@testset "QuasiPeriodicHelmholtz.jl" begin
    @testset "series expansions" begin

        @testset "hankelh1(1,x) series" begin
            x = abs(randn(Float64))

            h1x = hankelh1(1,x)

            h1x_split = sum(hankelh1_1(x,n)*cauchyweight2(n,0,x) for n in (0,:log,1,2))

            @test h1x ≈ h1x_split
        end

        @testset "hankelh1(2,x) series" begin
            x = abs(randn(Float64))

            h2x = hankelh1(2,x)
            h2x_split = sum(hankelh1_2(x,n)*cauchyweight2(n,0,x) for n in (0,:log,1,2))
            @test h2x ≈ h2x_split
        end

        @testset "hankelh1(1, x) / x series" begin
            x = abs(randn(Float64))

            h1x = hankelh1(1, x) / x
            h1x_split = sum(hinc(x,n)*cauchyweight2(n,0,x) for n in (0,:log,1,2))

            @test h1x ≈ h1x_split
        end
    end

    @testset "quasi-periodicity" begin
        θ = π/3
        d = 2.3
        k = 1

        x =  @SVector [-0.01,0]
        ξ =  @SVector [0,0]
        offset =  @SVector [d,0]
        @testset "tail" begin
            G = TailGreensFunction2D(;wavenumber=k,period=d,theta=θ)
            @test G(x,ξ;M=1000) ≈ exp(im*k*d*cos(θ))*G(x,ξ-offset;M=1000)
        end
    end

    @testset "published values" begin
        d = 1
        θ = π/4

        x =  @SVector [-0.01,0]
        ξ =  @SVector [0,0]

        @testset "naive" begin
            k = 2
            Gₐ = NaiveGreensFunction2D(;wavenumber=k,period=d,theta=θ)
            @test Gₐ(x,ξ;M=5000) ≈ -0.4634247358 - 0.3530848711im
            k = 10
            Gₐ = NaiveGreensFunction2D(;wavenumber=k,period=d,theta=θ)
            @test Gₐ(x,ξ;M=5000) ≈ -0.3547064117 - 0.1764507543im
        end
        @testset "tail" begin
            k = 2
            Gᵦ = TailGreensFunction2D(;wavenumber=k,period=d,theta=θ)
            @test Gᵦ(x,ξ;M=937) ≈ -0.4595298794 - 0.3509130869im
            k = 10
            Gᵦ = TailGreensFunction2D(;wavenumber=k,period=d,theta=θ)
            @test Gᵦ(x,ξ;M=257) ≈ -0.3538172307 - 0.1769332383im
        end
    end

end
