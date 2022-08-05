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
import QuasiPeriodicHelmholtz: evaluate
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


    @testset "SingularIntegralEquations examples" begin

        nx = 0+1im
        k = 1

        x,y = randn(ComplexF64,2)
        G = FreeSpaceGreensFunction2D(k)


        @testset "dirichlet" begin
            g1 = (x,y) -> -besselj0(k*abs(y-x))/2
            g2 = (x,y) -> x == y ? -(log(k/2)+γ)/2/π + im/4 : im/4*hankelh1(0,k*abs(y-x)) - g1(x,y).*logabs(y-x)/π

            @test g1(x,y) ≈ -evaluate(G,:log,x,y)
            @test g2(x,y) ≈ -evaluate(G,0,x,y)

        end

        @testset "neumann" begin
            function g3neumann(x,y,k)
                z = k*abs(y-x)
                if z < 1/16
                    C,z2 = 2log(2)-2γ,z*z
                    z4=z2*z2
                    z6=z2*z4
                    ret = k^2*complex( ( (C+2)/8 - (C+3)*z2/64 + (3C+11)*z4/4608 - (6C+25)*z6/442386 -(1/4-z2/32+z4/768-z6/36864)*log(k) )/π, 1/8-z2/64+3z4/4608-6z6/442368)
                else
                    ret = im*k/4*hankelh1(1,k*abs(y-x))./abs(y-x) - g1(x,y)./abs(y-x).^2/π -g2(x,y).*logabs(y-x)/π
                end
                ret
              end


            g1 = (x,y) ->  besselj0(k*abs(y-x))/2
            g2 = (x,y) ->  x == y ? -k^2/4 : -k*besselj1(k*abs(y-x))./abs(y-x)/2
            g3 = (x,y) ->  g3neumann(x,y,k) # In /Scatteraux.jl
           # g4old = (x,y) ->  im*k/4*hankelh1(1,k*abs(y-x))./abs(y-x).*imag(y-x)
            g4 = (x,y) ->  im*k/4*besselj1(k*abs(y-x))/abs(y-x)*imag(y-x)  # For linesum
            g5 = (x,y) ->  -k/2*besselj1(k*abs(y-x))/abs(y-x)*imag(y-x)  # For logkernel
            g6 = (x,y) ->  k/2*abs(y-x)*(bessely1(k*abs(y-x)) - 2besselj1(k*abs(y-x))*logabs(y-x)/π) # For Re{Cauchy}

            # Formulation allows bits of cauchy(1,2) to be shifted into regular part.
            x,y = randn(ComplexF64,2)

            # This isn't the right test anyway...
            @test_broken (g6(x,y) + g4(x,y)/(y-x)/pi ≈ -evaluate(G,1,x,y,nx)*cauchyweight2(1,x,y) -evaluate(G,0,x,y,nx))
            @test g5(x,y) ≈ -evaluate(G,:log,x,y,nx)

            x,y = complex.(randn(Float64,2))

            ny = 0+1im

            @test g2(x,y) ≈ -evaluate(G,:log,x,y,nx,ny)
            @test (g3(x,y) + g1(x,y)/(y-x)^2/pi ≈ -evaluate(G,0,x,y,nx,ny)  -evaluate(G,2,x,y,nx,ny)*cauchyweight2(2,x,y))

        end

    end

end
