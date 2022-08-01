abstract type AbstractQuasiPeriodicFunction2D <: AbstractGreensFunction2D end
angle(G::AbstractQuasiPeriodicFunction2D) = G.theta
period(G::AbstractQuasiPeriodicFunction2D) = G.period
freespace(G::AbstractQuasiPeriodicFunction2D) = FreeSpaceGreensFunction2D(wavenumber(G))

const QuasiPeriodic = AbstractQuasiPeriodicFunction2D

"""
    G = NaiveGreensFunction2D(wavenumber,period,theta)

    G(x,ξ;M)

Computes 2D quasi-periodic Greens function for Helmholtz equation with

    k = wavenumber(G)

that is quasi-periodic with f([x,y],[ξ,η]) = exp(im*k*d*cos(θ))*f([ξ+d,η]) where

    d = period(G)
    θ = theta(G)

by summation of the free space Greens Function for 2M+1 closest periods

"""
Base.@kwdef struct NaiveGreensFunction2D{T,S<:Real} <: AbstractQuasiPeriodicFunction2D
    wavenumber::T
    period::S
    theta::S
    function NaiveGreensFunction2D(wavenumber::T, period, theta) where {T}
        p, a = promote(period, theta)
        return new{T,typeof(p)}(wavenumber, p, a)
    end
end


"""
    G = TailGreensFunction2D(wavenumber,period,theta)

    G(x,ξ;M)

Computes 2D quasi-periodic Greens function for Helmholtz equation with

    k = wavenumber(G)

that is quasi-periodic with f(x,y) = exp(im*x*cos(θ))*f(x+d,y) where

    d = period(G)
    θ = angle(G)

by summation of the free space Greens Function for 2M-1 closest periods, and approximating the tail asymptotically to order M^{-5/2}

"""
Base.@kwdef struct TailGreensFunction2D{T,S<:Real} <: AbstractQuasiPeriodicFunction2D
    wavenumber::T
    period::S
    theta::S
    function TailGreensFunction2D(wavenumber::T, period, theta) where {T}
        p, a = promote(period, theta)
        return new{T,typeof(p)}(wavenumber, p, a)
    end
end

function evaluate(G::QuasiPeriodic, args...; M::Int=128)
    return evaluate(freespace(G),args...) + evaluate_tail(G, args...;M)
end

function evaluate(G::QuasiPeriodic, x,y; M::Int=128)
    return evaluate(freespace(G),x,y) + evaluate_tail(G, x,y;M)
end

function evaluate(G::QuasiPeriodic, x,y,dx; M::Int=128)
    return evaluate(freespace(G),x,y,dx) + evaluate_tail(G, x,y,dx;M)
end

function evaluate(G::QuasiPeriodic, x,y,dx,dy; M::Int=128)
    return evaluate(freespace(G),x,y,dx,dy) + evaluate_tail(G, x,y,dx,dy;M)
end


function evaluate(G::QuasiPeriodic,order::SymOrInt, args...; M::Int=128)
    return evaluate(freespace(G), order, args...) + evaluate_tail(G,args...;M)
end

function evaluate_tail(G::QuasiPeriodic, x, args...;M::Int=128)
    @unpack wavenumber, period, theta = G

    k, d, θ = wavenumber, period, theta
    kdc = k * d * cos(θ)

    F = FreeSpace(wavenumber)
    ret = sum(cis(m * kdc) * evaluate(F, x .+ SVector(d * m,0), args...) for m in 1:M)
    ret += sum(cis(m * kdc) * evaluate(F, x .+ SVector(d * m,0), args...) for m in -(1:M))
end

function symmetrictailend(G, M , args...)
    return tailend(G,true,M,args...) + tailend(G,false,M,args...)
end


# Computes tail-end sum

function tailend(G, s::Bool, M::Int, x::SVec2, y::SVec2)
    @unpack wavenumber, period, theta = G

    k, d, θ = wavenumber, period, theta
    u, v = x-y
    Δ = s ? 1 : -1

    kd, ku, kv = k * d, k * u, k * v
    kd², kv², ku² = kd^2, kv^2, ku^2
    kv⁴ = kv²^2

    α = kd * (1 + Δ * cos(θ))
    eα = cis(α)

    A = -(im + Δ * 4 * ku - 4im * kv²) / 8

    B = ((-9 + 24 * (2 * ku² - kv²)) + Δ * 24im * ku - 16 * kv⁴ - Δ * 96im * ku * kv²) / 128

    C = A / kd - eα / (2 * (1 - eα))

    D = B / kd² - 3 * A * eα / (2 * kd * (1 - eα)) + 3 * eα * (1 + eα) / (8 * (1 - eα)^2)

    λ = -(im / 4) * sqrt(2 / (π * kd)) * cis(Δ * ku - π / 4 + M * α) / (sqrt(M) * (1 - eα))

    return λ * (1 + C / M + D / (M^2))
end


function tailend(G, s::Bool, M::Int, x::SVec2, y::SVec2, dx::SVec2)
    @unpack wavenumber, period, theta = G

    k, d, θ = wavenumber, period, theta
    u, v = x-y
    nu, nv = dx

    Δ = s ? 1 : -1

    kd, ku, kv = k * d, k * u, k * v
    kd², kv², ku² = kd^2, kv^2, ku^2
    kv⁴ = kv²^2

    α = kd * (1 + Δ * cos(θ))
    eα = cis(α)

    A = (3im - Δ * 4 * ku + 4im * kv²) / 8

    B = (15 + 48 * ku² - 56 * kv² - 16 * kv⁴ - Δ * 72im * ku - Δ * 96im * ku * kv²) / 128

    C = (Δ * A * nu + k * nv * v) / kd - Δ * nu * eα / (2 * (1 - eα))

    D = (2 * kv * nv * (A - Δ * ku) + Δ * nu * (2 * B - kv²)) / (2 * kd^2) - 3 * A * eα / (2 * kd * (1 - eα)) - 3eα * (k * nv * v + Δ * A * nu) / (2kd * (1 - eα)) + 3 * Δ * nu * eα * (1 + eα) / (8 * (1 - eα)^2)

    λ = -(im * k / 4) * sqrt(2 / (π * kd)) * cis(Δ * ku - 3π / 4 + M * α) / (sqrt(M) * (1 - eα))

    return -λ * (Δ * nu + C / M + D / (M^2))
end


# TODO: compute asymptotic expansions for neumannkernel
