# Helper functions
# ψ = digamma
# bessely_coeff(k,n) = (ψ(k + 1) + ψ(k + n + 1)) / (factorial(k) * factorial(n + k))
# besselj_coeff(k,n) = inv(factorial(k) * gamma(n + k + 1))
# function besselj1_asymptotic(z, ::Val{N}) where N
#     n = 1
#     hz = z / 2
#     hz2 = hz^2
#     coeff(k) = convert(ComplexF64, besselj_coeff(big(k),big(n)))
#     ret = @evalpoly -hz2 (coeff(k) for k = 0:N)...
#     return (hz^n) * ret
# end
# function bessely1_asymptotic(z, ::Val{N}) where N
#     n = 1
#     hz = z / 2
#     hz2 = hz^2
#     coeff(k) = convert(ComplexF64, bessely_coeff(big(k),big(n)))
#     ret = @evalpoly -hz2 (coeff(k) for k = 0:N)...
#     # terms change for n !== 1
#     return -inv(hz) / π + (2 / π) * log(hz) * besselj(n, z) - (hz / π) * ret
# end

# besselj(1,x)/x = jinc(x / π) / 2

# hankelh1(1,x)/x

# seems a bit like hinc, so let's call it that...

const hinc_coeffs = (
    -0.15443132980306573,
    0.6727843350984671,
    0.18157516696085563,
    0.019182189839330562,
    0.001115359491966528,
    4.142247689271143e-5,
    1.071545914091181e-6,
    2.045286003593878e-8,
    3.002048746589188e-10,
    3.495928729692882e-12,
    3.309914735250272e-14,
    2.5986411321011286e-16,
    1.7195232826992565e-18,
    9.721207518823618e-21,
    4.750281743327667e-23,
    2.0264937604328578e-25
)
const hankelh1_2_coeffs = (
    0.17278433509846713,
    0.27981700058837794,
    0.05060212507354724,
    0.0041142157456438904,
    0.0001955383103894831,
    6.153702292307226e-6,
    1.382490703901454e-7,
    2.3332924714182106e-9,
    3.0703952724423275e-11,
    3.240877840449121e-13,
    2.8062045674315818e-15,
    2.0299018636752234e-17,
    1.2453360568075733e-19,
    6.562675770946843e-22,
    3.003191194935999e-24,
    1.204702020224618e-26
)

hinc(x) = hankelh1(1, x) / x

function hinc(x, order)
    @assert order in (0, :log, 1, 2)
    order == :log && return 2im * jinc(x / π) / 2
    order == 1 && return 0
    order == 2 && return -2im

    jincx = jinc(x / π) / 2
    x2 = -(x / 2)^2
    jincx - (2im / π) * (log(2)) * jincx - (im / (2π)) * (@evalpoly x2 hinc_coeffs...)
end

hankelh1_1(x) = hankelh1(1, x)

function hankelh1_1(x, order)
    @assert order in (0, :log, 1, 2)
    order == :log && return  2im * besselj(1, x)
    order == 1 && return -2im
    order == 2 && return 0

    jincx = jinc(x / π) / 2
    x2 = -(x / 2)^2
    x * (jincx - (2im / π) * (log(2)) * jincx - (im / (2π)) * (@evalpoly x2 hinc_coeffs...))
end

hankelh1_2(x) = hankelh1(2, x)

function hankelh1_2(x, order)
    @assert order in (0, :log, 1, 2)
    order == :log && return 2im * besselj(2, x)
    order == 1 && return 0
    order == 2 && return -4im

    besselj2x = besselj(2, x)
    x2 = -(x / 2)^2
    besselj2x - (2im / π) * (log(2)) * besselj2x - (im / π) + (x2 * im / π) * (@evalpoly x2 hankelh1_2_coeffs...)
end
