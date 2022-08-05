import QuasiPeriodicHelmholtz: TailGreensFunction2D, WaveguideNeumann, WaveguideDirichlet

using StaticArrays
d = 2
k = 4

GN = WaveguideNeumann(TailGreensFunction2D(k, 2d, π / 2), d)
GD = WaveguideDirichlet(TailGreensFunction2D(k, 2d, π / 2), d)

xx = LinRange(-3d, 3d, 400)
yy = LinRange(-3d, 3d, 400)

ξ = SVector(0.2, 1.2)

using ProgressMeter
zn = @showprogress [GN(SVector(x, y), ξ; M=64) for x in xx, y in yy]
zd = @showprogress [GD(SVector(x, y), ξ; M=64) for x in xx, y in yy]


using Plots
begin
    heatmap(xx, yy, real.(zn); c=:jet)
    plot!(collect(extrema(xx)), [d, d]; c=:black, label=:none)
    plot!(collect(extrema(xx)), [0, 0]; c=:black, label=:none)
    ymin, ymax = extrema(yy)
    for t in ((ymin+mod(ymin, d)):d:(ymax+mod(ymax, d)-1))
        plot!(collect(extrema(xx)), [t, t]; c=:black, ls=:dash, label=:none)
    end
    plot!(xlims=extrema(xx), ylims=extrema(yy))
    plot!(aspectratio=1)
end

using BenchmarkTools

@benchmark QuasiPeriodicHelmholtz.evaluate($G, SVector(complex(2.0), 3), SVector(4, 5.0); M=64)
