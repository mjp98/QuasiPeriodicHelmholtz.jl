module QuasiPeriodicHelmholtz

using SpecialFunctions
using UnPack
using DualNumbers
using StaticArrays
using LinearAlgebra

import SpecialFunctions: γ
import Base: angle

export cauchyweight2

include("util.jl")
include("hankelseries.jl")
include("greens.jl")
include("freespace.jl")
include("quasiperiodic.jl")
include("waveguide.jl")

end
