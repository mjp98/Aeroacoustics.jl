module Aeroacoustics

using ApproxFun
using IntervalSets
using HolomorphicFun
using LinearAlgebra
using Parameters
using SpecialFunctions
using StaticArrays
using StatsBase
using Unitful
using UnPack

# Geometric mean
import Base: *,+,-,/
import StatsBase: geomean
import ApproxFun: Fun, coefficients, evaluate, space

export PlaneWave, fourierx
export convectedgust
export reference_density, reference_length, reference_speed, reference_time

export characteristic_length

export SawtoothWave

# export convectedgust
include("util.jl")
include("units.jl")
include("turbulence/energyspectra.jl")


include("incident/planewave.jl")
include("incident/gust.jl")

include("boundary/waveform.jl")





end
