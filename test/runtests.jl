using Test
using Documenter
using TropicalSteadyStateBounds


DocMeta.setdocmeta!(TropicalSteadyStateBounds, :DocTestSetup,
    :(using TropicalSteadyStateBounds, Oscar, Catalyst),
    recursive = true,
)

doctest(TropicalSteadyStateBounds, manual = false)

include("test_root_bounds.jl")
include("test_crnt.jl")