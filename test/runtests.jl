using Test
using Documenter
using VerticalRootCounts


DocMeta.setdocmeta!(VerticalRootCounts, :DocTestSetup,
    :(using VerticalRootCounts, Oscar, Catalyst),
    recursive = true,
)

doctest(VerticalRootCounts, manual = false)

include("test_root_bounds.jl")
include("test_crnt.jl")