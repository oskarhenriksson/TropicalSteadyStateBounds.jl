

module TropicalSteadyStateBounds

    using Oscar
    using Catalyst
    import MixedSubdivisions: mixed_volume
    import ProgressMeter

    include("tropical_geometry.jl")
    include("transversal_presentation.jl")
    include("certify_genericity.jl")
    include("reaction_networks.jl")
    include("augmented_vertical_systems.jl")
    include("generic_root_counts.jl")
    include("positive_bounds.jl")
    include("toric_bounds.jl")
    include("misc.jl")

    function __init__()
        if isdefined(ProgressMeter, :ijulia_behavior)
            ProgressMeter.ijulia_behavior(:clear)
        end
    end

end



