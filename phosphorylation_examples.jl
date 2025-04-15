include("src/TropicalBounds.jl")

# Number of phosphorylation sites
m = 3

# Exponent matrix
A = matrix(ZZ, [transpose(collect(1:m))  transpose(collect(1:m))  transpose(collect(1:m)) 1 0 0;
    -transpose(collect(1:m))  -transpose(collect(0:m-1))  -transpose(collect(0:m-1)) 0 1 0;
    ones(Int,1,m) ones(Int,1,m) ones(Int,1,m)§ 0 0 1]
)

# Conserved quantities (rows correspond to S_tot, E_tot, F_tot)
L = matrix(QQ, [ones(Int, 1, m) ones(Int, 1, m) ones(Int, 1, m) 0 0 1;
    zeros(Int, 1, m) ones(Int, 1, m) zeros(Int, 1, m) 1 0 0;
    zeros(Int, 1, m) zeros(Int, 1, m) ones(Int, 1, m) 0 1 0]
)

# Toric steady state bound
@time toric_root_bound(A, L, check_transversality=false)

# This gives a too low answer!
@time toric_root_bound(A, L, check_transversality=true)