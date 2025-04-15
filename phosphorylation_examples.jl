include("src/TropicalBounds.jl")

m = 3

A = matrix(ZZ, [transpose(collect(1:m))  transpose(collect(1:m))  transpose(collect(1:m)) 1 0 0;
    -transpose(collect(1:m))  -transpose(collect(0:m-1))  -transpose(collect(0:m-1)) 0 1 0;
    ones(Int,1,m) ones(Int,1,m) ones(Int,1,m)§ 0 0 1]
)

L = matrix(QQ, [ones(Int, 1, m) ones(Int, 1, m) ones(Int, 1, m) 0 0 1;
    zeros(Int, 1, m) ones(Int, 1, m) zeros(Int, 1, m) 1 0 0;
    zeros(Int, 1, m) zeros(Int, 1, m) ones(Int, 1, m) 0 1 0]
)

# This gives the correct answer
@time toric_root_bound(A, L, check_transversality=false)

# This gives a too low answer
@time toric_root_bound(A, L, check_transversality=true)