using Revise
using Oscar

b = 3 #fix the total amounts
L = [1 1 -b] #conservation laws
C = [1 -1 -1] #stoichimetric matrix
M = [1 0 2;0 1 1] #exponent matrix

#Compute the positive tropicalization of the linear part 
d, n = size(L)  # Size of L
s, r = size(C)  # Size of C
N = [L zeros(Int64,d, r); zeros(Int64,s, n) C] 
#Compute a basis of ker(N)
K = transpose(kernel(matrix(QQ,N);side=:right)) 
#Tropicalize the linear part
TropL = tropical_linear_space(K)

#Compute the positive tropicalization of the binomial part
K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K,t)
R,x,y = polynomial_ring(K,"x"=>1:n,"y"=>1:r)

Ibin = ideal([yi-prod(x[1:end-1].^M[:,i]) for (i,yi) in enumerate(y)])
binomialShifted = [evaluate(g, vcat([yi*t^rand(1:9) for yi in y],x)) for g in gens(Ibin)]
TropHs = tropical_hypersurface.(binomialShifted,Ref(nu))
TropV = reduce(stable_intersection,vcat(TropHs,[TropL]))
tropicalPoints = [first(vertices(sigma)) for sigma in maximal_polyhedra(TropV)]
rootCount = sum([m for (sigma,m) in maximal_polyhedra_and_multiplicities(TropV) ])

#Compute the linear ideal
Ilin = ideal(C*y)+ideal(L*x)
#Check which intersection points are positive
isPos = [ Oscar.is_initial_positive(Ilin,nu,p) for p in tropicalPoints]
count(isPos)

P= maximal_polyhedra(TropV)[1]