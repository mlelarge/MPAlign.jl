#include("ER_graphs.jl")

using Memoize

@memoize function psi(PG::Pair_ER,k,d1,d2)
    exp(PG.λ*PG.s)*PG.s^k*(1-PG.s)^(d1+d2-2*k)/(PG.λ^k)
end

@memoize function psi(λ,s,k,d1,d2)
    exp(λ*s)*s^k*(1-s)^(d1+d2-2*k)/(λ^k)
end

@memoize function test_psi(k::Int,d1::Int,d2::Int)
    deg_min = min(d1,d2)+1
    factorial(d1-k)*factorial(d2-k)*factorial(k)/(factorial(d1)*factorial(d2)*deg_min)
end

function overlap(PG,max_index1,max_index2)
    sum([PG.perm_G2G1[i[2]]==i[1] for i in max_index1]), sum([PG.perm_G1G2[i[1]]==i[2] for i in max_index2])
end

"""
    eval_M(PG,M)

Inputs: a pair of networks `PG` and a matrix `M` of scores where M_{iu} is 
the score of pairing vertex i in G1 to vertex u in G2.
Creates the mapping from vertices in G1 to vertices in G2 by taking
the argmax on each row of `M`, then compute the overlap with the true permutation
(i.e. count the number of vertices in G1 matched correctly).
Returns also the sum of the max scores on each row.
Does a similar computation on the columns of `M`.
"""
function eval_M(PG,M)
    vmax2, max_index2 = findmax(M,dims=2)
    vmax1, max_index1 = findmax(M,dims=1)
    ov1, ov2 = overlap(PG,max_index1, max_index2)
    ov1, ov2, sum(vmax1), sum(vmax2)
end

"""
    match_edges(PG,map1,map2)

Inputs: a pair of networks `PG`, `map1` and `map2` such that map1[i] = u
for a vertex i in G1 and u in G2 and similarly map2[v] = j. 
"""
function match_edges(PG,map1,map2)
    sum(has_edge(PG.N2.G,map1[src(e)],map1[dst(e)]) for e in edges(PG.N1.G)), sum(has_edge(PG.N1.G,map2[src(e)],map2[dst(e)]) for e in edges(PG.N2.G))
end

function eval_edges(PG,M)
    _, max_index2 = findmax(M,dims=2)
    _, max_index1 = findmax(M,dims=1)
    map1 = [c[2] for c in max_index2]
    map2 = [c[1] for c in max_index1]
    match_edges(PG,map1,map2)
end