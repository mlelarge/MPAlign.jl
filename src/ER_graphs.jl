"""
    degrees(G::SimpleGraph)

Compute the sequence of degrees
"""
function degrees(G::SimpleGraph)
    degrees(CombinatorialAdjacency(adjacency_matrix(G)))
end

"""
    neighbors_edges(G::SimpleGraph)

Returns deg the sequence of degrees and neig a list of list of direcetd edges numbered from 1 to 2*ne(G)
neig[i] is the list of edges j->i
"""
function neighbors_edges(G::SimpleGraph)
    deg = degrees(G)
    neig = [Vector{Int64}() for _ in vertices(G)]
    c = 1
    for (i,d) in enumerate(deg)
        neig[i] = c:(c+d-1)
        c +=d
    end
    deg, neig
end

"""
    Network

`deg` is a `2*ne(G)` vectors with the degree of the j for each edge `j->i`
`edges_arr` is a `2*ne(G)x2` array containing all directed edges `j->i` where `i` is increasing.
`neig` is a list of vectors where `neig[e]` are all indices of edges `u->j` for `u!=i` where edge `e`
 is `j->i` 
"""
struct Network
    G::SimpleGraph
    deg::Vector{Int64}
    edges_arr::Array{Int64,2}
    neig::Vector{Vector{Int64}}
end

"""
    create_edges(G::SimpleGraph)

returns a 2*ne(G)x2 array containing all directed edges j->i where i is increasing.
"""
function create_edges(G::SimpleGraph)
    edges_array = Array{Int64,2}(undef,2*ne(G),2)
    c = 1
    for i in vertices(G)
        for j in G.fadjlist[i]
            edges_array[c,:] = [j i]
            c +=1
        end
    end
    edges_array
end

"""
    neig_remove(G::SimpleGraph,i::Int64,j::Int64)

returns the neighbors of i in G whitout j
"""
function neig_remove(G::SimpleGraph,i::Int64,j::Int64)
    neig = neighbors(G,i)
    if j ∉ neig
        error(j, "is not a neighbor of ", i, "in ", neig)
    end
    filter(s -> s!=j, neig)
end

"""
    create_dic(edges_arr::Array{Int64,2})

returns a dic such that dic[(i,j)] = index of edge i->j in edges_arr
"""
function create_dic(edges_arr::Array{Int64,2})
    ne, _ = size(edges_arr)
    local_dic = Dict{Tuple,Int64}()
    for i in 1:ne
        local_dic[(edges_arr[i,1],edges_arr[i,2])] = i
    end
    local_dic
end

"""
    create_neig_r(G::SimpleGraph,edges_arr::Array{Int64,2})

returns a list of vectors where neig_r[i] are all indices of edges u->v for u!=w where edge i is v->w
"""
function create_neig_r(G::SimpleGraph,edges_arr::Array{Int64,2})
    ne, _ = size(edges_arr)
    neig_r = [Vector{Int64}() for _ in 1:ne]
    local_dic = create_dic(edges_arr)
    for i in 1:ne
        neig_r[i] = map(s -> local_dic[(s,edges_arr[i,1])], neig_remove(G,edges_arr[i,1],edges_arr[i,2]))
    end
    neig_r
end

"""
    create_degrees(G::SimpleGraph,edges_arr::Array{Int64,2})

retunrs a list where degrees_G[i] is the degree of u where edge i is u->v
"""
function create_degrees(G::SimpleGraph,edges_arr::Array{Int64,2})
    ne, _ = size(edges_arr)
    degrees_G = zeros(ne)
    for i in 1:ne
        degrees_G[i] = degree(G,edges_arr[i,1])
    end
    degrees_G
end

function Network(G::SimpleGraph)
    edges_arr = create_edges(G)
    neig = create_neig_r(G,edges_arr)
    deg = create_degrees(G,edges_arr)
    return Network(G, deg, edges_arr, neig)
end

"""
    Pair_ER

for `i` in `G1`, `perm_G1G2[i]` is the corresponding vertex in `G2`. It is `0` if no corresponding vertex in `G2`.
`edge_G1G2` used for init of bias message passing (TBC).
"""
struct Pair_ER
    N1::Network
    N2::Network
    Ginter::SimpleGraph
    n::Int64
    λ::Float64
    s::Float64
    perm_G1G2::Array{Int64,1}
    perm_G2G1::Array{Int64,1}
    edge_G1G2::Array{Int64,2}
    dic_comb::Dict{Tuple{Int64, Int64}, Vector{Vector{Int64}}}
    dic_perm::Dict{Int64, Vector{Vector{Int64}}}
end

"""
    compute_ab(n,λ,s)

returns the solution of the system of equations:
a+(1-a)*b^2 = λ*s/n and (1-a)*b*(1-b) = λ*(1-s)/n
"""
function compute_ab(n,λ,s)
    λ*(s-λ/n)/(λ*(s-2)+n), λ*(1-s)/(n-λ)
end

"""
    create_cor_graphs(n,λ,s)

returns correlated ER graphs G1, G2 and their intersection. All graphs are aligned.
"""
function create_cor_graphs(n,λ,s)
    a,b = compute_ab(n,λ,s) 
    G_inter = erdos_renyi(n,a)
    G_add1 = erdos_renyi(n,b)
    G_add2 = erdos_renyi(n,b)
    G_1 = union(G_inter,G_add1)
    G_2 = union(G_inter,G_add2)
    return G_1, G_2, intersect(G_1,G_2)
end

function create_perm(n, vmap1, vmap2)
    revert_vmap2 = zeros(Int64,n)
    for (i,v) in enumerate(vmap2)
        revert_vmap2[v] = i
    end
    [revert_vmap2[i] for i in vmap1]
end

function Pair_ER(n,λ,s)
    G1,G2,Ginter = create_cor_graphs(n,λ,s)
    GG1, vmap1 = largest_comp(G1)
    GG2, vmap2 = largest_comp(G2)
    GGinter, vmapinter = largest_comp(Ginter)
    perm_GinterG1 = create_perm(n, vmapinter, vmap1)
    perm_GinterG2 = create_perm(n, vmapinter, vmap2)
    perm_G1G2 = create_perm(n, vmap1, vmap2)
    perm_G2G1 = create_perm(n, vmap2, vmap1)
    N1 = Network(GG1)
    N2 = Network(GG2)
    dic_G1 = create_dic(N1.edges_arr)
    dic_G2 = create_dic(N2.edges_arr)
    Ninter= Network(GGinter)
    ne , _ = size(Ninter.edges_arr)
    edge_G1G2 = Array{Int64,2}(undef,ne,2)
    for (i,e) in enumerate(eachrow(Ninter.edges_arr))
        #println(i,e)
        edge_G1G2[i,:] = [dic_G1[(perm_GinterG1[e[1]],perm_GinterG1[e[2]])], dic_G2[(perm_GinterG2[e[1]],perm_GinterG2[e[2]])]]
    end
    deg1 = unique(sort(N1.deg))
    deg2 = unique(sort(N2.deg))
    min_deg = min(deg1[end],deg2[end])
    max_deg = max(deg1[end],deg2[end])

    dic_comb = Dict([((d,k),collect(combinations(range(1,d),k))) for d in 0:max_deg for k in 0:d])
    dic_perm = Dict([(k,collect(permutations(range(1,k)))) for k in 0:min_deg])
    return Pair_ER(N1,N2,GGinter,n,λ,s,perm_G1G2,perm_G2G1,edge_G1G2,dic_comb,dic_perm)
end


"""
    largest_comp(G::SimpleGraph)

returns the largest connected component of `G` and 
for i a vertex of this component `vmap[i]` is the corresponding index in the original `G`.
"""
function largest_comp(G::SimpleGraph)
    cc = connected_components(G)
    i = argmax([length(c) for c in cc])
    sg, vmap = induced_subgraph(G,cc[i])
    return sg, vmap
end
