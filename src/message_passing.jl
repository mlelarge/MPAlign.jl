using Combinatorics
using Serialization

function init_mp(PG::Pair_ER)
    zeros(2*ne(PG.N1.G),2*ne(PG.N2.G))
end


function sum_combi(neig_i1::Array{Int64,1},neig_i2::Array{Int64,1},k::Int64,
    mp_previous::Array{Float64,2},gen_i1::Vector{Vector{Int64}},gen_i2::Vector{Vector{Int64}},perm_k::Vector{Vector{Int64}})
    if  k == 0
        return Float64(0)
    end
    sum_inj = Float64(0)
    for s1 in gen_i1
        for s2 in gen_i2
            for ps2 in perm_k
                sum_inj += prod(exp(mp_previous[neig_i1[s1[j]],neig_i2[s2[ps2[j]]]]) for j in 1:k)
            end
        end
    end
    log(sum_inj)
end


function m_passing(PG::Pair_ER,e1::Int64,e2::Int64,mp_previous::Array{Float64,2})
    # trees T_1(e1) and T_2(e2)
    dic_comb = PG.dic_comb
    dic_perm = PG.dic_perm
    d1 = PG.N1.deg[e1]
    d2 = PG.N2.deg[e2]
    deg_min = min(d1,d2)
    if deg_min == 1
        output = psi_log(PG,0,d1-1,d2-1)
    else
        neig_i1 = PG.N1.neig[e1]
        neig_i2 = PG.N2.neig[e2]
        gen_i1,gen_i2,perm_k = dic_comb[(d1-1,0)], dic_comb[(d2-1,0)], dic_perm[0]
        output = psi_log(PG,0,d1-1,d2-1) + sum_combi(neig_i1,neig_i2,0,mp_previous, gen_i1, gen_i2, perm_k)#Float64(0)
        for k in 1:(deg_min-1)
            gen_i1,gen_i2,perm_k = dic_comb[(d1-1,k)], dic_comb[(d2-1,k)], dic_perm[k]
            output = log_sum_exp(output, psi_log(PG,k,d1-1,d2-1) + sum_combi(neig_i1,neig_i2,k,mp_previous, gen_i1, gen_i2, perm_k))
        end
    end
    output
end

@inline function compute_next!(PG::Pair_ER,mp_previous::Array{Float64,2},mp_next::Array{Float64,2},ne1::Int64,ne2::Int64)
    @inbounds for e1 in 1:ne1
        @inbounds for e2 in 1:ne2
            mp_next[e1,e2] = m_passing(PG,e1,e2,mp_previous)
        end    
    end
    nothing
end


function run_bp(PG::Pair_ER,n_iter::Int64;mp_previous=nothing,verbose=false)
    if mp_previous === nothing
        mp_previous = init_mp(PG)
    end
    mp_next = similar(mp_previous)
    mp_opt = similar(mp_previous)
    (ne1,ne2) = size(mp_previous)
    res = zeros(6,n_iter)
    perf = 0
    last_iter = 0
    ov_opt = 0
    for i in 1:n_iter
        compute_next!(PG,mp_previous,mp_next,ne1,ne2)
        (mp_previous, mp_next) = (mp_next, mp_previous)
        M_loglr = create_matrix_lr(PG,mp_previous)
        #M_loglr = log.(M_loglr)
        ov1, ov2, v1,v2 = eval_M(PG, M_loglr)
        m1,m2 = eval_edges(PG, M_loglr)
        new_perf = (m1+m2)/2
        if new_perf > perf
            perf = new_perf
            mp_opt = deepcopy(mp_previous)
            last_iter = i
            ov_opt = (ov1+ov2)/2
        end
        res[1,i] = ov1
        res[2,i] = ov2
        res[3,i] = v1
        res[4,i] = v2
        res[5,i] = m1
        res[6,i] = m2
        if verbose
            println(i, " | ", m1, " | ", m2, " | " ,ov1," | ",ov2, " | ", v1, " | ", v2)
        end
        if max(v1,v2) == Inf
            break
        end
    end
    mp_opt#res, mp_opt, last_iter, perf, ov_opt/PG.n
end

function lr_root(PG::Pair_ER,mp,d1,d2,neig_i1,neig_i2)
    dic_comb = PG.dic_comb
    dic_perm = PG.dic_perm
    deg_min = min(d1,d2)
    gen_i1,gen_i2,perm_k = dic_comb[(d1,0)], dic_comb[(d2,0)], dic_perm[0]
    sum = psi_log(PG,0,d1,d2) + sum_combi(neig_i1,neig_i2,0,mp, gen_i1, gen_i2, perm_k)
    for k in 1:deg_min
        gen_i1,gen_i2,perm_k = dic_comb[(d1,k)], dic_comb[(d2,k)], dic_perm[k]
        sum = log_sum_exp(sum, psi_log(PG,k,d1,d2) + sum_combi(neig_i1,neig_i2,k,mp, gen_i1, gen_i2, perm_k))
    end
    sum
end

function create_matrix_lr(PG::Pair_ER,mp)
    v1 = vertices(PG.N1.G)
    v2 = vertices(PG.N2.G)
    matrix_lr = zeros(length(v1),length(v2))
    deg1, neig1 = neighbors_edges(PG.N1.G)
    deg2, neig2 = neighbors_edges(PG.N2.G)
    for i in v1
        for j in v2
            matrix_lr[i,j] = lr_root(PG,mp,deg1[i],deg2[j],neig1[i],neig2[j])
        end
    end
    matrix_lr
end

function init_mp_biased(PG::Pair_ER,bias=10)
    mp_init = ones(2*ne(PG.N1.G),2*ne(PG.N2.G))
    for (i,j) in eachrow(PG.edge_G1G2)
        mp_init[i,j] = bias
    end
    mp_init
end