function read_graphs(n, file) 
    G1 = SimpleGraph(n)
    G2 = SimpleGraph(n)
    const_G2 = 0
    for line in eachline(file)
        #println(line)
        if line == "graph2"
            const_G2 += 1
        elseif line == "matrice log likelihood ratio" 
            const_G2 += 1
        end
        if const_G2 == 0 
            u = split(line, " : ")
            if length(u) >1
                add_edge!(G1, parse(Int, u[1]),parse(Int, u[2]))
                #println(parse(Int, u[1]), "|",  parse(Int, u[2]))
            end
        elseif const_G2 == 1
            u = split(line, " : ")
            if length(u) >1
                add_edge!(G2, parse(Int, u[1]),parse(Int, u[2]))
                #println(parse(Int, u[1]), "|",  parse(Int, u[2]))
            end
        end
    end
    return G1, G2
end

function Pair_IO(G1,G2,n,λ,s)
    N1 = Network(G1)
    N2 = Network(G2)
    #neig1 = create_neig_r(G1)
    #neig2 = create_neig_r(G2)
    deg1 = unique(sort(N1.deg))
    deg2 = unique(sort(N2.deg))
    min_deg = min(deg1[end],deg2[end])
    max_deg = max(deg1[end],deg2[end])

    dic_comb = Dict([((d,k),collect(combinations(range(1,d),k))) for d in 0:max_deg for k in 0:d])
    dic_perm = Dict([(k,collect(permutations(range(1,k)))) for k in 0:min_deg])
    return Pair_ER(N1,N2,SimpleGraph(n),n,λ,s,1:n,1:n,zeros(Int64,n,2),dic_comb,dic_perm)
end

function create_matrix_check(data)
    n = length(data)
    M_data = zeros(n,n)
    for (i,l) in enumerate(data)
        for (j,e) in enumerate(l)
            M_data[i,j] = e
        end
    end
    M_data
end

function check_data(PG,M_data,n_iter)
    n,_ = size(M_data)
    mp_fin = run_bp(PG,n_iter)
    M_loglr = create_matrix_lr(PG,mp_fin)
    check = true
    for i in 1:n
        for j in 1:n
            if (M_data[i,j] - M_loglr[i,j]) > 0.00001
                #println(i,"|",j,"|", M_data[i,j],"|",M_loglr[i,j])
                check = false
            end
        end
    end
    check
end
