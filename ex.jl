using MPAlign

n = 200
λ = 3
s = 0.8
PG = Pair_ER(n,λ,s)

@time mp_fin = run_bp(PG,10,verbose=true)

M_loglr = create_matrix_lr(PG,mp_fin)

ov1, ov2, v1,v2 = eval_M(PG, M_loglr)
println(ov1,ov2)
