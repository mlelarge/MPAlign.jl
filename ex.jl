using MPAlign

n = 100
λ = 3
s = 0.55
PG = Pair_ER(n,λ,s)

@time mp_fin = run_bp(PG,5)

M_lr = create_matrix_lr(PG,mp_fin)
M_loglr = log.(M_lr)

ov1, ov2, v1,v2 = eval_M(PG, M_loglr)
println(ov1,ov2)
