module MPAlign
using Serialization
using Statistics
using Plots
using BenchmarkTools
using LambertW
using LightGraphs
using LightGraphs.LinAlg, LightGraphs.SimpleGraphs
import LightGraphs.LinAlg: degrees
using LinearAlgebra

include("ER_graphs.jl")
include("utils.jl")
include("message_passing.jl")

export Pair_ER
export run_bp
export create_matrix_lr
export eval_M
export match_edges
export eval_edges
export create_cor_graphs
export largest_comp

include("utils_IO.jl")
export read_graphs
export Pair_IO
export create_matrix_check
export check_data
end # module
