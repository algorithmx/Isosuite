using LinearAlgebra

ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"
using Isosuite

#fn = ARGS[1]
fn = "/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/try1/WO3/mp-19443.improved.cif"

##

w(s) = match(r"\s+[a-n]\s+(-)?\d+\.\d+\s+(-)?\d+\.\d+\s+(-)?\d+\.\d+",s).match
w.(extract_atom_config(fn))

WYCKPOS_IT(x::Int) = [k for k in keys(WYCKPOS_TABLE) if occursin("$x",k)]

WYCKPOS_IT(130)