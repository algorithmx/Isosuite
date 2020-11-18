using LinearAlgebra

ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"
using Isosuite

fn = ARGS[1]
fn = "/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/try1/WO3/mp-510417.improved.cif"
(nx, ny, nz) = (ARGS[2], ARGS[3], ARGS[4])
(nx, ny, nz) = (2, 2 ,2)

##

if isfile(fn) && !isdir(fn)
    fnenl = replace(fn,".cif"=>".enlarge($nx,$ny,$nz).cif")
    supercell(
        fn,
        (nx, ny, nz);
        SG_setting = QE_default_equivalent_settings_findsym
    ) ⇶ fnenl
end
