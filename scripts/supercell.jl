using LinearAlgebra

ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"
using Isosuite


fn = ARGS[1]
#fn = "/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/try1/WO3/mp-756478.improved.cif"
(nx, ny, nz) = (ARGS[2], ARGS[3], ARGS[4])
#(nx, ny, nz) = (2, 4 ,2)

##

if isfile(fn) && !isdir(fn)
    cif1 = supercell(
        fn,
        (nx, ny, nz);
        SG_setting = QE_default_equivalent_settings_findsym
    )
    fnenl = replace(fn,".cif"=>".enlarge($nx,$ny,$nz).cif")
    cif1 ⇶ fnenl
    fnimp = replace(fn,".cif"=>".enlarge($nx,$ny,$nz).improved.cif")
    improve_cif("", fnenl, SG_setting=QE_default_equivalent_settings_findsym) ⇶ fnimp
end
