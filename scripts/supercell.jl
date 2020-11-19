ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"
using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/Isosuite")

fn = ARGS[1]
(nx, ny, nz) = (ARGS[2], ARGS[3], ARGS[4])

##

if isfile(fn) && !isdir(fn)
    fnenl = replace(fn,".cif"=>".enlarge($nx,$ny,$nz).cif")
    supercell(
        fn,
        (nx, ny, nz);
        SG_setting = QE_default_equivalent_settings_findsym
    ) ⇶ fnenl
end
