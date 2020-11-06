using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/Isosuite")

using Isosuite
ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"
@inline decif(x) = first(SPLTD(x))

fn = ARGS[1]
if !isfile(fn) && ispath(fn)
    files = readdir(fn)
    for f in files
        f1 = replace(f,".cif"=>".improved.cif")
        improve_cif("",rstrip(fn,"/")*"/$f") ⇶ (rstrip(fn,"/")*"/$f1")
    end
elseif isfile(fn) && !ispath(fn)
    fn1 = replace(fn,".cif"=>".improved.cif")
    improve_cif("",fn) ⇶ fn1
end

exit()