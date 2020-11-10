using Isosuite
ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/" ##!! change

fn = ARGS[1]
if !isfile(fn) && isdir(fn)
    files = readdir(fn)
    for f in files
        f1 = replace(f,".cif"=>".improved.cif")
        improve_cif("",rstrip(fn,"/")*"/$f") ⇶ (rstrip(fn,"/")*"/$f1")
    end
elseif isfile(fn) && !isdir(fn)
    fn1 = replace(fn,".cif"=>".improved.cif")
    improve_cif("",fn) ⇶ fn1
end

exit()
