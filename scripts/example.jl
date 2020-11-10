using Isosuite
ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/" ##!! change
@inline decif(x) = first(SPLTD(x))

fn = ARGS[1]
println(ispath(fn))
println(isfile(fn))

if !isfile(fn) && ispath(fn)
    files = readdir(fn)
    for f in files
        f1 = replace(f,".cif"=>".improved.cif")
        old_fn = rstrip(fn,"/")*"/$f"
        improved = improve_cif("",old_fn)
        new_fn = rstrip(fn,"/")*"/$f1"
        write_to_file(join(improved,"\n"), new_fn)
    end
elseif isfile(fn) && !ispath(fn)
    fn1 = replace(fn,".cif"=>".improved.cif")
    improved = improve_cif("",fn)
    println(length(improved))
    write_to_file(join(improved,"\n"), fn1)
end

exit()