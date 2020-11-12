using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/Isosuite")

using Isosuite
ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/" ##!! change

#

#fn = ARGS[1]
fn = "/home/dabajabaza/jianguoyun/Workspace/ScF3/materials/R-3c/"

##

if !isfile(fn) && isdir(fn)
    fd = rstrip(fn,'/')
    files = [ff for ff in readdir(fn) if endswith(ff,".cif") && !occursin("improved",ff)]
    for f in files
        f1 = replace(f,".cif"=>".improved.cif")
        improve_cif("", "$fd/$f", SG_setting=QE_default_equivalent_settings_findsym) ⇶ "$fd/$f1"
    end
    try mkdir("$fd/improved") catch _ nothing end
    for ff in readdir(fn)
        if endswith(ff,".improved.cif") 
            try run(`mv $fd/$ff $fd/improved/`) catch _ nothing end
        end
    end
elseif isfile(fn) && !isdir(fn)
    fn1 = replace(fn,".cif"=>".improved.cif")
    improve_cif("", fn, SG_setting=QE_default_equivalent_settings_findsym) ⇶ fn1
end

##

exit()
