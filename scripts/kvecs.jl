ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"
using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/Isosuite")
using Isosuite

using JLD2

##

#kvecs = [(print("i=$i\n"); all_kvectors(i, unitcell_setting=QE_default_equivalent_settings_iso)) for i =1:230]
#@save "kvecs.jld2" kvecs

@load  "kvecs.jld2"  kvecs

##


IT_number_of_crystal_classes = [1:2, 3:15, 16:74, 75:142, 143:167, 168:194, 195:230]
@inline take0(dic) = [v[1] for (k,v) in dic if v[2]==0] |> sort
KVECS_BY_CRYSTAL_CLASS = [  union(unique([take0(kvecs[i]) for i in IT_number_of_crystal_classes[k]])...)
                            for k=1:length(IT_number_of_crystal_classes)  ]

##

jul(x) = eval(Meta.parse(x))