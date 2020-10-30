ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/" ##!! modify this line first!!!

## --------------------------------------------

module Isosuite

    using Printf

    if "ISOSUITE_FOLDER" ∉ keys(ENV)
        @warn "\"ISOSUITE_FOLDER\" ∉ keys(ENV) \nProgram folder not found in ENV variables.\nPlease define ENV[\"ISOSUITE_FOLDER\"]."
    end

    trslsh(l) = rstrip(l,['/'])

    export SPLTN, SPLTS, SPLTS1, SPLTC, SPLTEQ, SPLTA, STRPRM, write_to_file, >>>, ⦿
    
    export minimal_cif, extract_lattice_parameters

    export iso, irrep_matrix, irrep_names, all_kvectors, all_elements

    export findsym, findsym_from_cif, findsym_cifinput, findsym_input, extract_space_group
    export extract_cif, extract_wyckoff, extract_atoms, atom_to_wyckoff
    export generate_cif, improve_cif

    export smodes, input_smodes, translate_smodes_result

    include("cif_op.jl")
    
    include("string_op.jl")
    
    include("submit.jl")

    include("ISOTROPY.jl")

    include("SMODES.jl")

    include("FINDSYM.jl")

    # include("........jl")

end # module
