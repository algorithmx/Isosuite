#ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/" ##!! modify this line first!!!

## --------------------------------------------

module Isosuite

    global const _Isosuite_package_root_ = pwd()

    using LinearAlgebra
    using Printf

    if "ISOSUITE_FOLDER" ∉ keys(ENV)
        @warn "\"ISOSUITE_FOLDER\" ∉ keys(ENV) \nProgram folder not found in ENV variables.\nPlease assign the Isotropy Suite path to ENV[\"ISOSUITE_FOLDER\"].\nDon't forget the slash \"/\" at the end of the path."
    end

    trslsh(l) = rstrip(l,['/'])
    adslsh(l) = trslsh(l)*"/"

    export SPLTX, SPLTN, SPLTS, SPLTS1, SPLTC, SPLTEQ, SPLTA, SPLTD, STRPRM
    export trim_comments_pound, decif, write_to_file, ⇶, ⦿
    
    export minimal_cif

    export get_title_line, get_atom_frac_pos
    export extract_kw, extract_all_kw
    export get_symmetry_Int_Tables_number, get_String
    export get_cell_length_a, get_cell_length_b, get_cell_length_c
    export get_cell_angle_alpha, get_cell_angle_beta, get_cell_angle_gamma
    export get_cell_params
    export extract_atom_config, compute_chemical_formula_structural
    export swap_abc, swap_abc_by_perm
    export sort_atom_position_lines
    export download_cif, download_cif_conventional, download_cif_primitive
    export supercell, symmetry_operators, loop_sections

    export iso, irrep_matrix, irrep_names, all_kvectors, all_elements

    export findsym, findsym_from_cif, findsym_cifinput, findsym_input, extract_space_group
    export extract_cif, extract_wyckoff, extract_atoms, atom_to_wyckoff
    export generate_cif, improve_cif, improve_cif__findsym_cifinput

    export smodes, input_smodes, translate_smodes_result

    export comsubs_input

    export QE_default_equivalent_settings_findsym, default_settings_findsym
    export QE_default_equivalent_settings_iso, default_settings_iso

    export WYCKPOS_TABLE

    include("SG_settings.jl")

    include("parsers.jl")

    include("cif_op.jl")
    
    include("string_op.jl")

    include("wyckpos.jl")
    
    include("submit.jl")

    include("ISOTROPY.jl")

    include("SMODES.jl")

    include("FINDSYM.jl")

    include("COMSUBS.jl")

    # include("........jl")

end # module
