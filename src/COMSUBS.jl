
function comsubs_block(cif_fn)
    cif_new = improve_cif("", cif_fn, SG_setting=default_settings_findsym)
    IT = get_symmetry_Int_Tables_number(cif_new)
    (a, b, c, alpha, beta, gamma) = get_cell_params(cif_new)
    atom_pos = extract_atom_config(cif_new)
    # the improved cif from findsym has fixed order :
    # _atom_site_label
    # _atom_site_type_symbol
    # _atom_site_symmetry_multiplicity
    # _atom_site_Wyckoff_label
    # _atom_site_fract_x
    # _atom_site_fract_y
    # _atom_site_fract_z
    # _atom_site_occupancy
    # _atom_site_fract_symmform
    wk(s) = match(r"\s+[a-n]\s+(-)?\d+\.\d+\s+(-)?\d+\.\d+\s+(-)?\d+\.\d+",s).match
    am(s) = SPLTS(s)[2]
    wyck_lines = ["$a  $w" for (a,w) in zip(am.(atom_pos), wk.(atom_pos))]
    return [
        [ "$IT   ! space group symmetry ",
          (strip((@sprintf "%12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f"  a  b  c  alpha  beta  gamma)) 
          * "   ! lattice parameters: a,b,c,alpha,beta,gamma"),
          "$(length(wyck_lines))   ! number of Wyckoff positions"
        ];
        wyck_lines
    ]
end


function comsubs_input(
    cif1_fn::AbstractString, 
    cif2_fn::AbstractString;
    title = "generated by comsub_input",
    size = 4,
    strain = (0.95,1.05),
    shuffle = 0.5,
    neighbor = 0.0,
    subgroup = (1,230)
    )
    return [
        [title,];
        comsubs_block(cif1_fn);
        comsubs_block(cif2_fn);
        [
            (@sprintf "size         %d" size),
            (@sprintf "strain     %8.4f  %8.4f" strain[1] strain[2]),
            (@sprintf "shuffle    %8.4f" shuffle),
            (@sprintf "neighbor   %8.4f" neighbor),
            (@sprintf "subgroup     %d  %d" subgroup[1] subgroup[2]),
        ];
    ]
end
