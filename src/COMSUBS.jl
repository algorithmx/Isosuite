
function comsubs_block(
    cif_fn
    )
    IT = get_symmetry_Int_Tables_number(cif_fn)
    (a, b, c, alpha, beta, gamma) = get_cell_params(cif_fn)
    atom_pos = extract_atom_config(cif_fn)
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
    shuffle = 1.5,
    neighbor = 0.0,
    subgroup = (1,230)
    )
    return [
        [title,];
        comsubs_block(  cif1_fn  );
        comsubs_block(  cif2_fn  );
        [
            (@sprintf "size         %d" size),
            (@sprintf "strain     %8.4f  %8.4f" strain[1] strain[2]),
            (@sprintf "shuffle    %8.4f" shuffle),
            (@sprintf "neighbor   %8.4f" neighbor),
            (@sprintf "subgroup     %d  %d" subgroup[1] subgroup[2]),
        ];
    ]
end


function comsubs_output_section(res)
    p = findall(x->occursin("------------",x),res)
    if length(p)==0
        @info "Program comsubs didn't finished."
        return [res,]
    elseif length(p)==1
        @info "Program comsubs didn't find any common subgroups."
    end
    ps = [1, (p.+1)...]
    pe = [(p.-1)..., length(res)]
    return Vector{String}[res[i:j] for (i,j) ∈ zip(ps,pe)]
end


function comsubs_output_subgroup(sect)
    p1 = findfirst(x->occursin("Setting of crystal 1:",x), sect)
    p2 = findfirst(x->occursin("Setting of crystal 2:",x), sect)
    pm = findfirst(x->occursin("At midpoint:",x), sect)
    
    dic = Dict()
    dic["Subgroup"] = replace(sect[1], r"Subgroup\s+"=>"")

    for i=2:p1-1
        (k,v) = split(sect[i],":",keepempty=false)
        dic[strip(k)] = strip(v)
    end

    cryst1 = Dict()
    (k,v) = split(sect[p1+1],"=",keepempty=false)
    cryst1[strip(k)] = strip(v)
    for i=p1+2:p1+4
        (k,v) = split(sect[i],":",keepempty=false)
        cryst1[strip(k)] = strip(v)
    end
    cryst1["Wyckoff"] = sect[p1+5:p2-1]
    dic["Crystal 1"] = cryst1

    cryst2 = Dict()
    (k,v) = split(sect[p2+1],"=",keepempty=false)
    cryst2[strip(k)] = strip(v)
    for i=p2+2:p2+4
        (k,v) = split(sect[i],":",keepempty=false)
        cryst2[strip(k)] = strip(v)
    end
    cryst2["Wyckoff"] = sect[p2+5:pm-1]
    dic["Crystal 2"] = cryst2

    crystm = Dict()
    (k,v) = split(sect[pm+1],":",keepempty=false)
    crystm[strip(k)] = strip(v)
    crystm["Wyckoff"] = sect[pm+2:end]
    dic["Crystal m"] = crystm

    return dic
end
