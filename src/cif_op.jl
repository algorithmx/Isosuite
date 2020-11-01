function minimal_cif(
    title,
    latt_params::Tuple,
    atom_frac_positions::Vector
    )

    __cif__ = """TITLE

    _chemical_name_common                  'TITLE'
    _cell_length_a                         AAA
    _cell_length_b                         BBB
    _cell_length_c                         CCC
    _cell_angle_alpha                      aaa
    _cell_angle_beta                       bbb
    _cell_angle_gamma                      ggg
    _space_group_name_H-M_alt              'P 1'
    _space_group_IT_number                 1
    
    loop_
    _space_group_symop_operation_xyz
       'x, y, z'
    
    loop_
    _atom_site_label
    _atom_site_occupancy
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_adp_type
    _atom_site_B_iso_or_equiv
    _atom_site_type_symbol
    """
    #  U1     1.0     0.000000      0.000000      0.000000     Biso  1.000000 U


    (AAA,BBB,CCC,alpha,beta,gamma) = latt_params

    str4(x) = (@sprintf "%8.4f" x)
    str5(x) = (@sprintf "%9.5f" x)
    rules = [   "TITLE" => title,
                "AAA"=>str5(AAA), "BBB"=>str5(BBB), "CCC"=>str5(CCC), 
                "aaa"=>str5(alpha), "bbb"=>str5(beta), "ggg"=>str5(gamma)  ]

    cif_str = __cif__
    for p in rules
        cif_str = replace(cif_str, p)
    end

    @inline remove_num(s) = strip(replace(strip(s), r"\d+"=>""))
    frac_lines = String[]
    for (i,a) in enumerate(atom_frac_positions)
        line = (@sprintf  "%s  1.0   %9.5f   %9.5f  %9.5f  Biso  1.0   %s" a[1]*string(i) a[2]  a[3]  a[4] remove_num(a[1]))
        push!(frac_lines, line)
    end

    return cif_str * ⦿(frac_lines)
end


function extract_lattice_parameters(
    cif_lines::Vector{S}
    )::Tuple where {S<:AbstractString}
    kw = [  "_cell_length_a"
            "_cell_length_b"
            "_cell_length_c"
            "_cell_angle_alpha"
            "_cell_angle_beta"
            "_cell_angle_gamma" ]
    @inline get_line(lines,w) = cif_lines[findfirst(x->occursin(w,x), lines)]
    @inline get_number(l) = parse(Float64, SPLTS(strip(l))[2])
    return Tuple( collect((get_number(get_line(cif_lines,w)) for w in kw)) )
end


extract_lattice_parameters(cif_fn::S) where {S<:AbstractString} = extract_lattice_parameters(readlines(cif_fn))

extract_kw(cif_fn::S, kw) where {S<:AbstractString} = readlines(pipeline(`cat $cif_fn`, `grep $kw`))

function get_number(cif_fn, kw; default=0.0, parser=(x->parse(Float64,x)))
    l = extract_kw(cif_fn, kw)
    if l === nothing || length(l) == 0
        return default
    else
        return parser(last(SPLTS(first(l))))
    end
end

get_Int(cif_fn, kw) = get_number(cif_fn, kw; default=1, parser=(x->parse(Int,x)))

get_Float(cif_fn, kw) = get_number(cif_fn, kw; default=0.0, parser=(x->parse(Float64,x)))

get_String(cif_fn, kw) = get_number(cif_fn, kw; default=0.0, parser=(x->string(x)))

get_symmetry_Int_Tables_number(cif_fn) = get_Int(cif_fn, "symmetry_Int_Tables_number")

get_cell_length_a(cif_fn) = get_Float(cif_fn, "cell_length_a")
get_cell_length_b(cif_fn) = get_Float(cif_fn, "cell_length_b")
get_cell_length_c(cif_fn) = get_Float(cif_fn, "cell_length_c")

get_cell_angle_alpha(cif_fn) = get_Float(cif_fn, "cell_angle_alpha")
get_cell_angle_beta(cif_fn)  = get_Float(cif_fn, "cell_angle_beta")
get_cell_angle_gamma(cif_fn) = get_Float(cif_fn, "cell_angle_gamma")

