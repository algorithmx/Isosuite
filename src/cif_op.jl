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

function extract_kw(
    cif_fn::S, 
    kw
    ) where {S<:AbstractString}
    return readlines(pipeline(`cat $cif_fn`, `grep $kw`))
end

function extract_kw(
    cif_lines::Vector{S}, 
    kw
    ) where {S<:AbstractString}
    i=rand(0:99999)
    fn = "/tmp/extract_kw.$i.tmp"
    cif_lines >>> fn
    res = extract_kw(fn, kw)
    rm(fn)
    return res
end

function get_number(cif, kw::AbstractString; default=0.0, parser=(x->parse(Float64,x)))
    l = extract_kw(cif, kw)
    if l === nothing || length(l) == 0
        return default
    else
        return parser(last(SPLTS(first(l))))
    end
end

get_Int(cif, kw) = get_number(cif, kw; default=1, parser=(x->parse(Int,x)))

get_Float(cif, kw) = get_number(cif, kw; default=0.0, parser=(x->parse(Float64,x)))

get_String(cif, kw) = get_number(cif, kw; default=0.0, parser=(x->string(x)))

get_symmetry_Int_Tables_number(cif) = get_Int(cif, "symmetry_Int_Tables_number")

get_cell_length_a(cif) = get_Float(cif, "cell_length_a")
get_cell_length_b(cif) = get_Float(cif, "cell_length_b")
get_cell_length_c(cif) = get_Float(cif, "cell_length_c")

get_cell_angle_alpha(cif) = get_Float(cif, "cell_angle_alpha")
get_cell_angle_beta(cif)  = get_Float(cif, "cell_angle_beta")
get_cell_angle_gamma(cif) = get_Float(cif, "cell_angle_gamma")


function atom_config_pos(
    cif
    )::Int
    cif_lines = (cif isa AbstractString) ? readlines(cif) : cif[1:end]
    atm_lines = cif_lines[findall(x->occursin("_atom_site_",x), cif_lines)]
    pos = findlast(x->occursin(last(atm_lines),x), cif_lines)
    if pos === nothing
        @warn "extract_config($cif_fn) has got a cif file with wrong format."
        return -1
    end
    return pos+1
end


function extract_atom_config(
    cif
    )::Vector{String}
    cif_lines = (cif isa AbstractString) ? readlines(cif) : cif[1:end]
    pos = atom_config_pos(cif)
    return pos>0 ? cif_lines[pos:end] : String[] 
end


function compute_chemical_formula_structural(
    cif; 
    atom_type_pos=2, 
    multiplicity_pos=3
    )::Dict
    config_lines = extract_atom_config(cif)
    @inline atm_type(x) = x[atom_type_pos]
    @inline atm_mult(x) = x[multiplicity_pos]
    atm = zip(atm_type.(SPLTS.(config_lines)), atm_mult.(SPLTS.(config_lines)))
    return Dict(a=>sum([parse(Int,last(x)) for x in atm if first(x)==a]) for a in unique(first.(atm)))
end

function abc_sortperm(
    cif;
    tol = 1e-6
    )::Vector{Int64}

    close(x,y) = abs(x-y) < tol
    abc = [get_cell_length_a(cif), get_cell_length_b(cif), get_cell_length_c(cif)]
    angles = [get_cell_angle_alpha(cif), get_cell_angle_beta(cif) , get_cell_angle_gamma(cif)]
    d = Dict((2,3)=>1,(3,2)=>1,(1,2)=>3,(2,1)=>3,(1,3)=>2,(3,1)=>2)
    o = sortperm(abc)
    if close(abc[o[1]], abc[o[3]]) && close(abc[o[2]], abc[o[3]]) && close(abc[o[1]], abc[o[2]])
        a = sortperm(angles)
        if close(angles[a[1]],angles[a[2]]) && close(angles[a[1]],angles[a[3]]) && close(angles[a[2]],angles[a[3]])
            return [1,2,3]
        elseif close(angles[a[1]],angles[a[2]])
            return a
        elseif close(angles[a[1]],angles[a[3]])
            return a[[1,3,2]]
        elseif close(angles[a[2]],angles[a[3]])
            return a[[2,3,1]]
        else
            return a
        end
    elseif close(abc[o[1]], abc[o[2]])
        if angles[d[(o[1],o[3])]] < angles[d[(o[2],o[3])]]
            return o
        else
            return o[[2,1,3]]
        end
    elseif close(abc[o[2]], abc[o[3]])
        if angles[d[(o[1],o[2])]] < angles[d[(o[1],o[3])]]
            return o[[2,3,1]]
        else
            return o[[3,2,1]]
        end
    else
        return o
    end
end


function swap_abc(
    cif;
    id_xyz = 5,
    tol = 1e-6
    )

    perm_abc = abc_sortperm(cif)
    d = Dict((2,3)=>1,(3,2)=>1,(1,2)=>3,(2,1)=>3,(1,3)=>2,(3,1)=>2)
    perm_angles = Int64[ d[(perm_abc[2],perm_abc[3])], d[(perm_abc[3],perm_abc[1])], d[(perm_abc[1],perm_abc[2])] ]
    cif_lines0 = (cif isa AbstractString) ? readlines(cif) : cif
    
    cif_lines = cif_lines0[1:end]

    @inline il(kw) = findfirst(x->occursin(kw,x), cif_lines)
    abc_line_ids = Int64[il("cell_length_a"), il("cell_length_b"), il("cell_length_c")]
    αβγ_line_ids = Int64[il("cell_angle_alpha"), il("cell_angle_beta"), il("cell_angle_gamma")]
    abc0 = cif_lines[abc_line_ids]
    αβγ0 = cif_lines[αβγ_line_ids]
    cif_lines[abc_line_ids[1]] = abc0[perm_abc[1]]
    cif_lines[abc_line_ids[2]] = abc0[perm_abc[2]]
    cif_lines[abc_line_ids[3]] = abc0[perm_abc[3]]
    cif_lines[αβγ_line_ids[1]] = αβγ0[perm_angles[1]]
    cif_lines[αβγ_line_ids[2]] = αβγ0[perm_angles[2]]
    cif_lines[αβγ_line_ids[3]] = αβγ0[perm_angles[3]]

    pos = atom_config_pos(cif)
    atom_lines = cif_lines[pos:end]
    swap_components(lx) = lx[[collect(1:id_xyz-1); collect(id_xyz:id_xyz+2)[perm_abc]; collect(id_xyz+3:length(lx))]]
    swap_a_line(l) = join( swap_components(SPLTS(l)), "  " )
    if pos <= 0
        @warn "swap_abc() has got cif with wrong format. Did nothing."
        return cif_lines0
    else
        return String[cif_lines[1:pos-1]; swap_a_line.(cif_lines[pos:end])]
    end
end
