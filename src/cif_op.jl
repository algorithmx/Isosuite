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


function extract_all_kw(
    cif, 
    kw::AbstractString
    )
    @assert (cif isa AbstractString) || (cif isa AbstractVector)
    cif_lines = String[]
    if cif isa AbstractString
        if !isfile(cif)
            return String[]
        else
            cif_lines = readlines(cif)
        end
    else
        cif_lines = cif[1:end]
    end

    p = findall(x->occursin(kw,x), cif_lines)
    return p===nothing ? String[] : cif_lines[sort(p)]
end

extract_all_kw(cif,kw::AbstractVector) = [extract_all_kw(cif,k) for k in kw]


function extract_kw(
    cif, 
    kw::AbstractString
    )
    @assert (cif isa AbstractString) || (cif isa AbstractVector)
    cif_lines = String[]
    if cif isa AbstractString
        if !isfile(cif)
            return String[]
        else
            cif_lines = readlines(cif)
        end
    else
        cif_lines = cif[1:end]
    end

    p = findfirst(x->occursin(kw,x), cif_lines)
    return p===nothing ? String[] : cif_lines[p]
end

extract_kw(cif,kw::AbstractVector) = [extract_kw(cif,k) for k in kw]


function get_number(cif, kw::AbstractString; default=0.0, parser=(x->parse(Float64,x)))
    l = extract_kw(cif, kw)
    if length(l) == 0
        return default
    else
        return parser(last(SPLTS(l)))
    end
end


@inline function get_number(cif, kw::AbstractVector; default=0.0, parser=(x->parse(Float64,x)))
    return [get_number(cif,k,default=default,parser=parser) for k in kw]
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

get_cell_params(cif) = get_number(  cif, 
                                    ["cell_length_a",
                                        "cell_length_b",
                                        "cell_length_c",
                                        "cell_angle_alpha",
                                        "cell_angle_beta",
                                        "cell_angle_gamma",
                                    ],
                                    parser=(x->parse(Float64,x))  )
                                                

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
    params = get_cell_params(cif)
    abc = [params[1],params[2],params[3]]
    angles = [params[4],params[5],params[6]]
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


function swap_xyz(l, perm)::String
    cmpnt = SPLTS(l)
    xyz = SPLTA(last(cmpnt))
    perm_inv = Dict(perm[i]=>i for i=1:3)
    rules = [["x","y","z"][perm[i]]=>["L","M","N"][i] for i=1:3]
    replall(x) = replace(replace(replace(x,rules[1]),rules[2]),rules[3])
    xyz_new = String[strip(replall(xyz[perm[i]])) for i=1:3]
    rules1 = [["L","M","N"][i]=>["x","y","z"][i] for i=1:3]
    replall1(x) = replace(replace(replace(x,rules1[1]),rules1[2]),rules1[3])
    return first(cmpnt) * "   " * replall1(join(xyz_new,","))
end


function swap_abc_by_perm(
    cif,
    perm_abc;
    id_xyz = 5,
    tol = 1e-6
    )

    d = Dict((2,3)=>1,(3,2)=>1,(1,2)=>3,(2,1)=>3,(1,3)=>2,(3,1)=>2)
    perm_angles = Int64[ d[(perm_abc[2],perm_abc[3])], d[(perm_abc[3],perm_abc[1])], d[(perm_abc[1],perm_abc[2])] ]
    cif_lines0 = (cif isa AbstractString) ? readlines(cif) : cif

    cif_lines = cif_lines0[1:end]

    @inline il(kw) = findfirst(x->occursin(kw,x), cif_lines)
    abc_line_ids = Int64[il("cell_length_a"), il("cell_length_b"), il("cell_length_c")]
    αβγ_line_ids = Int64[il("cell_angle_alpha"), il("cell_angle_beta"), il("cell_angle_gamma")]
    abc0 = cif_lines[abc_line_ids]
    αβγ0 = cif_lines[αβγ_line_ids]
    cif_lines[abc_line_ids[1]] = (@sprintf "_cell_length_a      %s" last(SPLTS(abc0[perm_abc[1]])))
    cif_lines[abc_line_ids[2]] = (@sprintf "_cell_length_b      %s" last(SPLTS(abc0[perm_abc[2]])))
    cif_lines[abc_line_ids[3]] = (@sprintf "_cell_length_c      %s" last(SPLTS(abc0[perm_abc[3]])))
    cif_lines[αβγ_line_ids[1]] = (@sprintf "_cell_angle_alpha   %s" last(SPLTS(αβγ0[perm_angles[1]])))
    cif_lines[αβγ_line_ids[2]] = (@sprintf "_cell_angle_beta    %s" last(SPLTS(αβγ0[perm_angles[2]])))
    cif_lines[αβγ_line_ids[3]] = (@sprintf "_cell_angle_gamma   %s" last(SPLTS(αβγ0[perm_angles[3]])))

    p_symm_op = findlast(x->occursin("space_group_symop",x), cif_lines)
    p = p_symm_op+1
    while !occursin("loop", cif_lines[p])
        cif_lines[p] = swap_xyz(cif_lines[p], perm_abc)
        p += 1
    end

    p_remove = findlast(x->occursin("atom_site_fract_symmform",x),cif_lines)
    cif_lines[p_remove] = ""
    @debug "swap_abc() : \natom_site_fract_symmform has been removed."

    pos = atom_config_pos(cif)
    atom_lines = cif_lines[pos:end]
    dp = pos-p_remove-1
    @debug "swap_abc() :  dp = $dp"

    ids(l) = [i for i in [collect(1:id_xyz-1); collect(id_xyz:id_xyz+2)[perm_abc]; collect(id_xyz+3:length(l))] if i!=length(l)-dp]
    swap_components(lx) = lx[ids(lx)]
    swap_a_line(l) = join( swap_components(SPLTS(l)), "  " )
    if pos <= 0
        @warn "swap_abc() has got cif with wrong format. Did nothing."
        return cif_lines0
    else
        return String[cif_lines[1:pos-1]; swap_a_line.(cif_lines[pos:end])] |> STRPRM
    end
end


swap_abc(
    cif;
    id_xyz = 5,
    tol = 1e-6
    ) =  swap_abc_by_perm(cif, abc_sortperm(cif), id_xyz=id_xyz, tol=tol)

#

function sort_atom_position_lines(
    cif
    )
    # section _atom_site_ from cif
    _atom_site_ = extract_all_kw(cif, "_atom_site_")
    id_xyz = findfirst(x->occursin("_atom_site_fract_x",x), _atom_site_)
    @assert findfirst(x->occursin("_atom_site_fract_z",x), _atom_site_) == id_xyz+2
    id_type  = findfirst(x->occursin("_atom_site_type_symbol",x), _atom_site_)
    id_label = findfirst(x->occursin("_atom_site_label",x), _atom_site_)

    # local functions
    @inline pf(s) = (abs(parse(Float64,s)<1e-8) ? 0.0 : parse(Float64,s)) 
    num2str(x) = (@sprintf "%10.6%" x)
    @inline correct_sign(x) = String[x[1:id_xyz-1]; num2str.(pf.(x[id_xyz:id_xyz+2])) ; x[id_xyz+3:end]]
    sortbyxyz(V) = V[sortperm(V,by=x->x[id_xyz:id_xyz+2])]
    @inline correct_label(x,lb) = String[x[1:id_label-1]; [lb,]; x[id_label+1:end]]
    @inline joinS(x) = join(x, "   ")

    cif_lines = (cif isa AbstractString) ? readlines(cif) : cif[1:end]
    pos_lines = SPLTS.(cif_lines)
    atm_unique = unique(map(x->x[id_type], pos_lines))
    pos_by_atm = [sortbyxyz([correct_sign(p) for p in pos if last(p)==a]) for a in atm_unique]
    pos_line_final = vcat([[correct_label(atm_group[i],"$(atm_group[i][id_type])$i") 
                            for i=1:length(atm_group)] 
                                for atm_group in pos_by_atm]...)  .|> joinS

    @debug "sort_atom_position_lines : "
    println.(lines)
    println("---")
    println.(pos_line_final)
    return pos_line_final
end
