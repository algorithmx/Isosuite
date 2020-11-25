function get_title_line(cif)
    cif_lines = ( (cif isa AbstractString) ? (isfile(cif) ? readlines(cif) : String[]) : cif ) |> STRPRM |> trim_comments_pound
    if length(cif_lines) == 0
        @info "get_title_line(): \nEmpty content."
        return ""
    else
        return cif_lines[1]
    end
end


function loop_sections(cif)
    cif_lines = (cif isa Vector) ? cif : readlines(cif)
    lpos0 = findall(x->occursin("loop_",x),cif_lines)
    lpos1 = [[1]; lpos0]
    lpos2 = [lpos0; length(cif_lines)+1]
    return [cif_lines[a:b-1] |> STRPRM  for (a,b) in zip(lpos1,lpos2)]
end



function extract_all_kw(
    cif, 
    kw::Union{AbstractString,Regex}
    )
    @assert (cif isa AbstractString) || (cif isa AbstractVector)
    cif_lines = String[]
    if cif isa AbstractString
        if !isfile(cif)
            @warn "extract_all_kw() : \n$cif file not found."
            return String[]
        else
            cif_lines = readlines(cif)
        end
    else
        cif_lines = cif[1:end]
    end

    p = findall(x->occursin(kw,x), cif_lines)
    if p===nothing
        @info "extract_all_kw() : \n$kw not found."
        return String[]
    else
        return cif_lines[sort(p)]
    end
end

extract_all_kw(cif,kw::AbstractVector) = [extract_all_kw(cif,k) for k in kw]


function extract_kw(
    cif, 
    kw::Union{AbstractString,Regex}
    )::String
    @assert (cif isa AbstractString) || (cif isa AbstractVector)
    cif_lines = String[]
    if cif isa AbstractString
        if !isfile(cif)
            @warn "extract_kw() : \n$cif file not found."
            return ""
        else
            cif_lines = readlines(cif)
        end
    else
        cif_lines = cif[1:end]
    end
    
    p = findfirst(x->occursin(kw,x), cif_lines)
    if p===nothing
        @info "extract_kw() : \n$kw not found."
        return ""
    else
        return cif_lines[p]
    end
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
    if pos <=0
        @warn "extract_atom_config() has pos < 0."
    end
    return pos>0 ? cif_lines[pos:end] : String[] 
end


##! TODO deal with Fddd (UO3) a,c,-b
function apply_space_group__transform_Pp_abc() end


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


function swap_abc(
    cif;
    tol = 1e-6
    )
    _atom_site_ = extract_all_kw(cif, "_atom_site_")
    id_xyz = findfirst(x->occursin("_atom_site_fract_x",x), _atom_site_)
    return  swap_abc_by_perm(cif, abc_sortperm(cif), id_xyz=id_xyz, tol=tol)
end


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
    num2str(x) = (@sprintf "%10.6f" x)
    @inline correct_sign(x) = String[x[1:id_xyz-1]; num2str.(pf.(x[id_xyz:id_xyz+2])) ; x[id_xyz+3:end]]
    sortbyxyz(V) = V[sortperm(V,by=x->x[id_xyz:id_xyz+2])]
    @inline correct_label(x,lb) = String[x[1:id_label-1]; [lb,]; x[id_label+1:end]]
    @inline joinS(x) = join(x, "   ")

    pos_lines = SPLTS.(extract_atom_config(cif))
    atm_unique = unique(map(x->x[id_type], pos_lines))
    pos_by_atm = [  sortbyxyz([correct_sign(p) for p in pos_lines if p[id_type]==a])
                    for a in atm_unique ]

    pos_line_final = vcat([[correct_label(atm_group[i],"$(atm_group[i][id_type])$i") 
                            for i=1:length(atm_group)] 
                                for atm_group in pos_by_atm]...)  .|> joinS

    cif_lines = (cif isa AbstractString) ? readlines(cif) : cif[1:end]
    atm_p = atom_config_pos(cif)
    return [cif_lines[1:atm_p-1] ; pos_line_final]
end


function get_atom_frac_pos(
    cif
    )
    # section _atom_site_ from cif
    _atom_site_ = extract_all_kw(cif, "_atom_site_")
    id_xyz = findfirst(x->occursin("_atom_site_fract_x",x), _atom_site_)
    @assert findfirst(x->occursin("_atom_site_fract_z",x), _atom_site_) == id_xyz+2
    id_type  = findfirst(x->occursin("_atom_site_type_symbol",x), _atom_site_)

    # local functions
    @inline pf(s) = (abs(parse(Float64,s)<1e-10) ? 0.0 : parse(Float64,s)) 

    pos_lines = SPLTS.(extract_atom_config(cif))
    return map(x->[x[id_type], pf(x[id_xyz]), pf(x[id_xyz+1]), pf(x[id_xyz+2])], pos_lines) |> sort
end # |> atom_list


function symmetry_operators(cif)
    ops   = String[]
    sect  = loop_sections(cif)
    symm0 = [i for i=1:length(sect) if "_space_group_symop_operation_xyz" in sect[i]]
    if length(symm0)==0
        @info "symmetry_operators():\nNo symmetry operators found in cif."
        return ops
    end
    symm = sect[first(symm0)]
    ops  = [replace(
                replace(
                    replace(x, r"^\d+\s+"=>""), 
                    r"\s+"=>""), 
                "'"=>""
            ) for x in symm if !occursin("_",x)]
    return ops
end # |> symm_ops


function extend_positions(atom_list, symm_ops)
    if length(symm_ops)==0
        return atom_list
    elseif length(symm_ops)==1
        if symm_ops[1]=="x,y,z"
            return atom_list
        else
            @error "extend_positions():\nWrong format of symmetry operators.\nGot len=1 with symm_ops[1]=$(symm_ops[1]).\nNo operation performed."
            return atom_list
        end
    end
    close(lst, p) = any([norm([l[2:4]...].-[p...])<1e-5 for l in lst])
    mod1(l) = map(x->round(mod(x+2,1),digits=6), l)    
    atms = unique(first.(atom_list))
    ext_pos = []
    for a in atms
        ext_pos_a = []
        apos = findall(x->x[1]==a, atom_list)
        for p in apos
            line = atom_list[p]
            for sop in symm_ops
                new_p = mod1( eval(Meta.parse("ft(x,y,z)=[$sop]; ft($(line[2]),$(line[3]),$(line[4]))")) )
                if !close(ext_pos_a, new_p)
                    push!(ext_pos_a, (line[1], new_p...))
                end
            end
        end
        ext_pos = [ext_pos; ext_pos_a]
    end
    return sort(ext_pos)
end


## ---------------------------------------------------------

function download_cif(url::AbstractString, fn="")
    outp_fn = length(fn)>0 ? fn : last(SPLTX(url, "/"))
    try
        run(`wget $url -O $outp_fn`)
    catch
        @warn "download_cif() :\nwget command exited with error."
    end
    return
end


download_cif_conventional(mp_number::Int) = download_cif("https://materialsproject.org/materials/mp-$(mp_number)/cif?type=conventional_standard&download=true", "mp-$(mp_number).cif")

download_cif_primitive(mp_number::Int) = download_cif("https://materialsproject.org/materials/mp-$(mp_number)/cif?type=primitive&download=true", "mp-$(mp_number).cif")


## ---------------------------------------------------------


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
    _symmetry_Int_Tables_number            IIITTT
    
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


    (AAA,BBB,CCC,alpha,beta,gamma) = latt_params[1:6]
    IT = (length(latt_params)==7) ? latt_params[7] : 1

    str4(x) = (@sprintf "%8.4f" x)
    str5(x) = (@sprintf "%9.5f" x)
    rules = [   "TITLE" => title,
                "AAA"=>str5(AAA), "BBB"=>str5(BBB), "CCC"=>str5(CCC), 
                "aaa"=>str5(alpha), "bbb"=>str5(beta), "ggg"=>str5(gamma),
                "IIITTT"=>string(IT)  ]

    cif_str = __cif__
    for p in rules
        cif_str = replace(cif_str, p)
    end

    @inline remove_num(s) = strip(replace(strip(string(s)), r"\d+"=>""))
    frac_lines = String[]
    for (i,a) in enumerate(atom_frac_positions)
        line = (@sprintf  "%s  1.0   %9.5f   %9.5f  %9.5f  Biso  1.0   %s" a[1]*string(i) a[2]  a[3]  a[4] remove_num(a[1]))
        push!(frac_lines, line)
    end

    return cif_str * ⦿(frac_lines)
end


## ---------------------------------------------------------


function supercell(
    cif,
    nnn;
    SG_setting = default_settings_findsym
    )
    symm_ops = symmetry_operators(cif)
    atom_list = get_atom_frac_pos(cif)
    atom_list_ext = extend_positions(atom_list, symm_ops)
    atom_list_ext_enlarge = []
    (nx, ny, nz) = nnn
    for i=1:nx
        for j=1:ny
            for k=1:nz
                for l in atom_list_ext
                    push!(atom_list_ext_enlarge, (l[1], (l[2]+(i-1))/nx, (l[3]+(j-1))/ny, (l[4]+(k-1))/nz))
                end
            end
        end
    end
    (a, b, c, alpha, beta, gamma) = get_cell_params(cif)
    IT = get_symmetry_Int_Tables_number(cif)
    cif_enlarge = minimal_cif(
        get_title_line(cif) * " super cell ($nx, $ny, $nz)",
        (nx*a, ny*b, nz*c, alpha, beta, gamma, IT),
        atom_list_ext_enlarge
    )
    write_to_file(cif_enlarge, "tmp.cif")
    return cif_enlarge |> SPLTN
end


## ------------------------------------------------------------------
