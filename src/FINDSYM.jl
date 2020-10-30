using Printf

a_number(t) = (@sprintf "%10.8f" t)

six_numbers(t) = (@sprintf "%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f" t[1] t[2] t[3] t[4] t[5] t[6])

three_numbers_from_4tuple(t) = (@sprintf "%12.8f  %12.8f  %12.8f" t[2] t[3] t[4])

count_occurance(a,l) = sum(first.(l).==a)

n_times_atom(a,l) = "$(count_occurance(a,l))*$a"

unique_atoms(l) = sort(unique(first.(l)))

atom_list(l) = join([n_times_atom(a,l) for a in unique_atoms(l)], "  ")

atom_pos(l) = vcat([three_numbers_from_4tuple.(l[findall(first.(l).==a)]) for a in unique_atoms(l)]...)

#atm = [("Sc",0.0,0.0,0.0),("F",0.5,0.0,0.0),("F",0.0,0.5,0.0),("F",0.0,0.0,0.5)]

#atom_list(atm)
#atom_pos(atm)


function findsym_input(
    TITLE,
    latticeParameters::Tuple,
    atomList::Vector;
    latticeTolerance = 0.00001,
    atomicPositionTolerance = 0.00001,
    )
    @assert all(length.(atomList).==4)
    @assert all(t<:AbstractString for t ∈ typeof.(first.(atomList)))
    l0 = [
        "!useKeyWords",
        "!title",
        "$TITLE",
        "!latticeTolerance",
        "$(a_number(latticeTolerance))",
        "!atomicPositionTolerance",
        "$(a_number(atomicPositionTolerance))",
        "!latticeParameters",
        "$(six_numbers(latticeParameters))",
        "!atomCount",
        "$(length(atomList))",
        "!atomType",
        "$(atom_list(atomList))",
        "!atomPosition",
    ]
    return join([l0; atom_pos(atomList)], "\n")
end


##!! TODO  modify tolerance lines
# submit scripts to the program `findsym_cifinput`
function findsym_from_cif(
    fn_cif::AbstractString;
    latt_tol = 1e-5,
    pos_tol = 1e-5) 
    res0 = findsym_cifinput(fn_cif)
    ##!! TODO  modify tolerance lines below
    res1 = res0 |> STRPRM |> trim_comments_pound
    findsym(res1)
end


function generate_cif(
    title::AbstractString, 
    scripts::Vector{S}
    ) where {S<:AbstractString}
    cif_lines = scripts |> findsym |> extract_cif
    @assert occursin("data_findsym-output", cif_lines[1])
    cif_lines[1] = title
    return STRPRM(cif_lines)
end

generate_cif(title::AbstractString, script::AbstractString)  = generate_cif(title, SPLTN(script))

improve_cif(title::AbstractString, old_cif_fn::AbstractString) = generate_cif(title, findsym_cifinput(old_cif_fn))

improve_cif(title::AbstractString, old_cif_fn::AbstractString, new_cif_fn::AbstractString) = generate_cif(title, findsym_cifinput(old_cif_fn)) >>> new_cif_fn


function extract_cif(outp)
    cif_line = findfirst(x->occursin(r"\#\s+CIF\s+file",x), outp)
    cif_line_end = findfirst(x->occursin(r"\#\s+end\s+of\s+cif",x), outp)
    @assert cif_line_end !== nothing
    return outp[(cif_line+1):(cif_line_end-1)] |> STRPRM |> trim_comments_pound
end


function extract_atoms(outp)
    atm_line = findfirst(x->occursin(r"Type\s+of\s+each\s+atom",x), outp)
    atm_line_end = findfirst(x->occursin("Tolerance",x), outp)
    @assert atm_line_end !== nothing
    ATOM_LINES = STRPRM(outp[(atm_line+1):(atm_line_end-1)])

    @assert occursin(r"((\d+\s*\*\s*)?\w+\s*)+",ATOM_LINES[1])
    @assert occursin(r"[Pp]osition",ATOM_LINES[2]) && occursin(r"[Oo]cc",ATOM_LINES[2])
    lines = [SPLTS(l) for l in ATOM_LINES[3:end]]
    # example output 
    # Dict( 4 => ("O", [0.75, 0.25, 0.00291], 1.0),
    #       2 => ("O", [0.0, 0.0, 0.73138], 1.0),
    #       3 => ("O", [0.25, 0.25, 0.00291], 1.0),
    #       1 => ("U", [0.0, 0.0, 0.27215], 1.0))
    return Dict(parse(Int,l[1])=>(string(l[2]),parse.(Float64,l[3:end-1]), parse(Float64,l[end])) for l in lines)
end


function extract_wyckoff(outp)
    wyc_line = findfirst(x->occursin("Atomic positions and occupancies in terms",x), outp)
    wyc_line_end = findfirst(x->occursin(r"\#\s*CIF",x), outp)
    @assert wyc_line_end !== nothing
    WYCKOFF_LINES = outp[(wyc_line+1):(wyc_line_end-2)]

    WYCKS = []
    @inline get_wyckoff_symbol(l) = string(SPLTS(l)[3])
    @inline get_atom_label(l) = string(strip(SPLTS(l)[4],['(',')']))
    @inline parse_xyz(ls) = [
        (t1=findfirst(m->occursin("x",m),ls); (t1===nothing ? 0.0 : parse(Float64,SPLTEQ(ls[t1])[2]))),
        (t2=findfirst(m->occursin("y",m),ls); (t2===nothing ? 0.0 : parse(Float64,SPLTEQ(ls[t2])[2]))),
        (t3=findfirst(m->occursin("z",m),ls); (t3===nothing ? 0.0 : parse(Float64,SPLTEQ(ls[t3])[2]))),
    ]
    @inline get_xyz(ls) = length(ls)==1 ? [0.0,0.0,0.0] : parse_xyz(ls[2:end])
    @inline get_atom_id(l) = parse(Int, first(SPLTS(l)))
    BS = findall(x->occursin(r"Wyckoff\s+",x), WYCKOFF_LINES)
    BN = [BS[2:end].-1; length(WYCKOFF_LINES)]
    for (s,f) ∈ zip(BS,BN)
        @info WYCKOFF_LINES[s]
        wls = SPLTA(WYCKOFF_LINES[s])
        symb = get_wyckoff_symbol(wls[1])
        atml = get_atom_label(wls[1])
        xyz = get_xyz(wls)
        atm_ids = [get_atom_id(l) for l in WYCKOFF_LINES[s+1:f] if startswith(strip(l), r"\d+")]
        push!(WYCKS, symb=>(atm_ids,xyz,atml))
    end
    # example output
    # ["a" => ([1], [0.0, 0.0, 0.0], "U1"), "d" => ([2, 3], [0.0, 0.0, 0.0], "O1"), "b" => ([4], [0.0, 0.0, 0.0], "O2")]
    return WYCKS
end


function atom_to_wyckoff(outp)
    atoms = extract_atoms(outp)
    unique_atoms = unique(first.(values(atoms)))
    wycks = extract_wyckoff(outp)
    @info string(wycks)
    @inline atm_id(a) = [k for (k,v) ∈ atoms if v[1]==a]
    @inline id_wyck(i) = first([((w,ids_xyz_atml[2:end]...),i) for (w,ids_xyz_atml) ∈ wycks if (i ∈ ids_xyz_atml[1])])
    return Dict(a=>id_wyck.(atm_id(a)) for a in unique_atoms)
end


function extract_space_group(outp)
    for l in outp
        if occursin(r"Space\s+Group\s+\d+", l)
            return parse(Int, split(l," ",keepempty=false)[3])
        end
    end
    return 0
end

