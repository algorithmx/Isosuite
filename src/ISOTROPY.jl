# to download ISOTROPY suite, visit
# https://stokes.byu.edu/iso/isotropy.php

import LinearAlgebra.transpose
transpose(s::String) = s

##

@inline starts_at(key, str) = findfirst(key,str)[1]

function split_content_line( 
    ps::Vector{Int}, 
    line::String 
    )::Vector{String}
    S_E = zip(ps,[ps[2:end];[length(line)+1]])
    return String[strip(line[s:e-1]) for (s,e) ∈ S_E]
end

function parse_results(
    result_lines::Vector{S}, 
    keywords::Vector{String}
    ) where {S<:AbstractString}
    @info "raw lines : "
    println.(result_lines)
    if length(result_lines) == 0
        @warn "empty result."
        return nothing
    end
    @inline rm_continue(lines) = String[l for l in lines if !occursin("RETURN",l) && all(!occursin(kw,l) for kw in keywords)]
    title     = result_lines[1]
    contents  = rm_continue(result_lines[2:end])
    @info "result lines : "
    @info "title : "
    println(title)
    @info "contents : "
    println.(contents)
    println("[END]")
    positions = [starts_at(k,title) for k ∈ keywords]
    return hcat( [split_content_line(positions,l) 
                 for l ∈ contents if length(strip(l)) > 0]... ) |> transpose |> Matrix
end


function parse_results1(
    result_lines::Vector{S}, 
    keywords::Vector{String}
    ) where {S<:AbstractString}
    if length(result_lines) == 0
        @warn "empty result."
        return nothing
    end
    @inline rm_continue(lines) = String[l for l in lines if !occursin("RETURN",l) && all(!occursin(kw,l) for kw in keywords)]
    title     = result_lines[1]
    contents  = rm_continue(result_lines[2:end])
    return hcat( [SPLTS(l) for l ∈ contents if length(strip(l)) > 0]... ) |> transpose |> Matrix
end

##

function parse_number(ns)
    n = try
        parse(Int64,ns)
    catch
        parse(Float64,ns)
    end
    return n
end


function parse_matrix(res)::Matrix
    hcat( [parse_number.(split(r, ' ', keepempty=false)) 
          for r ∈ res]...) |> transpose
end

## ------------------------- recipies -----------------------------

function all_elements(
    parent::Int64;
    setting="MILLER-LOVE",
    )

    comms = ["VALUE PARENT $parent",
             "SETTING $setting",
             "SHOW PARENT",
             "SHOW ELEMENTS",
             "DISPLAY PARENT",
             "QUIT"]

    s = iso(comms)

    p = parse_results(s, ["Parent", "Elements"])

    el = (p===nothing) ? String[] : split(join(p[:,2]," "),", ",keepempty=false)

    return Dict("Parent"   => p[1,1],
                "Elements" => String.(el))
end



function all_kvectors(
    parent::Int64;
    setting="MILLER-LOVE",
    )
    comms = ["VALUE PARENT $parent",
             "SETTING $setting",
             "SHOW KPOINT",
             "SHOW KDEGREE",
             "DISPLAY KPOINT",
             "\n",
             "\n",
             "\n",
             "QUIT"]

    s = iso(comms)

    p = parse_results1(s, ["k vector", "k degree"])

    return Dict(p[k,1]=>(p[k,2],parse(Int,p[k,3])) for k = 1:size(p,1))
end


function irrep_names(
    parent::Int64, 
    kpoint::String;
    setting="MILLER-LOVE",
    )

    comms = ["VALUE PARENT $parent",
             "SETTING $setting",
             "VALUE KPOINT $kpoint",
             "SHOW IRREP",
             "DISPLAY IRREP",
             "QUIT"]
    
    s = iso(comms)

    p = parse_results(s, ["Irrep"])

    return Dict("Parent"  => parent,
                "Kpoint"  => kpoint, 
                "Irreps"  => (p===nothing) ? String[] : p)
end


function irrep_matrix(
    parent::Int64, 
    ir::String,
    elem::String;
    setting="MILLER-LOVE",
    )
    
    comms = ["VALUE PARENT $parent",
             "SETTING $setting",
             "VALUE IRREP $ir",
             "VALUE ELEMENT $elem",
             "SHOW IRREP",
             "SHOW CHARACTER",
             "SHOW MATRIX",
             "DISPLAY IRREP",
             "QUIT"]
    
    s = iso(comms)
    
    p = parse_results(s, ["Irrep", "Element", "Char", "Matrix"])

    return Dict("Parent"    => parent,
                "Irrep"     => (p===nothing) ? ir      : p[1,1], 
                "Element"   => (p===nothing) ? elem    : p[1,2], 
                "Character" => (p===nothing) ? nothing : parse_number(p[1,3]),
                "Matrix"    => (p===nothing) ? nothing : parse_matrix(p[:,4]))
end


## ----------------------------------------------

##

# W1 = irrep_matrix(225,"W1","E 0 0 0")

##

# irrep_names(221, "X")

##

# kp221 = all_kvectors(221)