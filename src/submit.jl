

## -----------------------------------------------

function write_to_file(str::S, fn) where {S<:AbstractString}
    open(fn, "w") do fout
        write(fout, str)
    end 
end

write_to_file(lines::Vector{S}, fn) where {S<:AbstractString} = write_to_file(⦿(lines), fn)

⇶(lines, fn) = write_to_file(lines, "$fn")


## -----------------------------------------------
# two modes of communicating with external program

split_outp(outp::String, delims::Vector) = String[ st for st ∈ split(outp,delims,keepempty=false) ]

function submit_script_lines(
    scripts::Vector{S}, 
    program::Cmd;
    delims = ['\n',]
    ) where {S<:AbstractString}
    @inline split_output(o) = split_outp(o, delims)
    ISO = open(program, read=true, write=true)
    for s ∈ scripts
        write(ISO, "$s\n")
    end
    res = String(Char.(read(ISO))) |> split_output
    if length(res) == 1 && (occursin("Error", res[1]) || occursin("error", res[1]))
        @warn "submit_script_lines() error : " * res[1]
        return String[]
    else
        return res
    end 
end


function run_program_with_file(
    fn::AbstractString, 
    program::Cmd;
    delims = ['\n',]
    )
    @assert isfile(fn)
    res = readlines(`$program $fn`)
    if length(res) <=3 || any([occursin("Error",r) for r in res]) || any([occursin("error",r) for r in res])
        println.(res)
        @warn "run_program_with_file() error."
        return String[]
    else
        return res
    end 
end


function run_program_with_tmp_in(
    scripts::Vector{S}, 
    program::Cmd;
    delims = ['\n',]
    ) where {S<:AbstractString}

    tmp_f = "/tmp/tmp_isosuite_$(rand(10000:90000)).in"

    write_to_file(trim_comments_pound(scripts), tmp_f)

    run_program_with_file(tmp_f, program; delims=delims)

end


function run_program_with_pipeline(
    fn::AbstractString, 
    program::Cmd
    )
    @assert isfile(fn)
    res = readlines(pipeline(`cat $fn`,program))
    if length(res) <=3 || any([occursin("Error",r) for r in res]) || any([occursin("error",r) for r in res])
        println.(res)
        @warn "run_program_with_file() error."
        return String[]
    else
        return res
    end 
end


function run_program_with_pipeline_tmp_in(
    scripts::Vector{S}, 
    program::Cmd
    ) where {S<:AbstractString}
    
    tmp_f = "/tmp/tmp_isosuite_$(rand(10000:90000)).in"

    write_to_file(trim_comments_pound(scripts), tmp_f)

    run_program_with_pipeline(tmp_f, program)

end



## -----------------------------------------------
# CALL PROGRAMS

@inline remove_welcome(outp::Vector, n) = outp[n:end]

# submit scripts to the program `iso`
function iso(scripts::Vector{S}) where {S<:AbstractString}
    ENV["ISODATA"] = adslsh(ENV["ISOSUITE_FOLDER"])
    res = submit_script_lines(scripts, `$(adslsh(ENV["ISOSUITE_FOLDER"]))./iso`, delims=['\n','*'])
    return remove_welcome(res, 7)
end

# submit scripts to the program `smodes`
function smodes(scripts::Vector{S}) where {S<:AbstractString}
    ENV["ISODATA"] = adslsh(ENV["ISOSUITE_FOLDER"])
    res = submit_script_lines( scripts, `$(adslsh(ENV["ISOSUITE_FOLDER"]))./smodes` )
    return res
end

# submit scripts to the program `findsym`
function findsym(scripts::Vector{S}) where {S<:AbstractString}
    ENV["ISODATA"] = adslsh(ENV["ISOSUITE_FOLDER"])
    res = run_program_with_tmp_in(scripts, `$(adslsh(ENV["ISOSUITE_FOLDER"]))./findsym`)
    return remove_welcome(res,4)
end

function findsym_cifinput(fn_cif::AbstractString)
    ENV["ISODATA"] = adslsh(ENV["ISOSUITE_FOLDER"])
    res = run_program_with_file(fn_cif, `$(adslsh(ENV["ISOSUITE_FOLDER"]))./findsym_cifinput`)
    return res
end

function comsubs(scripts::Vector{S}) where {S<:AbstractString}
    ENV["ISODATA"] = adslsh(ENV["ISOSUITE_FOLDER"])
    res = run_program_with_pipeline_tmp_in(scripts, `$(adslsh(ENV["ISOSUITE_FOLDER"]))./comsubs`)
    return res
end

