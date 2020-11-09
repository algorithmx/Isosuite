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
