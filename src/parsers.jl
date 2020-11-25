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

parse_n_float64(x,n) = (try parse.(Float64, split(strip(x)," ",keepempty=false)) catch _ Vector([NaN for i=1:n]) end)

parse_3_float64(x) = parse_n_float64(x,3)

parse_6f(x) = parse_n_float64(x,6)