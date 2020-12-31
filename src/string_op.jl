SPLTX(sss::AbstractString, X)  = string.(strip.(split(sss,X,keepempty=false)))
SPLTFORT(sss::AbstractString)  = string.(strip.(split(sss,[' ' , ',' , ':' , ';' , '\n'],keepempty=false)))
SPLTX1(sss::AbstractString, X)  = string.(strip.(split(sss,X,limit=2,keepempty=false)))
SPLTFORT1(sss::AbstractString)  = string.(strip.(split(sss,[' ' , ',' , ':' , ';' , '\n'],limit=2,keepempty=false)))
SPLTN(sss::AbstractString)  = SPLTX(sss,"\n")
SPLTS(sss::AbstractString)  = SPLTX(sss," ")
SPLTS1(sss::AbstractString) = SPLTX1(sss," ")
SPLTC(sss::AbstractString)  = SPLTX(sss,":")
SPLTD(sss::AbstractString)  = SPLTX(sss,".")
SPLTA(sss::AbstractString)  = SPLTX(sss,",")
SPLTEQ(sss::AbstractString)  = SPLTX(sss,"=")

STRPRM(lines) = String[strip(string(l)) for l in lines if length(strip(string(l)))>0]
trim_comments_pound(lines) = String[l for l in STRPRM(lines) if !startswith(l,"#")]
⦿(sss::Vector{S}) where {S<:AbstractString} = join(sss,'\n')

@inline decif(x) = first(SPLTD(x))