module mna

using Reduce
@force using Reduce.Algebra

export nodes, volts, analyze, transferFunction

"""
    nodes(netlist::Dict)::Array{Int64,1}

Takes a netlist dictionary and returns a list of all
nodes in the netlist.
"""
function nodes(netlist::Dict)::Array{Int64,1}
    return ( values(netlist) 
                 |> x -> reduce((n,p) ->  [n...,p[2:(end-1)]...], x, init=[])
                 |> unique
                 |> x -> filter((n) -> n != 0, x))
end

"""
    volts(netlist::Dict)::Array{Any,1}

Takes a netlist dictionary and returns all the voltages
in the netlist.
"""
function volts(netlist::Dict)::Array{Any,1}
    return ( keys(netlist)
           |> x -> filter( c -> in( first(netlist[c])
                                  , [:V :VCVS :VCCS :CCVS :CCCS] )
                         , x)
           |> x -> identity.(x) )
end

"""
    analyze(netlist::Dict)::Dict

Performs modified nodal analysis (mna) on the given
`netlist` and returns the vector containing node voltages
and branch currents, as well as the transfer function H(ω).
"""
function analyze(netlist::Dict, ω = nothing)

    Reduce.rounded(true);
    Algebra.scientific_notation(2,2);

    m = length(nodes(netlist));
    n = length(volts(netlist));
    M = RExpr(zeros(m+n , m+n)) |> parse |> Reduce.mat;
    y = RExpr(zeros(m+n , 1)) |> parse |> Reduce.mat;

    v = volts(netlist);

    expr = function(type, value)
        if type in [:R :Z :VCCS]
            return (1 / value)
        elseif type == :C
            return (im * (ω == nothing ? :s : ω) * value)
        elseif type == :L
            return (1 / (im * (ω == nothing ? :s : ω) * value))
        else
            return value
        end
    end

    for (id, component) in netlist
        type = first(component);
        nodes = component[2:(end-1)];
        value = expr(first(component), last(component));

        if in(type, [:R :G :L :C :Z :Y])
            idx = ( repeat([nodes], 2) 
                  |> x -> Iterators.product(x...) 
                  |> collect 
                  |> x -> reshape(x, (sum(size(x)),1))
                  |> x -> filter(y -> !in(0,y), x) );
            for (i,j) in idx
                M[i,j] = M[i,j] + (i == j ? value : -value);
            end
        elseif type == :I
            iv = ( identity.(zip(nodes, [-value,value])) 
                  |> x -> filter(y -> y[1] != 0, x) );
            for i in iv
                y[i[1]] = y[i[1]] + i[2];
            end
        elseif type == :V
            i_idx = m + findfirst(x -> x == id, v);
            y[i_idx] = y[i_idx] + value;
            m_idx = ( identity.(zip(nodes, [1,-1])) 
                    |> x -> filter(y -> y[1] != 0, x) );
            for m in m_idx
                M[m[1], i_idx] = M[m[1], i_idx] + m[2];
                M[i_idx, m[1]] = M[i_idx, m[1]] + m[2];
            end
        elseif type == :VCCS
            i_idx = m + findfirst(x -> x == id, v);
            b_idx = ( identity.(zip(nodes[1:2], [-1 1])) 
                    |> x -> filter(y -> y[1] != 0, x) );
            c_idx = ( identity.(zip(nodes[3:4], [1 -1])) 
                    |> x -> filter(y -> y[1] != 0, x) );

            M[i_idx, i_idx] = M[i_idx, i_idx] + value;
            for b in b_idx
                M[b[1], i_idx] = M[b[1], i_idx] + b[2];
            end
            for c in c_idx
                M[i_idx, c[1]] = M[i_idx, c[1]] + c[2];
            end
        elseif type == :VCVS
            i_idx = m + findfirst(x -> x == id, v);
            b_idx = ( identity.(zip(nodes[3:4], [1 -1])) 
                    |> x -> filter(y -> y[1] != 0, x) );
            c_idx = ( identity.(zip(nodes, [-value value 1 -1])) 
                    |> x -> filter(y -> y[1] != 0, x) );
            for b in b_idx
                M[b[1], i_idx] = m[b[1], i_idx] + b[2];
            end
            for c in c_idx
                M[i_idx, c[1]] = m[i_idx, c[1]] + c[2];
            end
        elseif type == :CCVS
            nothing
        elseif type == :CCCS
            nothing
        else
            nothing
        end
    end

    if ω == nothing
        x = RExpr(M \ y);
    else
        x = convert(Array{Complex, 2}, M) \ convert(Array{Complex, 1}, y);
    end

    return x
end

"""
    transferFunciton(::RExpr) => ::Function

Takes an RExpr with only one symbolic variable `s` and
returns a function that evaluates the small signal model
at the given angle frequency ω.
"""
function transferFunction(x::RExpr)
    return (ω) -> Algebra.sub(:(s = $ω), x) |> parse |> eval
end

end # module
