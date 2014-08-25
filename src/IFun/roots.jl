## Root finding


function complexroots(cin::Vector)
    c=chop(cin,10eps())
    if c == [] || length(c) == 1
        return []
    elseif length(c) == 2
        return [-c[1]/c[2]]
    else 
        n=length(c)-1;
        
        I = [ones(Int64,n),2:n-1,2:n];
        J=[1:n,3:n,1:n-1];
        V = [-c[end-1]/(2c[end]),.5-c[end-2]/(2c[end]),-c[end-3:-1:1]/(2c[end]),.5*ones(n-2),.5*ones(n-2),1];
        A=full(sparse(I,J,V));
        
        return eigvals(A)
    end
end


function complexroots(f::IFun)
    fromcanonical(f,complexroots(f.coefficients))
end



function roots(f::IFun)
    irts=map(real,filter!(x->abs(x)<=1.+10eps(),filter!(isreal,complexroots(f.coefficients))))
    
    map!(x->x>1.?1.:x,irts)
    map!(x->x<-1.?-1.:x,irts)
        
    if length(irts)==0
        Float64[]
    else
        fromcanonical(f,irts)
    end
end


##TODO: allow using routines for complex domains below
function Base.maximum(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    maximum(f[pts])
end

function Base.minimum(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    minimum(f[pts])
end

function Base.indmax(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    pts[indmax(f[pts])]
end

function Base.indmin(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    pts[indmin(f[pts])]
end

