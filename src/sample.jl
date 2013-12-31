bisectioninv(f::IFun,x::Real) = first(bisectioninv(f,[x]))

bisectioninv(c::Vector{Float64},xl::Vector{Float64})=bisectioninv(c,xl,Array(Float64,length(xl)),Array(Float64,length(xl)),Array(Float64,length(xl)))


function bisectioninv(c::Vector{Float64},x::Float64)
    a = -1.;
    b = 1.;
    
#     a_v = unsafe_view(a);
#     b_v = unsafe_view(b);    
    
    for k=1:47  #TODO: decide 47
        m=.5*(a+b);
        val = clenshaw(c,m);

            (val<= x) ? (a = m) : (b = m)
    end
    .5*(a+b)    
end

function bisectioninv(c::Vector{Float64},xl::Vector{Float64},bk::Vector{Float64},bk1::Vector{Float64},bk2::Vector{Float64}) 
    n = length(xl);
    a = -ones(n);
    b = ones(n);
    
#     a_v = unsafe_view(a);
#     b_v = unsafe_view(b);    
    
    for k=1:47  #TODO: decide 47
        m=.5*(a+b);
        vals = clenshaw(c,m,bk,bk1,bk2);
        
        for j = 1:n
            (vals[j] <= xl[j]) ? (a[j] = m[j]) : (b[j] = m[j])
        end
    end
    m=.5*(a+b)    
end


## Sampling

function sample(f::IFun,n::Integer)
    cf = cumsum(f);
    cf = cf - cf[f.domain.a];
    cf = cf/cf[f.domain.b];
    
    fromcanonical(f,bisectioninv(cf.coefficients,rand(n)))
end


sample(f::IFun)=sample(f,1)[1]
