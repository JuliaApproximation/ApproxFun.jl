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


complexroots{D<:IntervalDomainSpace}(f::Fun{D})=fromcanonical(f,complexroots(f.coefficients))

#function roots(f::Fun)
#    irts=map(real,filter!(x->abs(x)<=1.+10eps(),filter!#(isreal,complexroots(f.coefficients))))
#    
#    map!(x->x>1.?1.:x,irts)
#    map!(x->x<-1.?-1.:x,irts)
#        
#    if length(irts)==0
#        Float64[]
#    else
#        fromcanonical(f,irts)
#    end
#end


function roots( f::Fun{ChebyshevSpace} )
# FIND THE ROOTS OF AN IFUN.  

    d = domain(f)
    c = f.coefficients
    vscale = maxabs( values( f ) )
    if vscale == 0
        warn("Tried to take roots of a zero function.  Returning [].")
        ##TODO: could be complex domain, in which case type should be Complex{Float64}
        return Float64[]
    end    
    
    hscale = maximum( [first(d), last(d)] ) 
    htol = eps(2000.)*max(hscale, 1)  # TODO: choose tolerance better
    r = rootsunit_coeffs(c./vscale, htol)
    # Map roots from [-1,1] to domain of f: 
    return fromcanonical(d, r)
end


function ColleagueEigVals( c::Array{Float64,1} )
#TODO: This is command isn't typed correctly
# COMPUTE THE ROOTS OF A LOW DEGREE POLYNOMIAL BY USING THE COLLEAGUE MATRIX: 
    n = length(c) - 1
    c = -1/2 * c[1:n] / c[n+1]
    c[n-1] += 0.5
    oh = 0.5 * ones(n-1)
    A = diagm(oh, 1) + diagm(oh, -1)
    A[n-1, n] = 1.0
    A[:, 1] = flipud(c)
    # TODO: can we speed things up because A is upper-Hessenberg
    # Standard colleague matrix (See [Good, 1961]):
    return eigvals( A )
end

function PruneOptions( r, htol::Float64 )
# ONLY KEEP ROOTS IN THE INTERVAL
    
    # Remove dangling imaginary parts:
    r = real( r[ abs(imag(r)) .< htol ] )
    # Keep roots inside [-1 1]:
    r = sort( r[ abs(r) .<= 1+htol ] )
    # Put roots near ends onto the domain:
    r = min( max( r, -1 ), 1)
    
    # Return the pruned roots: 
    return r
end

rootsunit_coeffs(c::Array{Float64,1}, htol::Float64)=rootsunit_coeffs(c,htol,ClenshawPlan(Float64,length(c)))
function rootsunit_coeffs(c::Array{Float64,1}, htol::Float64,clplan::ClenshawPlan{Float64})
# Computes the roots of the polynomial given by the coefficients c on the unit interval.

    
    # If recursive subdivision is used, then subdivide [-1,1] into [-1,splitPoint] and [splitPoint,1]. 
    const splitPoint = -0.004849834917525;
    
    # Simplify the coefficients by chopping off the tail:
    c = chop(c,eps()*norm(c, 1))
    n = length(c)

    if n == 0
        
        # EMPTY FUNCTION
        r = Float64[]

    elseif n == 1
        
        # CONSTANT FUNCTION 
        r = ( c[1] == 0.0 ) ? [0.0] : Float64[]

    elseif n == 2

        # LINEAR POLYNOMIAL
        r = -c[1]/c[2];
        r = ( (abs(imag(r))>htol) | (abs(real(r))>(1+htol)) ) ? Float64[] : Float64[max(min(real(r),1),-1)]
        
    elseif n <= 50

        # COLLEAGUE MATRIX
        # The recursion subdividing below will keep going until we have a piecewise polynomial 
        # of degree at most 50 and we likely end up here for each piece. 
    
        # Adjust the coefficients for the colleague matrix:
        r = ColleagueEigVals( c )

        # Prune roots depending on preferences: 
        r = PruneOptions( r, htol )

    else 

        #  RECURSIVE SUBDIVISION: 
        # If n > 50, then split the interval [-1,1] into [-1,splitPoint], [splitPoint,1]. 
        # Find the roots of the polynomial on each piece and then concatenate. Recurse if necessary.  
        
        # Evaluate the polynomial at Chebyshev grids on both intervals:
        #(clenshaw! overwrites points)
        v1 = clenshaw!( c, points([-1,splitPoint],n),clplan) 
        v2 = clenshaw!( c, points([splitPoint,1] ,n),clplan) 
        
        # Recurse (and map roots back to original interval):
        p = plan_chebyshevtransform( v1 )
        r = [ (splitPoint - 1)/2 + (splitPoint + 1)/2*rootsunit_coeffs( chebyshevtransform(v1,p), 2*htol,clplan) ; 
                 (splitPoint + 1)/2 + (1 - splitPoint)/2*rootsunit_coeffs( chebyshevtransform(v2,p), 2*htol,clplan) ]

    end
    
    # Return the roots: 
    r
end



for op in (:(Base.maximum),:(Base.minimum),:(Base.extrema),:(Base.maxabs),:(Base.minabs))
    @eval begin
        function ($op)(f::Fun)
            d=domain(f)
            # the following avoids warning when diff(f)==0
            if length(f) <=2
                pts=[d.a,d.b]
            else
                pts=[d.a,d.b,roots(diff(f))]
            end
                
            $op(f[pts])
        end
    end
end


for op in (:(Base.indmax),:(Base.indmin))
    @eval begin
        function ($op)(f::Fun)
            # the following avoids warning when diff(f)==0
            d=domain(f)
            if length(f) <=2
                pts=[d.a,d.b]
            else
                pts=[d.a,d.b,roots(diff(f))]
            end
                
            pts[$op(f[pts])]
        end
    end
end

for op in (:(Base.findmax),:(Base.findmin))
    @eval begin
        function ($op)(f::Fun)
            # the following avoids warning when diff(f)==0
            d=domain(f)
            if length(f) <=2
                pts=[d.a,d.b]
            else
                pts=[d.a,d.b,roots(diff(f))]
            end
            ext,ind = $op(f[pts])
	    ext,pts[ind]
        end
    end
end




## Fourier



## Root finding

function complexroots{T<:Number}(coefficients::ShiftVector{T})
    c=chop(coefficients.vector,10eps())
    n=length(c)-1
    A=zeros(T,n,n)
    A[:,end]=c[1:end-1]/c[end]
    for k=2:n
        A[k,k-1]=1.
    end
    eigvals(A)
end

complexroots(f::Fun{LaurentSpace})=fromcanonical(f,tocanonical(Circle(),complexroots(deinterlace(f.coefficients))))


function roots(f::Fun{LaurentSpace})
    irts=filter!(x->abs(abs(x)-1)<=100.eps(),complexroots(deinterlace(f.coefficients)))
    if length(irts)==0
        Complex{Float64}[]
    else
        rts=fromcanonical(f,tocanonical(Circle(),irts))
        if isa(domain(f),PeriodicInterval)
            real(rts)  # Make type safe?
        else
            rts
        end
    end
end



function roots{P<:PiecewiseSpace}(f::Fun{P})
    rts=[map(roots,vec(f))...]
    k=1
    while k < length(rts)
        if isapprox(rts[k],rts[k+1])
            rts=rts[[1:k,k+2:end]]
        else
            k+=1
        end
    end
    
    rts
end