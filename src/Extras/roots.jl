## Root finding


complexroots{D<:IntervalSpace}(f::Fun{D})=fromcanonical(f,colleague_eigvals(f.coefficients))

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

function roots{S,T}(f::Fun{S,T})
    f2=Fun(f,domain(f)) # default is to convert to Chebyshev/Fourier
    if space(f2)==space(f)
        error("roots not implemented for "*string(typeof(f)))
    else
        roots(f2)
    end
end


function roots( f::Fun{Chebyshev} )
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



function colleague_matrix{T<:Number}( c::Vector{T} )
#TODO: This is command isn't typed correctly
# COMPUTE THE ROOTS OF A LOW DEGREE POLYNOMIAL BY USING THE COLLEAGUE MATRIX: 
    n = length(c) - 1
    A=zeros(T,n,n)
    
    for k=1:n-1
        A[k+1,k]=0.5
        A[k,k+1]=0.5        
    end
    for k=1:n
        A[1,end-k+1]-=0.5*c[k]/c[end]
    end
    A[n,n-1]=one(T)
    
    # TODO: can we speed things up because A is lower-Hessenberg
    # Standard colleague matrix (See [Good, 1961]):
    A
end


colleague_eigvals( c::Vector )=hesseneigvals(colleague_matrix(c))

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

rootsunit_coeffs{T<:Number}(c::Vector{T}, htol::Float64)=rootsunit_coeffs(c,htol,ClenshawPlan(T,length(c)))
function rootsunit_coeffs{T<:Number}(c::Vector{T}, htol::Float64,clplan::ClenshawPlan{T})
# Computes the roots of the polynomial given by the coefficients c on the unit interval.

    
    # If recursive subdivision is used, then subdivide [-1,1] into [-1,splitPoint] and [splitPoint,1]. 
    const splitPoint = -0.004849834917525
    
    # Simplify the coefficients by chopping off the tail:
    c = chop(c,eps()*norm(c, 1))
    n = length(c)

    if n == 0
        
        # EMPTY FUNCTION
        r = Float64[]

    elseif n == 1
        
        # CONSTANT FUNCTION 
        r = ( c[1] == 0.0 ) ? Float64[0.0] : Float64[]

    elseif n == 2

        # LINEAR POLYNOMIAL
        r1 = -c[1]/c[2];
        r = ( (abs(imag(r1))>htol) | (abs(real(r1))>(1+htol)) ) ? Float64[] : Float64[max(min(real(r1),1),-1)]
        
    elseif n <= 50

        # COLLEAGUE MATRIX
        # The recursion subdividing below will keep going until we have a piecewise polynomial 
        # of degree at most 50 and we likely end up here for each piece. 
    
        # Adjust the coefficients for the colleague matrix
        # Prune roots depending on preferences: 
        r = PruneOptions( colleague_eigvals( c ), htol )::Vector{Float64}

    else 

        #  RECURSIVE SUBDIVISION: 
        # If n > 50, then split the interval [-1,1] into [-1,splitPoint], [splitPoint,1]. 
        # Find the roots of the polynomial on each piece and then concatenate. Recurse if necessary.  
        
        # Evaluate the polynomial at Chebyshev grids on both intervals:
        #(clenshaw! overwrites points, which only makes sence if c is real)

        v1 = isa(c,Vector{Float64})?clenshaw!( c, points([-1,splitPoint],n),clplan):clenshaw( c, points([-1,splitPoint],n),clplan)  
        v2 = isa(c,Vector{Float64})?clenshaw!( c, points([splitPoint,1] ,n),clplan):clenshaw( c, points([splitPoint,1] ,n),clplan) 
        
        # Recurse (and map roots back to original interval):
        p = plan_chebyshevtransform( v1 )
        r = [ (splitPoint - 1)/2 + (splitPoint + 1)/2*rootsunit_coeffs( chebyshevtransform(v1,p), 2*htol,clplan) ; 
                 (splitPoint + 1)/2 + (1 - splitPoint)/2*rootsunit_coeffs( chebyshevtransform(v2,p), 2*htol,clplan) ]

    end
    
    # Return the roots: 
    r
end


function extremal_args{S<:PiecewiseSpace}(f::Fun{S})
    return cat(1,[extremal_args(fp) for fp in vec(f)]...)
end

function extremal_args{S<:IntervalSpace,T}(f::Fun{S,T})
    d=domain(f)
    da,db=first(d),last(d)
    if length(f) <=2 #TODO this is only relevant for Polynomial bases 
        pts=[da,db]
    else
        pts=cat(1, da, db, roots(differentiate(f)))
    end
    return pts
end 

extremal_args{S<:PeriodicSpace,T}(f::Fun{S,T})=roots(differentiate(f))

for op in (:(Base.maximum),:(Base.minimum),:(Base.extrema),:(Base.maxabs),:(Base.minabs))
    @eval begin
        function $op{S<:RealSpace,T<:Real}(f::Fun{S,T})
            pts = extremal_args(f)
                
            $op(f[pts])
        end
    end
end

for op in (:(Base.maxabs),:(Base.minabs))
    @eval begin
        function $op{S,T}(f::Fun{S,T})
            # complex spaces/types can have different extrema
            pts = extremal_args(abs(f))
                
            $op(f[pts])
        end        
    end
end



for op in (:(Base.indmax),:(Base.indmin))
    @eval begin
        function $op{S<:RealSpace,T<:Real}(f::Fun{S,T})
            # the following avoids warning when differentiate(f)==0
            pts = extremal_args(f)
            pts[$op(f[pts])]
        end
    end
end

for op in (:(Base.findmax),:(Base.findmin))
    @eval begin
        function $op{S<:RealSpace,T<:Real}(f::Fun{S,T})
            # the following avoids warning when differentiate(f)==0
            pts = extremal_args(f)
            ext,ind = $op(f[pts])
	    ext,pts[ind]
        end
    end
end




## Fourier



## Root finding


function companion_matrix{T}(c::Vector{T})
    n=length(c)-1
    
    if n==0
        T[]
    else
        A=zeros(T,n,n)
        for k=1:n
            A[k,end]=-c[k]/c[end]
        end
        for k=2:n
            A[k,k-1]=one(T)
        end
        A
    end
end


if isdir(Pkg.dir("AMVW"))
    require("AMVW")
    function complexroots{T<:Number}(coefficients::ShiftVector{T})
        c=chop(coefficients.vector,10eps())
        
        # Only use special routine for large roots
        if length(c)≥450 || (isa(eltype(c),Complex) && length(c)≥250)
            Main.AMVW.rootsAMVW(c)
        else
            hesseneigvals(companion_matrix(chop(coefficients.vector,10eps())))
        end
    end
else
    complexroots{T<:Number}(coefficients::ShiftVector{T})=hesseneigvals(companion_matrix(chop(coefficients.vector,10eps())))
end

complexroots(f::Fun{Laurent})=mappoint(Circle(),domain(f),complexroots(deinterlace(f.coefficients)))
complexroots(f::Fun{Fourier})=complexroots(Fun(f,Laurent))

roots{S<:MappedSpace}(f::Fun{S})=fromcanonical(f,roots(Fun(coefficients(f),space(f).space)))

function roots(f::Fun{Laurent})
    irts=filter!(z->abs(abs(z)-1.)/length(f) < 100eps(),complexroots(Fun(f.coefficients,Laurent(Circle()))))
    if length(irts)==0
        Complex{Float64}[]
    else
        rts=fromcanonical(f,tocanonical(Circle(),irts))
        if isa(domain(f),PeriodicInterval)
            sort!(real(rts))  # Make type safe?
        else
            rts
        end
    end
end


roots(f::Fun{Fourier})=roots(Fun(f,Laurent))

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