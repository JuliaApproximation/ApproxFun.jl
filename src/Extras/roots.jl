## Root finding for Chebyshev expansions
#
#  Contains code that is based in part on Chebfun v5's chebfun/@chebteck/roots.m,
# which is distributed with the following license:

# Copyright (c) 2015, The Chancellor, Masters and Scholars of the University
# of Oxford, and the Chebfun Developers. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Oxford nor the names of its
#       contributors may be used to endorse or promote products derived from
#       this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function complexroots{C<:Chebyshev}(f::Fun{C})
    if length(f)==0 || (length(f)==1 && isapprox(f.coefficients[1],0))
        warn("Tried to take roots of a zero function.  Returning [].")
        Complex128[]
    elseif length(f)==1
        Complex128[]
    elseif length(f)==2
        Complex128[-f.coefficients[1]/f.coefficients[2]]
    else
        fromcanonical(f,colleague_eigvals(f.coefficients))
    end
end

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

function roots(f::Fun)
    f2=Fun(f,domain(f)) # default is to convert to Chebyshev/Fourier
    if space(f2)==space(f)
        error("roots not implemented for "*string(typeof(f)))
    else
        roots(f2)
    end
end


function roots{C<:Chebyshev}( f::Fun{C} )
# FIND THE ROOTS OF AN IFUN.

    d = domain(f)
    c = f.coefficients
    vscale = maxabs(values(f))
    if vscale == 0
        warn("Tried to take roots of a zero function.  Returning [].")
        return eltype(domain(f))[]
    end

    hscale = maximum( [abs(first(d)), abs(last(d))] )
    htol = eps(2000.)*max(hscale, 1)  # TODO: choose tolerance better

    if eltype(f) == BigFloat
        r = rootsunit_coeffs(convert(Vector{Float64},c./vscale), Float64(htol))
        # Map roots from [-1,1] to domain of f:
        rts = fromcanonical(d,r)
        fp = differentiate(f)
        while norm(f(rts)) > 1000eps(eltype(f))
            rts .-=f(rts)./fp(rts)
        end
    elseif eltype(f) == Complex{BigFloat}
        r = rootsunit_coeffs(convert(Vector{Complex{Float64}},c./vscale), Float64(htol))
        # Map roots from [-1,1] to domain of f:
        rts = fromcanonical(d,r)
        fp = differentiate(f)
        while norm(f(rts)) > 1000eps(eltype(f))
            rts .-=f(rts)./fp(rts)
        end
    else
        cvscale=c./vscale
        r = rootsunit_coeffs(cvscale, Float64(htol))

        # Check endpoints, as these are prone to inaccuracy
        # which can be deadly.
        if (isempty(r) || !isapprox(last(r),1.)) && abs(sum(cvscale)) < htol
            push!(r,1.)
        end
        if (isempty(r) || !isapprox(first(r),-1.)) && abs(alternatingsum(cvscale)) < htol
            insert!(r,1,-1.)
        end


        # Map roots from [-1,1] to domain of f:
        rts = fromcanonical(d,r)
    end
    return rts
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


function colleague_balance!(M)
    n=size(M,1)
    if n<2
        return M
    end

    d=1+abs(M[1,2])
    M[1,2]/=d;M[2,1]*=d

    for k=3:n
        M[k-1,k]*=d;M[k,k-1]/=d
        d+=abs(M[1,k])
        M[1,k]/=d;M[k-1,k]/=d;M[k,k-1]*=d
    end
    M
end


# colleague_eigvals( c::Vector )=hesseneigvals(colleague_balance!(colleague_matrix(c)))

colleague_eigvals( c::Vector )=eigvals(colleague_matrix(c))

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

rootsunit_coeffs{T<:Number}(c::Vector{T}, htol::Float64)=rootsunit_coeffs(c,htol,ClenshawPlan(T,Chebyshev(),length(c),length(c)))
function rootsunit_coeffs{S,T<:Number}(c::Vector{T}, htol::Float64,clplan::ClenshawPlan{S,T})
# Computes the roots of the polynomial given by the coefficients c on the unit interval.


    # If recursive subdivision is used, then subdivide [-1,1] into [-1,splitPoint] and [splitPoint,1].
    const splitPoint = -0.004849834917525

    # Simplify the coefficients by chopping off the tail:
    nrmc=norm(c, 1)
    @assert nrmc > 0  # There's an error if we make it this far
    c = chop(c,eps()*nrmc)
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

    elseif n <= 70

        # COLLEAGUE MATRIX
        # The recursion subdividing below will keep going until we have a piecewise polynomial
        # of degree at most 50 and we likely end up here for each piece.

        # Adjust the coefficients for the colleague matrix
        # Prune roots depending on preferences:
        r = PruneOptions( colleague_eigvals(c), htol )::Vector{Float64}
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

function extremal_args(f::Fun)
    d=domain(f)
    if isa(d,PeriodicInterval)
        roots(differentiate(f))
    elseif isa(d,PeriodicDomain)  # avoid complex domains
        S=typeof(space(f))
        fromcanonical(f,extremal_args(Fun(f.coefficients,S(canonicaldomain(f)))))
    else
        dab=∂(domain(f))
        if length(f) <=2 #TODO this is only relevant for Polynomial bases
            dab
        else
            [dab;roots(differentiate(f))]
        end
    end
end

for op in (:(Base.maximum),:(Base.minimum),:(Base.extrema),:(Base.maxabs),:(Base.minabs))
    @eval begin
        function $op{S<:RealSpace,T<:Real}(f::Fun{S,T})
            pts = extremal_args(f)

            $op(f(pts))
        end
    end
end

for op in (:(Base.maxabs),:(Base.minabs))
    @eval begin
        function $op{S,T}(f::Fun{S,T})
            # complex spaces/types can have different extrema
            pts = extremal_args(abs(f))

            $op(f(pts))
        end
    end
end



for op in (:(Base.indmax),:(Base.indmin))
    @eval begin
        function $op{S<:RealSpace,T<:Real}(f::Fun{S,T})
            # the following avoids warning when differentiate(f)==0
            pts = extremal_args(f)
            # the extra real avoids issues with complex round-off
            pts[$op(real(f(pts)))]
        end

        function $op{S,T}(f::Fun{S,T})
            # the following avoids warning when differentiate(f)==0
            pts = extremal_args(f)
            fp=f(pts)
            @assert norm(imag(fp))<100eps()
            pts[$op(real(fp))]
        end
    end
end

for op in (:(Base.findmax),:(Base.findmin))
    @eval begin
        function $op{S<:RealSpace,T<:Real}(f::Fun{S,T})
            # the following avoids warning when differentiate(f)==0
            pts = extremal_args(f)
            ext,ind = $op(f(pts))
	    ext,pts[ind]
        end
    end
end



## Root finding for Laurent expansion


function companion_matrix{T}(c::Vector{T})
    n=length(c)-1

    if n==0
        zeros(T,0,0)
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
    function complexroots(coefficients::Vector)
        c=chop(coefficients,10eps())

        # Only use special routine for large roots
        if length(c)≥70
            Main.AMVW.rootsAMVW(c)
        else
            hesseneigvals(companion_matrix(c))
        end
    end
else
    complexroots(coefficients::Vector)=hesseneigvals(companion_matrix(chop(coefficients,10eps())))
end

complexroots(neg::Vector,pos::Vector)=complexroots([flipdim(chop(neg,10eps()),1);pos])
complexroots{DD}(f::Fun{Laurent{DD}})=mappoint(Circle(),domain(f),complexroots(f.coefficients[2:2:end],f.coefficients[1:2:end]))
complexroots{DD}(f::Fun{Taylor{DD}})=mappoint(Circle(),domain(f),complexroots(f.coefficients))



function roots{DD}(f::Fun{Laurent{DD}})
    irts=filter!(z->in(z,Circle()),complexroots(Fun(f.coefficients,Laurent(Circle()))))
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


roots{D}(f::Fun{Fourier{D}})=roots(Fun(f,Laurent))

function roots{P<:PiecewiseSpace}(f::Fun{P})
    rts=mapreduce(roots,vcat,vec(f))
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


## Root finding for JacobiWeight expansion

# Add endpoints for JacobiWeight
# TODO: what about cancellation?
function roots{S<:JacobiWeight,T}(f::Fun{S,T})
    sp=space(f)
    d=domain(sp)
    rts=roots(Fun(f.coefficients,sp.space))
    if sp.α > 0
        rts=[first(d);rts]
    end
    if sp.β > 0
        rts=[rts;last(d)]
    end
    rts
end
