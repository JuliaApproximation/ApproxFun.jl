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



#function roots(f::IFun)
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


function roots( f::IFun )
# main() for root finding. 

domain = f.domain
c = f.coefficients
v = ichebyshevtransform( flipud( c ) )
vscale = maxabs( v, 1 )
#hscale = maxabs(domain)   TODO: fix this. 
hscale = 1

# Trivial case for f constant: 
if length( c ) == 1
    # Return a root at centre of domain if zero function, otherwise no root:
    r = ( c[1] == 0 ) ? 0 : [] 
else
    # Get scaled coefficients for the recursive call:
    htol = 100*eps(Float64)*max(hscale, 1)
    r = rootsunit_coeffs(c, htol)
end
   # Map roots from [-1,1] to domain of f: 
    return fromcanonical(domain, r)
end


function ColleagueEigVals( c )
# INPUT Chebyshev coefficients. OUTPUT: roots (via eigvals of colleague matrix).
    n = length(c)
    c = -1/2 * c[1:n-1] / c[n]
    n = n - 1
    c[n-1] = c[n-1] + 0.5
    oh = 0.5 * ones(n-1)
    A = diagm(oh, 1) + diagm(oh, -1)
    A[n-1, n] = 1
    A[:, 1] = flipud(c)
    # Standard colleague (See [Good, 1961]):
    r = eigvals(A)
    return r
end

function PruneOptions(r, all, prune, htol)
# Sometimes roots are pruned out because we don't want complex ones for example.
# This command deals with pruning. 
    # Clean the roots up a bit:
    if ( all==0 )
        # Remove dangling imaginary parts:
        mask = abs(imag(r)) .< htol
        r = real( r[mask] )
        # Keep roots inside [-1 1]:
        r = sort( r[abs(r) .<= 1 + htol] )
        # Correct roots over ends:
        if ( ~isempty(r) )
            r[1] = max(r[1], -1)
            r[end] = min(r[end], 1)
        end
    elseif ( prune==1 )
        rho = sqrt(eps(Float64))^(-1/n)
        rho_roots = abs(r + sqrt(r.^2 - 1))
        rho_roots[rho_roots .< 1] = 1./rho_roots[rho_roots .< 1]
        r = r[rho_roots .<= rho]
    end
    return r
end

function rootsunit_coeffs(c, htol)
# Computes the roots of the polynomial given by the coefficients c on the unit interval.
    
    splitPoint = -0.004849834917525;
    splitInteger = 513;

    # TODO: How do we allow users to submit preferences to the roots command. 
    all = 0
    recurse = 1
    prune = 0
    
    n = length(c)
    # Simplify the coefficients:
    tailMmax = eps(Float64)*norm(c, 1)
    # Find the final coefficient about tailMax:
    idx = findfirst( flipud( abs(c[:]) .> tailMmax ) ) 
    n = n - idx + 1
    
    # Truncate the coefficients (rather than alias):
    if ( n > 1 && n < length(c) )
        c = c[1:n]; 
    end
    
if ( n == 0 )
    r = []
elseif ( n == 1 )
     # Return a root at centre of domain if zero function, otherwise no root:
    r = ( c[1] == 0 ) ? 0 : [] 
elseif ( n == 2 )
    # Is the root in [-1,1]?
    r = -c[1]/c[2];
    if ( all==0 )
        r = ( (abs(imag(r)) .> htol) | (r < -(1 + htol)) | (r > (1 + htol)) ) ? [] : max( min(real(r), 1), -1)
    end
    # Is n small enough for the roots to be calculated directly?
elseif ( n <= 50 )
        # The recursion subdividing below will keep going until we have a piecewise polynomial of degree at most 50 and we likely end up here for each piece. 
    
        # Adjust the coefficients for the colleague matrix:
        r = ColleagueEigVals(c);

        # Prune roots depending on preferences: 
        r = PruneOptions( r, all, prune, htol);

#  RECURSE SUBDIVISION: 
elseif ( n <= splitInteger )
     # For small n (n <= 513) it is faster (at least in MATLAB) to form evaluation matrices explicitly and keep them in memory rather than doing a Clenshaw evaluation every time. 
        
    # Have we assembled the matrices TLEFT and TRIGHT. If we have do not make them again: 
    if ( ~isdefined(:Tleft) )
        # TODO: This can all be done with chebyshevtransforms. Change. 

        # Create the coefficients for TLEFT using the FFT directly: 
        x = points([-1,splitPoint],splitInteger)
        x = x[:]
        Tleft = ones(splitInteger,splitInteger) 
        Tleft[:,2] = x
        for k = 3:splitInteger
                Tleft[:,k] = 2 * x .* Tleft[:,k-1] - Tleft[:,k-2]
        end
        Tleft = [ Tleft[splitInteger:-1:2,:] ; Tleft[1:(splitInteger-1),:] ]
        Tleft = real(fft(Tleft,1) / (splitInteger-1))
        Tleft = triu( [ 0.5*Tleft[1,:] ; Tleft[2:(splitInteger-1),:] ;0.5*Tleft[splitInteger,:] ] )
            
        # Create the coefficients for TRIGHT much in the same way:
            x = points([splitPoint,1],splitInteger)
            x = x[:];
            Tright = ones(splitInteger,splitInteger)
            Tright[:,2] = x
        for k = 3:513
                Tright[:,k] = 2 * x .* Tright[:,k-1] - Tright[:,k-2]
        end
        Tright = [ Tright[splitInteger:-1:2,:] ; Tright[1:(splitInteger-1),:] ]
        Tright = real(fft(Tright,1) / (splitInteger-1))
        Tright = triu( [ 0.5*Tright[1,:] ; Tright[2:(splitInteger-1),:] ; 0.5*Tright[splitInteger,:] ] )
    end
    # Compute the new coefficients:
        cleft = Tleft[1:n,1:n] * c
        cright = Tright[1:n,1:n] * c
        
    # Recurse:
    r = [ (splitPoint - 1)/2 + (splitPoint + 1)/2*rootsunit_coeffs(cleft, 2*htol) ;
         (splitPoint + 1)/2 + (1 - splitPoint)/2*rootsunit_coeffs(cright, 2*htol) ]
else 
    # In n is too large then don't form evaluation matrices (TLEFT and TRIGHT) explicitly, just do Clenshaw.  
    # Evaluate the polynomial on both intervals:
    pts = [ points([ -1, splitPoint ],n) ; points([ splitPoint, 1 ],n) ]
    v = clenshaw(c, pts[:])
    # Get the coefficients on the left:
    cleft = chebyshevtransform(v[1:n])

    # Get the coefficients on the right:
    cright = chebyshevtransform(v[n+1:2*n])
    
    # Recurse:
    r = [ (splitPoint - 1)/2 + (splitPoint + 1)/2*rootsunit_coeffs(cleft, 2*htol) ;
            (splitPoint + 1)/2 + (1 - splitPoint)/2*rootsunit_coeffs(cright, 2*htol) ]
end
    return r
end
