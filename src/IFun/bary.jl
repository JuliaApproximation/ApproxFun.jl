## Barycentric formula

export bary,barysum

function bary(v::Vector{Float64},pt::Vector{Float64},x::Float64)
  n=length(v)  
  @assert n == length(pt)
  vv = unsafe_view(v)
  
  pts = unsafe_view(pt)
  
  
  retd = .5/(x-pts[1])
  retn = v[1]*retd
  
  for i = 2:2:n-1
    cd = 1./(x-pts[i])
    retd -= cd
    retn -= vv[i]*cd
  end
  
  for i = 3:2:n-1
    cd = 1./(x-pts[i])    
    retd += cd
    retn += vv[i]*cd
  end

  cd = .5*(-1.)^(n-1)/(x-pts[n])  
  retd += cd
  retn += vv[n]*cd
  
  retn/retd
end
  

bary(v::Vector{Float64},x::Float64)=bary(v,chebyshevpoints(length(v)),x)
bary(v::Vector{Float64},d::IntervalDomain,x::Float64)=bary(v,tocanonical(d,x))

function randomadaptivebary(f::Function)
  r=rand()
  fr=f(r)
  err=1.
  logn=1
  
  tol=200eps()
  
  vals=Float64[]
  
  while err > tol
      n=2^logn + 1
      pts=chebyshevpoints(n)
      vals=f(pts)
      err=abs(bary(vals,pts,r)-fr)/n
      logn+=1
  end
  
  vals
end