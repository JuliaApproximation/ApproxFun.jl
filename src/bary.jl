## Barycentric formula

export bary,barysum

function barysum(v,pt,x)
  n=length(v)  
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
  
  
function bary(v::Array{Float64,1},x::Float64)
  barysum(v,chebyshevpoints(length(v)),x)
end

function bary(v::Array{Float64,1},d,x::Float64)
  bary(v,tocanonical(d,x))
end

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
      err=abs(barysum(vals,pts,r)-fr)/n
      logn+=1
  end
  
  vals
end