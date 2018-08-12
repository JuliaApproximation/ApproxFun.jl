## Barycentric formula

export bary,barysum

function bary(v::AbstractVector{Float64},pts::AbstractVector{Float64},x::Float64)
  n=length(v)
  @assert n == length(pts)


  retd = 1/(2*(x-pts[1]))
  retn = v[1]*retd

  for i = 2:2:n-1
    @inbounds cd = 1/(x-pts[i])
    retd -= cd
    @inbounds retn -= v[i]*cd
  end

  for i = 3:2:n-1
    @inbounds cd = 1/(x-pts[i])
    retd += cd
    @inbounds retn += v[i]*cd
  end

  cd = .5*(-1.)^(n-1)/(x-pts[n])
  retd += cd
  retn += v[n]*cd

  retn/retd
end


bary(v::AbstractVector{Float64},x::Float64)=bary(v,chebyshevpoints(length(v);kind=2),x)

function randomadaptivebary(f::Function)
  r=rand()
  fr=f(r)
  err=1.
  logn=1

  tol=200eps()

  vals=Float64[]

  while err > tol
      n=2^logn + 1
      pts=chebyshevpoints(n;kind=2)
      vals=f(pts)
      err=abs(bary(vals,pts,r)-fr)/n
      logn+=1
  end

  vals
end
