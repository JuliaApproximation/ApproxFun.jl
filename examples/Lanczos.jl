using ApproxFun

function lanczos(w,N)
    x = Fun(identity,domain(w))
    P = Array(IFun,N + 1)
    β = Array(Float64,N)
    γ = Array(Float64,N)

    P[1] = Fun([1./sqrt(sum(w))],domain(x))

    v = x.*P[1]
    β[1] = sum(w.*v.*P[1])

    v = v - β[1]*P[1]
    γ[1] = sqrt(sum(w.*v.^2))

    P[2] = v/γ[1]

    for k = 2:N
        v = x.*P[k] - γ[k-1]*P[k-1]
        β[k] = sum(w.*v.*P[k])
        v = v - β[k]*P[k]
        γ[k] = sqrt(sum(w.*v.^2))
        P[k+1] = v/γ[k]
    end
    
    β,γ
end



 # Semicircle law
x=Fun(x->x,[-2.,2.])
w=sqrt(4-x^2)/(2π)

lanczos(w,5)


# Marchenco–Pastur law

r = .5
lmax = (1+sqrt(r))^2
lmin = (1-sqrt(r))^2

x= Fun(identity,[lmin,lmax])

w=sqrt((lmax-x)*(x-lmin))/(π*x)

lanczos(w,5)


# Wachter law

a = 5; b= 10;
c = sqrt(a/(a+b)*(1-1/(a+b)))
d = sqrt(1/(a+b)*(1-a/(a+b)))

lmax = (c+d)^2
lmin = (c-d)^2

x = Fun(identity,[lmin,lmax])

w = (a+b)*sqrt((x-lmin).*(lmax-x))/(2π*x*(1-x))

lanczos(w,5)