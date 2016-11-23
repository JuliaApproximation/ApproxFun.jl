using ApproxFun

 # Semicircle law
x=Fun(identity,-2..2)
w=sqrt(4-x^2)/(2π)

lanczos(w,5)


# Marchenco–Pastur law

r = .5
lmax = (1+sqrt(r))^2
lmin = (1-sqrt(r))^2

x= Fun(identity,lmin..lmax)

w=sqrt((lmax-x)*(x-lmin))/(π*x)

lanczos(w,5)


# Wachter law

a = 5; b= 10;
c = sqrt(a/(a+b)*(1-1/(a+b)))
d = sqrt(1/(a+b)*(1-a/(a+b)))

lmax = (c+d)^2
lmin = (c-d)^2

x = Fun(identity,lmin..lmax)

w = (a+b)*sqrt((x-lmin).*(lmax-x))/(2π*x*(1-x))

lanczos(w,5)
