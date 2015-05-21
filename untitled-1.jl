using ApproxFun,SO
    import ApproxFun:plot,diagop,columnspace
    setplotter("GLPlot")
d=Disk()
    # initial condition
    u0   = Fun((x,y)->exp(-50x^2-40(y-.1)^2)+.5exp(-30(x+.5)^2-40(y+.2)^2),d)
    B= [dirichlet(d) ,neumann(d)]

    L=-lap(d)^2

    h    = 0.0005
    bcs=[0.,0.]

    uin=[ProductFun(u0),ProductFun(u0)]
    op=L


    nt=size(uin[1],2)
    @time SBE  = discretize([B,I-h^2*op],domain(uin[1]),nt)            # backward euler for first 2 time steps
    @time SBDF = discretize([B,I-4.0/9.0*h^2*op],domain(uin[1]),nt)    # BDF formula for subsequent itme steps

    A=I-h^2*op
    C=I-4.0/9.0*h^2*op


u1,u2=uin
#    glp=plot(pad(u1,80,80))
    u3 =SBE\[bcs,2u2-u1]
    u4 =SBE\[bcs,2u3-u2]
    u5e =SBE\[bcs,2u4-u3]
    u5  = SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)]
#    u4,u3,u2,u1  = chop(SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],1000eps()),u4,u3,u2
#    plot(pad(u4,80,80),glp...)


plot([real(u5.coefficients[1]),real(u5e.coefficients[1])])


Gxin,Gyin,F=ApproxFun.pde_standardize_rhs(SBDF,[bcs,1/9.0*(24u4-22u3+8u2-u1)])
    F=ProductFun(F)
n = length(OS.Rdiags)
F=pad(F,size(F,1),n)
Gx=pad(coefficients(Gxin).',:,n)

rangespace(SBDF)

cfs=Fun(1/9.0*(24u4-22u3+8u2-u1)).coefficients
A,B=domainspace(SBDF),rangespace(SBDF)

g=ProductFun(Fun(cfs,A))
rcfs=Fun{typeof(columnspace(B,1)),eltype(cfs)}[Fun(g.coefficients[k],columnspace(B,k)) for k=1:length(g.coefficients)]

pf=ProductFun(rcfs,B)
f=Fun(pf)
pf2=ProductFun(f)

(pf.coefficients[1]        -pf2.coefficients[1]).coefficients |>norm
pf2.coefficients[1]           |>space
(pf-pf2 )|>ApproxFun.coefficients|>norm

pf.coefficients[1][0.1]-pf2.coefficients[1][0.1]

plot(pf2.coefficients[1])
pad(ApproxFun.coefficients(pf),100,100)-pad(ApproxFun.coefficientmatrix(f),100,100)|>norm



@code_typed ProductFun(f)



ProductFun(.coefficients[1]|>plot
ProductFun(Fun(Fun(ProductFun(rcfs,B)           ).coefficients,B)).coefficients[1]|>plot

ProductFun(rcfs,B).coefficients[1]|>plot


Fun(f.coefficients[1]     ,ApproxFun.columnspace(rangespace(SBDF),1))|>plot

@code_typed coefficients(Fun(f).coefficients,domainspace(SBDF),rangespace(SBDF)   )

ProductFun(1/9.0*(24u4-22u3+8u2-u1),rangespace(SBDF)   ).coefficients[1]|>plot
ProductFun(Fun(1/9.0*(24u4-22u3+8u2-u1),rangespace(SBDF)   )).coefficients[1]|>plot

ApproxFun.columnspace(rangespace(op),1)

Fun(1/9.0*(24u4-22u3+8u2-u1).coefficients[1],rangespace(op))|>ApproxFun.plot
F.coefficients[1]           |>ApproxFun.plot


OS=SBDF;
k=1;op=OS.Rdiags[k]
    f= F.coefficients[1]
    [OS.Bx[k];op]\Any[0.+0.im,0.+0.im,f]|>ApproxFun.plot

@code_typed pdesolve(SBDF,[bcs,1/9.0*(24u4-22u3+8u2-u1)])

SBDF.Bx[1][2]*u5.coefficients[1]
SBDF.Bx[1][2][1:10]
ApproxFun.plot([u5.coefficients[1]',u5e.coefficients[1]'])
ApproxFun.plot([u5.coefficients[1],u5e.coefficients[1]])

diagop(B,1)\(2u3-u2).coefficients[1]|>ApproxFun.plot

f=(2u3-u2).coefficients[1]
    [ldirichlet(space(f));lneumann(space(f));diagop(A,1)]\Any[0.,0.,f]|>ApproxFun.plot

f=1/9.0*(24u4-22u3+8u2-u1).coefficients[1]
    ([ldirichlet(space(f));lneumann(space(f));diagop(C,1)]\Any[0.,0.,f])'|>plot

f=(2u3-u2).coefficients[1]
    [ldirichlet(space(f));lneumann(space(f));diagop(C,1)]\Any[0.,0.,f]|>ApproxFun.plot





(2u3-u2).coefficients[1]|>ApproxFun.plot
(1/9.0*(24u4-22u3+8u2-u1)).coefficients[1]|>ApproxFun.plot
diagop(A,1)[1:10,1:10]|>chopm

diagop(B,1)[1:10,1:10]|>chopm


u4.coefficients[1]|>plot









diagop(h^2*op,1)[1:10,1:10]|>chopm
ApproxFun.columnspace(domainspace(op),1)
ApproxFun.columnspace(rangespace(op),1)


diagop(I-4.0/9.0*h^2*op,1)[1:10,1:10]|>chopm

I-4.0/9.0*h^2*op|>ApproxFun.introspect


glp=plot(u4)
u4,u3,u2,u1  = chop(SBDF\[bcs,1/9.0*(24u4-22u3+8u2-u1)],1000eps()),u4,u3,u2
glp=plot(u4)

u4=timeevolution(2,B,L,u0,h,10);


# domain is a disk
d=Disk()
# initial condition
u0   = Fun((x,y)->exp(-40(x-.1)^2-40(y+.2)^2),d)
# dirichlet boundary conditions
B=dirichlet(d)

Δ=lap(d)
h=0.0001 # time step
timeevolution(B,Δ,u0,h)

op=Δ

bcs=[0.]
uin=ProductFun(u0)
glp=plot(pad(uin,80,80))
nt=size(uin,2)
d=domain(uin)
@time SBE   = discretize([B,I-h*op],d,nt)            # backward euler for first 2 time steps
@time SBDF2 = discretize([B,I-2.0/3.0*h*op],d,nt)    # BDF formula for subsequent itme steps
@time SBDF3 = discretize([B,I-6.0/11.0*h*op],d,nt)    # BDF formula for subsequent itme steps
@time SBDF4 = discretize([B,I-12.0/25.0*h*op],d,nt)    # BDF formula for subsequent itme steps

@time u1=uin
@time u2=SBE\[bcs,u1]
@time plot(pad(u2,80,80),glp...)#updates window
@time u3=SBDF2\[bcs,1/3.0*(4u2-u1)]
@time plot(pad(u3,80,80),glp...)#updates window
@time u4=SBDF3\[bcs,1/11.0*(18u3-9u2+2u1)]
@time plot(pad(u4,80,80),glp...)#updates window

for k=1:10
    @time u4,u3,u2,u1  = chop(SBDF4\[bcs,1/25.0*(48u4-36u3+16u2-3u1)],1000eps()),u4,u3,u2
    @show size(u4)
#    plot(pad(u4,80,80),glp...)#updates window
end

size(chop!(u4,10000eps()))