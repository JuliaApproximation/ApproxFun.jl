using ApproxFun

#The following solves the Airy ODE with dirichlet boundary conditions

	x=Fun(x->x,[-1000.,15.])   # Fun corresponding to multiplication by x, on [-100,15]
    d=domain(x)                # The domain
    D=Derivative(d)            # The derivative operator
    B=dirichlet(d)             # Dirichlet boundary conditions, [u(-100),u(15)]


    #Construct operator

    A=[B;D^2-x]                # This is dirichlet conditions and u'' - x u
    b=[airyai(d.a),0.]         # We want it to equal airyai(-100) at -100, and 0 at
                               # 10, with 0 rhs
    #Solve ODE

    u=A\b                      # u satisfies A*u = b, or in other words,
                               # B*u = [airyai(-100.),0.] and (D2 - x)*u = 0.

    # Check the accuracy
    norm(u - Fun(airyai,d))


## We now solve with Neumann conditions

    B=neumann(d)
    A=[B;D^2-x]
    b=[airyaiprime(d.a),0.,0.]

    u=A\b


    # Check the accuracy
    norm(u - Fun(airyai,d))
