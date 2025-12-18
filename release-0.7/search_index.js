var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ApproxFun.jl-Documentation-1",
    "page": "Home",
    "title": "ApproxFun.jl Documentation",
    "category": "section",
    "text": "ApproxFun is a package for approximating and manipulating functions, and for solving differential and integral equations.  "
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "A basic approach of computational mathematics that ApproxFun exploits is expansion in a basis $ f(x) \\approx \\sum_{k=1}^n c_k \\psi_k(x) $ Some traditional examples of bases psi_1(x)psi_2(x)ldots areTaylor series: 1zz^2ldots\nFourier series (for periodic functions on 0..2π): 1sin x cos x sin 2 x ldots\nChebyshev series (for non-periodic functions on -1..1): 1xcos 2 hboxacos x cos 3 hboxacos x ldotsIn ApproxFun, functions are represented by a Fun with two components: space, which dictates the basis and coefficients which is a finite vector of coefficients.  Note that each Fun can have a different length vector of coefficients, allowing for approximation of many different functions to high accuracy.  The approximation by a Fun can be determined by a variety of methods:(1) Explicitly specifying the coefficients:julia> f = Fun(Taylor(),[1,2,3]) # Represents 1 + 2z + 3z^2\nFun(Taylor(🕒),[1.0,2.0,3.0])\n\njulia> f(1.0)\n6.0(2) Constructors take in a Function and adaptively determine the     number of coefficients.  For example,julia> Fun(exp)\nFun(Chebyshev(【-1.0,1.0】),[1.26607,1.13032,0.271495,0.0443368,0.00547424,0.000542926,4.49773e-5,3.19844e-6,1.99212e-7,1.10368e-8,5.5059e-10,2.49797e-11,1.03911e-12,3.99195e-14])determines that f can be approximated to roughly machine precision using 14 coefficients.  See Constructors for more information.(3) Manipulation of Funs give new Funs, where the number of coefficients is determined from the input.  The simplest example is addition, which for compatible bases is just padding the vectors to the same length and adding.  julia> a = Fun(cos,Chebyshev()); ncoefficients(a)\n13\n\njulia> b = Fun(x->cos(10cos(x^2)),Chebyshev()); ncoefficients(b)\n51\n\njulia> ncoefficients(a+b)\n51On the other hand, multiplication results in an approximation with more coefficients than either a or b, so that the result approximates the true a*b to roughly machine accuracy:julia> ncoefficients(a*b)\n63\n\njulia> a(0.1)*b(0.1) - (a*b)(0.1)\n1.1102230246251565e-16The example of multiplication highlights the importance of adaptivity: if with a fixed discretization size, operations like multiplication would lose accuracy when the true function is no longer resolved by the discretization.  More complicated examples are solving differential equations, where the coefficients of the solution can be determined adaptively, see Equations.ApproxFun supports a number of different spaces, as described in Spaces.  A key component of ApproxFun is support for interaction between different spaces.  This is crucial for efficient solution of differential equations, where linear operators are described as acting between different spaces, see Operators.  "
},

{
    "location": "index.html#Contents-1",
    "page": "Home",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"usage/constructors.md\",\n         \"usage/domains.md\",\n         \"usage/spaces.md\",\n         \"usage/operators.md\",\n         \"usage/equations.md\",\n         \"faq.md\",\n         \"library.md\"]"
},

{
    "location": "usage/domains.html#",
    "page": "Domains",
    "title": "Domains",
    "category": "page",
    "text": ""
},

{
    "location": "usage/domains.html#Domains-1",
    "page": "Domains",
    "title": "Domains",
    "category": "section",
    "text": "Domain is an abstract type whose subtypes represent oriented domains on which we wish to approximate functions.   Examples include Interval, Ray, Line and Arc.   Periodic domains include PeriodicInterval, PeriodicLine and Circle."
},

{
    "location": "usage/domains.html#Relationship-with-spaces-1",
    "page": "Domains",
    "title": "Relationship with spaces",
    "category": "section",
    "text": "Every domain d has a default space, constructed via Space(d).  For example, the default space for Interval() is Chebyshev(Interval()), which is efficient for representing smooth functions.  On the other hand, the default space for PeriodicInterval() is Fourier(Interval()), which uses trigonometric polynomials to approximate periodic functions.  "
},

{
    "location": "usage/domains.html#Manipulating-domains-1",
    "page": "Domains",
    "title": "Manipulating domains",
    "category": "section",
    "text": "Domains can be manipulated to make more complicated domains.  For example, you can take the union of an interval and a circleInterval() ∪ Circle(3,0.5)    # equivalent to union(Interval(),Circle(3,0.5))and the following creates a rectangle [0,1]^2:rect=Interval(0,1)^2Some other set operations are partially implemented:Interval(0,2) ∩ Interval() # returns Interval(0,1)"
},

{
    "location": "usage/spaces.html#",
    "page": "Spaces",
    "title": "Spaces",
    "category": "page",
    "text": ""
},

{
    "location": "usage/spaces.html#Spaces-1",
    "page": "Spaces",
    "title": "Spaces",
    "category": "section",
    "text": "A Space is an abstract type whose subtypes indicate which space a function lives in. This typically corresponds to the span of a (possibly infinite) basis."
},

{
    "location": "usage/spaces.html#Chebyshev-space-1",
    "page": "Spaces",
    "title": "Chebyshev space",
    "category": "section",
    "text": "The default space in ApproxFun is Chebyshev, which represents expansions in Chebyshev polynomials:f(x) = sum_k=0^infty f_k T_k(x)where T_k(x) = cos k rm acos x.   Note that there is an intrinsic link between Chebyshev and CosSpace:  g(theta) = f(cos theta) = sum_k=0^infty f_k cos k thetaIn other words:DocTestSetup = quote\n    using ApproxFun\nendjulia> f=Fun(exp,Chebyshev());\n\njulia> g=Fun(CosSpace(),f.coefficients); # specify the coefficients directly\n\njulia> f(cos(0.1))\n2.70473560723178\n\njulia> g(0.1)\n2.7047356072317794"
},

{
    "location": "usage/spaces.html#Ultraspherical-spaces-1",
    "page": "Spaces",
    "title": "Ultraspherical spaces",
    "category": "section",
    "text": "A key tool for solving differential equations are the ultraspherical spaces, which can be defined by the span of derivatives of Chebyshev polynomials.  Note that Ultraspherical(1) corresponds to the Chebyshev basis of the second kind: U_k(x) = sin (k+1) rm acos x over sin rm acos x.   The relationship with Chebyshev polynomials follows from trigonemetric identities: T_k(x) = k U_k-1(x).  Converting between ultraspherical polynomials (with integer orders) is extremely efficient: it requires O(n) operations, where n is the number of coefficients."
},

{
    "location": "usage/spaces.html#Fourier-and-Laurent-spaces-1",
    "page": "Spaces",
    "title": "Fourier and Laurent spaces",
    "category": "section",
    "text": "There are several different spaces to represent functions on periodic domains, which are typically a PeriodicInterval, Circle or PeriodicLine.  CosSpace represents expansion in cosine series:f(theta) = sum_k=0^infty f_k cos k thetaSinSpace represents expansion in sine series:f(theta) = sum_k=0^infty f_k sin (k+1) thetaTaylor represents expansion with only non-negative complex exponential terms:f(theta) = sum_k=0^infty f_k rm e^rm i k thetaHardy{false} represents expansion with only negative complex exponential terms:f(theta) = sum_k=0^infty f_k rm e^-rm i (k+1) thetaFourier represents functions that are sums of sines and cosines.  Note that if a function has the formf(theta) = f_0 + sum_k=1^infty f_k^rm c cos k theta + f_k^rm s sin kthetathen the coefficients of the resulting Fun are order as f_0f_1^rm sf_1^rm c. For example:julia> f = Fun(Fourier(),[1,2,3,4]);\n\njulia> f(0.1)\n4.979356652307978\n\njulia> 1 + 2sin(0.1) + 3cos(0.1) + 4sin(2*0.1)\n4.979356652307979Laurent represents functions that are sums of complex exponentials.  Note that if a function has the formf(theta) = sum_k=-infty^infty f_k rm e^rm i k thetathen the coefficients of the resulting Fun are order as f_0f_-1f_1. For example:julia> f = Fun(Laurent(),[1,2,3,4]);\n\njulia> f(0.1)\n9.895287137755096 - 0.694843906533417im\n\njulia> 1 + 2exp(-im*0.1) + 3exp(im*0.1) + 4exp(-2im*0.1)\n9.895287137755094 - 0.6948439065334167im"
},

{
    "location": "usage/spaces.html#Modifier-spaces-1",
    "page": "Spaces",
    "title": "Modifier spaces",
    "category": "section",
    "text": "Some spaces are built out of other spaces.  A simple example is JacobiWeight(β,α,space) which weights space, which is typically Chebyshev() or Jacobi(b,a),  by a Jacobi weight (1+x)^b*(1-x)^a.@meta  DocTestSetup = nothing"
},

{
    "location": "usage/spaces.html#Unset-space-1",
    "page": "Spaces",
    "title": "Unset space",
    "category": "section",
    "text": "UnsetSpace is a special space that is used as a stand in when a space has not yet been determined, particularly by operators.  "
},

{
    "location": "usage/constructors.html#",
    "page": "Constructors",
    "title": "Constructors",
    "category": "page",
    "text": ""
},

{
    "location": "usage/constructors.html#Constructors-1",
    "page": "Constructors",
    "title": "Constructors",
    "category": "section",
    "text": "Funs in ApproxFun are instances of Julia types with one field to store coefficients and another to describe the function space. Similarly, each function space has one field describing its domain, or another function space. Let\'s explore:DocTestSetup = quote\n    using ApproxFun\nendjulia> x = Fun(identity,-1..1);\n\njulia> f = exp(x);\n\njulia> g = f/sqrt(1-x^2);\n\njulia> space(f)   # Output is pretty version of Chebyshev(Interval(-1.0,1.0))\nChebyshev(【-1.0,1.0】)\n\njulia> space(g)   # Output is pretty version of  JacobiWeight(-0.5,-0.5,Interval(-1.0,1.0))\n(1-x^2)^-0.5[Chebyshev(【-1.0,1.0】)]The absolute value is another case where the space of the output is inferred from the operation:julia> f = Fun(x->cospi(5x),-1..1);\n\njulia> g = abs(f);\n\njulia> space(f)\nChebyshev(【-1.0,1.0】)\n\njulia> space(g)\nChebyshev(【-1.0,-0.9000000000000002】)⨄Chebyshev(【-0.9000000000000002,-0.6999999999999996】)⨄Chebyshev(【-0.6999999999999996,-0.5000000000000001】)⨄Chebyshev(【-0.5000000000000001,-0.30000000000000043】)⨄Chebyshev(【-0.30000000000000043,-0.09999999999999962】)⨄Chebyshev(【-0.09999999999999962,0.10000000000000053】)⨄Chebyshev(【0.10000000000000053,0.29999999999999966】)⨄Chebyshev(【0.29999999999999966,0.500000000000001】)⨄Chebyshev(【0.500000000000001,0.6999999999999998】)⨄Chebyshev(【0.6999999999999998,0.9000000000000006】)⨄Chebyshev(【0.9000000000000006,1.0】)"
},

{
    "location": "usage/constructors.html#Convenience-constructors-1",
    "page": "Constructors",
    "title": "Convenience constructors",
    "category": "section",
    "text": "The default space is Chebyshev, which can represent non-periodic functions on intervals.  Each Space type has a default domain: for Chebyshev this is -1..1, for Fourier and Laurent this is -π..π.  Thus the following are synonyms:Fun(exp,Chebyshev(Interval(-1,1)))\nFun(exp,Chebyshev(Interval()))\nFun(exp,Chebyshev(-1..1))\nFun(exp,Chebyshev())\nFun(exp,-1..1)\nFun(exp,Interval())\nFun(exp,Interval(-1,1))\nFun(exp)If a function is not specified, then it is taken to be identity.  Thus we have the following synonyms:x = Fun(identity,-1..1)\nx = Fun(-1..1)\nx = Fun(identity)\nx = Fun()"
},

{
    "location": "usage/constructors.html#Specifying-coefficients-explicitly-1",
    "page": "Constructors",
    "title": "Specifying coefficients explicitly",
    "category": "section",
    "text": "It is sometimes necessary to specify coefficients explicitly.  This is possible via specifying the space followed by a vector of coefficients:julia> f = Fun(Taylor(),[1,2,3]);  # represents 1 + 2z + 3z^2\n\njulia> f(0.1)\n1.23\n\njulia> 1+2*0.1+3*0.1^2\n1.23In higher dimensions, ApproxFun will sum products of the 1D basis functions. So if T_i(x) is the ith basis function, then a 2D function can be approximated as the following: f(x  y) = sum_i j c_ij  T_i(x)  T_j(y)The products will be ordered lexicographically by the degree of the polynomial, i.e. in the order T_0(x)  T_0(y)  T_0(x)  T_1(y)   T_1(x)  T_0(y)   T_0(x)  T_2(y)   T_1(x)  T_1(y)   T_2(x)  T_0(y)    . For example, if we are in the two dimensional CosSpace space and we have coefficients c_1 c_2 c_3, then $ f(x, y) = c_1 \\cos(0 x) \\cos(0 y) + c_2 \\cos(0 x) \\cos(1 y) + c_3 \\cos(1 x) \\cos(0 y). $This is illustrated in the following code:julia> f = Fun(CosSpace()^2, [1,2,3])\nFun(CosSpace(【0.0,6.283185307179586❫)⊗CosSpace(【0.0,6.283185307179586❫),[1.0,2.0,3.0])\n\njulia> f(1,2)\n1.7886132445101346\n\njulia> 1cos(0*1)*cos(0*2) + 2cos(0*1)*cos(1*2) + 3cos(1*1)*cos(0*2)\n1.7886132445101346"
},

{
    "location": "usage/constructors.html#Using-ApproxFun-for-“manual”-interpolation-1",
    "page": "Constructors",
    "title": "Using ApproxFun for “manual” interpolation",
    "category": "section",
    "text": "The ApproxFun package for Julia implements all of the necessary operations for Chebyshev interpolation and operations (like differentiation or integration) on Chebyshev interpolants.Normally, you give it a function f and a domain d, and construct the Chebyshev interpolant by fc = Fun(f, d). The ApproxFun package figures out the necessary number of Chebyshev points (i.e., the polynomial order) required to interpolate f to nearly machine precision, so that subsequent operations on fc can be viewed as \"exact\".However, in cases where the function to be interpolated is extremely expensive, and possibly even is evaluated by an external program, it is convenient to be able to decide on the desired Chebyshev order in advance, evaluate the function at those points \"manually\", and then construct the Chebyshev interpolant. However, this procedure isn\'t documented in the ApproxFun manual. In this notebook, we show how to do that for the example f(x) = exp(2x) on the domain 0..1.DocTestSetup = nothing"
},

{
    "location": "usage/operators.html#",
    "page": "Operators",
    "title": "Operators",
    "category": "page",
    "text": ""
},

{
    "location": "usage/operators.html#Operators-1",
    "page": "Operators",
    "title": "Operators",
    "category": "section",
    "text": "DocTestSetup = quote\n    using ApproxFun\nendLinear operators between two spaces in ApproxFun are represented by subtypes of Operator.  Every operator has a domainspace and rangespace. That is, if a Fun f has the space domainspace(op), thenop*f is a Fun with space rangespace(op).Note that the size of an operator is specified by the dimension of the domain and range space.  "
},

{
    "location": "usage/operators.html#Calculus-operators-1",
    "page": "Operators",
    "title": "Calculus operators",
    "category": "section",
    "text": "Differential and integral operators are perhaps the most useful type of operators in mathematics.  Consider the derivative operator on CosSpace:julia> D = Derivative(CosSpace())\nConcreteDerivative:CosSpace(【0.0,6.283185307179586❫)→SinSpace(【0.0,6.283185307179586❫)\n 0.0  -1.0\n       0.0  -2.0\n             0.0  -3.0\n                   0.0  -4.0\n                         0.0  -5.0\n                               0.0  -6.0\n                                     0.0  -7.0\n                                           0.0  -8.0\n                                                 0.0  -9.0\n                                                       0.0  ⋱\n                                                            ⋱\n\njulia> f = Fun(θ->cos(cos(θ)),CosSpace());\n\njulia> fp = D*f;\n\njulia> fp(0.1) ≈ f\'(0.1) ≈ sin(cos(0.1))*sin(0.1)\ntrueHere, we specified the domain space for the derivative operator, and it automatically determined the range space:DocTestSetup = quote\n    using ApproxFun\n    D = Derivative(CosSpace())\n    f = Fun(θ->cos(cos(θ)),CosSpace())\n    fp = D*f\nendjulia> rangespace(D) == space(fp) == SinSpace()\ntrueOperators can be identified with infinite-dimensional matrices, whose entries are given by the canonical bases in the domain and range space.  In this case, the relevant formula is D cos k theta = -k sin k theta That is, the (k,k+1)th entry is as follows:julia> k,j = 5,6;\n\njulia> ej = Fun(domainspace(D),[zeros(j-1);1]);\n\njulia> D[k,j] ≈ (D*ej).coefficients[k] ≈ -k\ntrueThe Chebyshev space has the property that its derivatives are given by ultraspherical spaces:julia> Derivative(Chebyshev())\nConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n 0.0  1.0                                           \n      0.0  2.0                                      \n           0.0  3.0                                 \n                0.0  4.0                            \n                     0.0  5.0                       \n                          0.0  6.0                  \n                               0.0  7.0             \n                                    0.0  8.0        \n                                         0.0  9.0   \n                                              0.0  ⋱\n                                                   ⋱"
},

{
    "location": "usage/operators.html#Functionals-1",
    "page": "Operators",
    "title": "Functionals",
    "category": "section",
    "text": "A particularly useful class of operators are _functionals_, which map from functions to scalar numbers.  These are represented by operators of size 1 × ∞: that is, infinite-dimensional analogues of row vectors.As an example, the evaluation functional f(0) on CosSpace has the form:julia> B = Evaluation(CosSpace(),0)\nConcreteEvaluation:CosSpace(【0.0,6.283185307179586❫)→ConstantSpace(Point(0))\n 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  ⋯\n\njulia> B*f ≈ f(0)\ntrueAs can be seen from the output, rangespace(B) is a ConstantSpace(Point(0)), a one-dimensional space used to represent scalars whose domain is a single point, 0.Closely related to functionals are operators with finite-dimensional range. For example, the Dirichlet operator represents the restriction of a space to its boundary. In the case, of Chebyshev(), this amounts to evaluation at the endpoints ±1:julia> B = Dirichlet(Chebyshev())\nConcreteDirichlet:Chebyshev(【-1.0,1.0】)→2-element ArraySpace:\n ConstantSpace(Point(-1.0))\n ConstantSpace(Point(1.0))\n 1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  ⋯\n 1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  ⋯\n\njulia> size(B)\n(2, ∞)\n\njulia> B*Fun(exp)\nFun(2-element ArraySpace:\n ConstantSpace(Point(-1.0))\n ConstantSpace(Point(1.0)) ,[0.367879, 2.71828])\n\njulia> B*Fun(exp) ≈ Fun([exp(-1),exp(1)])\ntrue"
},

{
    "location": "usage/operators.html#Algebraic-manipulation-of-operators-1",
    "page": "Operators",
    "title": "Algebraic manipulation of operators",
    "category": "section",
    "text": "Operators can be algebraically manipulated, provided that the domain and range spaces are compatible, or can be made compatible.   As a simple example, we can add the second derivative of a Fourier space to the identity operator:julia> D2 = Derivative(Fourier(),2)\nDerivativeWrapper:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)\n 0.0   0.0                                                      \n 0.0  -1.0   0.0                                                \n       0.0  -1.0   0.0                                          \n             0.0  -4.0   0.0                                    \n                   0.0  -4.0   0.0                              \n                         0.0  -9.0   0.0                        \n                               0.0  -9.0    0.0                 \n                                     0.0  -16.0    0.0          \n                                            0.0  -16.0    0.0   \n                                                   0.0  -25.0  ⋱\n                                                           ⋱   ⋱\n\njulia> D2 + I\nPlusOperator:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)\n 1.0  0.0                                                     \n 0.0  0.0  0.0                                                \n      0.0  0.0   0.0                                          \n           0.0  -3.0   0.0                                    \n                 0.0  -3.0   0.0                              \n                       0.0  -8.0   0.0                        \n                             0.0  -8.0    0.0                 \n                                   0.0  -15.0    0.0          \n                                          0.0  -15.0    0.0   \n                                                 0.0  -24.0  ⋱\n                                                         ⋱   ⋱When the domain and range space are not the same, the identity operator becomes a conversion operator.  That is, to represent D+I acting on the Chebyshev space, we would do the following:julia> D = Derivative(Chebyshev())\nConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n 0.0  1.0                                           \n      0.0  2.0                                      \n           0.0  3.0                                 \n                0.0  4.0                            \n                     0.0  5.0                       \n                          0.0  6.0                  \n                               0.0  7.0             \n                                    0.0  8.0        \n                                         0.0  9.0   \n                                              0.0  ⋱\n                                                   ⋱\n\njulia> C = Conversion(Chebyshev(),Ultraspherical(1))\nConcreteConversion:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n 1.0  0.0  -0.5                                             \n      0.5   0.0  -0.5                                       \n            0.5   0.0  -0.5                                 \n                  0.5   0.0  -0.5                           \n                        0.5   0.0  -0.5                     \n                              0.5   0.0  -0.5               \n                                    0.5   0.0  -0.5         \n                                          0.5   0.0  -0.5   \n                                                0.5   0.0  ⋱\n                                                      0.5  ⋱\n                                                           ⋱\n\njulia> D + C\nPlusOperator:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n 1.0  1.0  -0.5                                             \n      0.5   2.0  -0.5                                       \n            0.5   3.0  -0.5                                 \n                  0.5   4.0  -0.5                           \n                        0.5   5.0  -0.5                     \n                              0.5   6.0  -0.5               \n                                    0.5   7.0  -0.5         \n                                          0.5   8.0  -0.5   \n                                                0.5   9.0  ⋱\n                                                      0.5  ⋱\n                                                           ⋱ApproxFun can automatically determine the spaces, so if one writes D + I it will translate it to D + C.  Now consider the Fredholm integral operator of the second kind:L u = u + rm e^x int_-1^1 u(x) rm dxWe can construct this usingjulia> x = Fun();\n\njulia> Σ = DefiniteIntegral(Chebyshev())\nConcreteDefiniteIntegral:Chebyshev(【-1.0,1.0】)→ConstantSpace\n 2.0  0.0  -0.666667  0.0  -0.133333  0.0  -0.0571429  0.0  -0.031746  0.0  ⋯\n\njulia> L = I + exp(x)*Σ\nLowRankPertOperator:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)\n 3.53213     0.0  -0.844044     0.0  …  0.0  -0.0401926    0.0  ⋯\n 2.26064     1.0  -0.753545     0.0     0.0  -0.0358831    0.0  ⋱\n 0.542991    0.0   0.819003     0.0     0.0  -0.0086189    0.0  ⋱\n 0.0886737   0.0  -0.0295579    1.0     0.0  -0.00140752   0.0  ⋱\n 0.0109485   0.0  -0.00364949   0.0     0.0  -0.000173785  0.0  ⋱\n 0.00108585  0.0  -0.000361951  0.0  …  0.0  -1.72358e-5   0.0  ⋱\n 8.99546e-5  0.0  -2.99849e-5   0.0     0.0  -1.42785e-6   0.0  ⋱\n 6.39687e-6  0.0  -2.13229e-6   0.0     1.0  -1.01538e-7   0.0  ⋱\n 3.98425e-7  0.0  -1.32808e-7   0.0     0.0   1.0          0.0  ⋱\n 2.20735e-8  0.0  -7.35785e-9   0.0     0.0  -3.50374e-10  1.0  ⋱\n  ⋮           ⋱     ⋱            ⋱   …   ⋱     ⋱            ⋱   ⋱\n\njulia> u = cos(10x^2);\n\njulia> (L*u)(0.1)\n1.3777980523127336\n\njulia> u(0.1) + exp(0.1)*sum(u)\n1.3777980523127336Note that DefiniteIntegral is a functional, i.e., a 1 × ∞ operator.  when multiplied on the left by a function, it automatically constructs the operator rm e^x int_-1^1 f(x) dx viaDocTestSetup = quote\n    using ApproxFun\n    x = Fun()\n    Q = DefiniteIntegral(Chebyshev())\nendjulia> M = Multiplication(exp(x),ConstantSpace())\nConcreteMultiplication:ConstantSpace→Chebyshev(【-1.0,1.0】)\n 1.26607    \n 1.13032    \n 0.271495   \n 0.0443368  \n 0.00547424\n 0.000542926\n 4.49773e-5\n 3.19844e-6\n 1.99212e-7\n 1.10368e-8\n  ⋮         \n\njulia> M*Q\nTimesOperator:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)\n 2.53213     0.0  -0.844044     0.0  …  0.0  -0.0401926    0.0  ⋯\n 2.26064     0.0  -0.753545     0.0     0.0  -0.0358831    0.0  ⋱\n 0.542991    0.0  -0.180997     0.0     0.0  -0.0086189    0.0  ⋱\n 0.0886737   0.0  -0.0295579    0.0     0.0  -0.00140752   0.0  ⋱\n 0.0109485   0.0  -0.00364949   0.0     0.0  -0.000173785  0.0  ⋱\n 0.00108585  0.0  -0.000361951  0.0  …  0.0  -1.72358e-5   0.0  ⋱\n 8.99546e-5  0.0  -2.99849e-5   0.0     0.0  -1.42785e-6   0.0  ⋱\n 6.39687e-6  0.0  -2.13229e-6   0.0     0.0  -1.01538e-7   0.0  ⋱\n 3.98425e-7  0.0  -1.32808e-7   0.0     0.0  -6.32421e-9   0.0  ⋱\n 2.20735e-8  0.0  -7.35785e-9   0.0     0.0  -3.50374e-10  0.0  ⋱\n  ⋮           ⋱     ⋱            ⋱   …   ⋱     ⋱            ⋱   ⋱Note that Q*exp(x) applies the operator to a function.  To construct the operator that first multiplies by exp(x), use Q[exp(x)].  This is equivalent to Q*Multiplication(exp(x),Chebyshev())."
},

{
    "location": "usage/operators.html#Operators-and-space-promotion-1",
    "page": "Operators",
    "title": "Operators and space promotion",
    "category": "section",
    "text": "It is often more convenient to not specify a space explicitly, but rather infer it when the operator is used.  For example, we can construct Derivative(), which has the alias 𝒟, and represents the first derivative on any space:julia> f = Fun(cos,Chebyshev(0..1)); (𝒟*f)(0.1)\n-0.09983341664681705\n\njulia> f = Fun(cos,Fourier()); (𝒟*f)(0.1)\n-0.09983341664682804Behind the scenes, Derivative() is equivalent to Derivative(UnsetSpace(),1). When multiplying a function f, the domain space is promoted before multiplying, that is, Derivative()*f is equivalent to Derivative(space(f))*f.  This promotion of the domain space happens even when operators have spaces attached. This facilitates the following construction:julia> D = Derivative(Chebyshev());\n\njulia> D^2\nConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(2,【-1.0,1.0】)\n 0.0  0.0  4.0                                           \n      0.0  0.0  6.0                                      \n           0.0  0.0  8.0                                 \n                0.0  0.0  10.0                           \n                     0.0   0.0  12.0                     \n                           0.0   0.0  14.0               \n                                 0.0   0.0  16.0         \n                                       0.0   0.0  18.0   \n                                             0.0   0.0  ⋱\n                                                   0.0  ⋱\n                                                        ⋱Note that rangespace(D) ≠ Chebyshev(), hence the operators are not compatible. Therefore, it has thrown away its domain space, and thus this is equivalent to Derivative(rangespace(D))*D.DocTestSetup = nothing"
},

{
    "location": "usage/operators.html#Concatenating-operators-1",
    "page": "Operators",
    "title": "Concatenating operators",
    "category": "section",
    "text": "The concatenation functions vcat, hcat and hvcat are overriden for operators to represent the resulting combined operator, now with a rangespace or domainspace that is an ArraySpace. "
},

{
    "location": "usage/equations.html#",
    "page": "Linear Equations",
    "title": "Linear Equations",
    "category": "page",
    "text": ""
},

{
    "location": "usage/equations.html#Linear-equations-1",
    "page": "Linear Equations",
    "title": "Linear equations",
    "category": "section",
    "text": "Linear equations such as ordinary and partial differential equations,  fractional differential equations and integral equations can be solved using ApproxFun. This is accomplished using A\\b where A is an Operator and b is a Fun.  As a simple example, consider the equationu(theta) + cu(theta) = costhetawhere we want a solution that is periodic on 02pi).  This can be solved succinctly as follows:DocTestSetup = quote\n    using ApproxFun\nendjulia> b = Fun(cos,Fourier());\n\njulia> c = 0.1; u = (𝒟+c*I)\\b;\n\njulia> u(0.6)\n0.64076835137228\n\njulia> (c*cos(0.6)+sin(0.6))/(1+c^2)  # exact solution\n0.6407683513722804Recall that 𝒟 is an alias to Derivative() == Derivative(UnsetSpace(),1).As another example, consider the Fredholm integral equationu + rm e^x int_-1^1 cos x  u(x) rm dx = cos rm e^xWe can solve this equation as follows:julia> Σ = DefiniteIntegral(Chebyshev()); x=Fun();\n\njulia> u = (I+exp(x)*Σ[cos(x)])\\cos(exp(x));\n\njulia> u(0.1)\n0.21864294855628802Note that we used the syntax op[f::Fun], which is a shorthand for op*Multiplication(f)."
},

{
    "location": "usage/equations.html#Boundary-conditions-1",
    "page": "Linear Equations",
    "title": "Boundary conditions",
    "category": "section",
    "text": "Incorporating boundary conditions into differential equations is important so that the equation is well-posed.  This is accomplished via combining operators and functionals (i.e., 1 × ∞ operators).  As a simple example, consider the first order initial value problemu = t u qquadhboxandqquad u(0) = 1To pose this in ApproxFun, we want to find a u such that Evaluation(0)*u == 1 and (𝒟 - t)*u == 0.  This is accomplished via:julia> t = Fun(0..1);\n\njulia> u = [Evaluation(0); 𝒟 - t]  \\ [1;0];\n\njulia> u(0)\n0.9999999999999996\n\njulia> norm(u\'-t*u)\n1.2016080299388273e-16Behind the scenes, the Vector{Operator{T}} representing the functionals and operators are combined into a single InterlaceOperator.A common usage is two-point boundary value problems. Consider the singularly perturbed boundary value problem:epsilon u-xu+u = u qquad u(-1) = 1quad u(1) = 2This can be solved in ApproxFun via:julia> x = Fun();\n\njulia> u = [Evaluation(-1);\n            Evaluation(1);\n            1/70*𝒟^2-x*𝒟+I] \\ [1,2,0];\n\njulia> u(0.1)\n0.049999999999960326Note in this case the space is inferred from the variable coefficient x.This ODE can also be solved using the Dirichlet operator:julia> x = Fun();\n\njulia> u = [Dirichlet();\n            1/70*𝒟^2-x*𝒟+I] \\ [[1,2],0];\n\njulia> u(0.1)\n0.04999999999996019"
},

{
    "location": "usage/equations.html#Systems-of-equations-1",
    "page": "Linear Equations",
    "title": "Systems of equations",
    "category": "section",
    "text": "Systems of equations can be handled by creating a matrix of operators and functionals.  For example, we can solve the systembeginalign*\n    u - u + 2v = rm e^x  cr\n    v + v = cos(x) cr\n    u(-1) = u(-1) = v(-1) = 0\nendalign*using the following code:julia> x = Fun(); B = Evaluation(Chebyshev(),-1);\n\njulia> A = [B      0;\n            B*𝒟    0;\n            0      B;\n            𝒟^2-I  2I;\n            I      𝒟+I];\n\njulia> u,v = A\\[0;0;0;exp(x);cos(x)];\n\njulia> u(-1),u\'(-1),v(-1)\n(-4.163336342344337e-17,-2.7755575615628914e-16,-2.220446049250313e-16)\n\njulia> norm(u\'\'-u+2v-exp(x))\n5.981056979045254e-16\n\njulia> norm(u + v\'+v-cos(x))\n2.3189209621240424e-16In this example, the automatic space detection failed and so we needed to specify explicitly that the domain space for B is Chebyshev()."
},

{
    "location": "usage/equations.html#QR-Factorization-1",
    "page": "Linear Equations",
    "title": "QR Factorization",
    "category": "section",
    "text": "Behind the scenes, A\\b where A is an Operator is implemented via an adaptive QR factorization.  That is, it is equivalent to qrfact(A)\\b.  (There is a subtly here in space inferring: A\\b can use     both A and b to determine the domain space, while qrfact(A) only     sees the operator A.)       Note that qrfact adaptively caches a partial QR Factorization as it is applied to different right-hand sides, so the same operator can be inverted much more efficiently in subsequent problems."
},

{
    "location": "usage/equations.html#Partial-differential-equations-1",
    "page": "Linear Equations",
    "title": "Partial differential equations",
    "category": "section",
    "text": "Partial differential operators are also supported.  Here\'s an example of solving the Poisson equation with zero boundary conditions:d = (-1..1)^2\nx,y = Fun(d)\nf = exp.(-10(x+0.3)^2-20(y-0.2)^2)  # use broadcasting as exp(f) not implemented in 2D\nA = [Dirichlet(d);Δ]              # Δ is an alias for Laplacian()\n@time u = A \\ [zeros(∂(d));f]     #4s for ~3k coefficientsUsing a QR Factorization reduces the cost of subsequent calls substantially:QR = qrfact(A)\n@time QR \\ [zeros(∂(d));f]   # 4s\ng = exp.(-10(x+0.2)^2-20(y-0.1)^2)\n@time QR \\ [zeros(∂(d));g]  # 0.09sMany PDEs have weak singularities at the corners, in which case it is beneficial to specify a tolerance to reduce the time:\\(A,[zeros(∂(d));f]; tolerance=1E-6)"
},

{
    "location": "usage/equations.html#Nonlinear-equations-1",
    "page": "Linear Equations",
    "title": "Nonlinear equations",
    "category": "section",
    "text": "There is preliminary support for nonlinear equations, via Newton iteration in function space.  Here is a simple two-point boundary value problem:beginalign*\n    epsilon u + 6(1-x^2)u +u^2=1 cr\n    u(-1)=u(1)=0\nendalign*This can be solved usingx = Fun()\nN = u->[u(-1.)-c;u(1.);ε*u\'\'+6*(1-x^2)*u\'+u^2-1.]\nu = newton(N,u0)"
},

{
    "location": "faq.html#",
    "page": "Frequently Asked Questions",
    "title": "Frequently Asked Questions",
    "category": "page",
    "text": ""
},

{
    "location": "faq.html#Frequently-Asked-Questions-1",
    "page": "Frequently Asked Questions",
    "title": "Frequently Asked Questions",
    "category": "section",
    "text": ""
},

{
    "location": "faq.html#Approximating-functions-1",
    "page": "Frequently Asked Questions",
    "title": "Approximating functions",
    "category": "section",
    "text": ""
},

{
    "location": "faq.html#How-do-I-interpolate-a-function-at-a-specified-grid?-1",
    "page": "Frequently Asked Questions",
    "title": "How do I interpolate a function at a specified grid?",
    "category": "section",
    "text": "In the case where the grid is specified by points(space,n), you can apply the default transform to data:DocTestSetup = quote\n    using ApproxFun\nendjulia> S = Chebyshev(1..2);\n\njulia> p = points(S,20); # the default grid\n\njulia> v = exp.(p);      # values at the default grid\n\njulia> f = Fun(S,ApproxFun.transform(S,v));\n\njulia> f(1.1)\n3.0041660239464347\n\njulia> exp(1.1)\n3.0041660239464334ApproxFun has no inbuilt support for interpolating functions at other sets of points, but this can be accomplished manually by evaluating the basis at the set of points and using \\:julia> S = Chebyshev(1..2);\n\njulia> n = 50;\n\njulia> p = linspace(1,2,n);   # a non-default grid\n\njulia> v = exp.(p);           # values at the non-default grid\n\njulia> V = Array(Float64,n,n); # Create a Vandermonde matrix by evaluating the basis at the grid\n\njulia> for k = 1:n\n           V[:,k] = Fun(S,[zeros(k-1);1]).(p)\n       end\n\njulia> f = Fun(S,V\\v);\n\njulia> f(1.1)\n3.0041660228311926\n\njulia> exp(1.1)\n3.0041660239464334Note that an evenly spaced grid suffers from instability for large n.  The easiest way around this is to use least squares with more points than coefficients, instead of interpolation:julia> S = Chebyshev(1..2);\n\njulia> n = 100; m = 50;\n\njulia> p = linspace(1,2,n);   # a non-default grid\n\njulia> v = exp.(p);           # values at the non-default grid\n\njulia> V = Array(Float64,n,m); # Create a Vandermonde matrix by evaluating the basis at the grid\n\njulia> for k = 1:m\n           V[:,k] = Fun(S,[zeros(k-1);1]).(p)\n       end\n\njulia> f = Fun(S,V\\v);\n\njulia> f(1.1)\n3.004166023946434\n\njulia> exp(1.1)\n3.0041660239464334We can use this same approach for multivariate functions:julia> S = Chebyshev(0..1)^2;\n\njulia> n = 1000; m = 50;\n\njulia> srand(0); x = rand(n); y = rand(n);\n\njulia> v = exp.(x .* cos(y));  # values at the non-default grid\n\njulia> V = Array(Float64,n,m); # Create a Vandermonde matrix by evaluating the basis at the grid\n\njulia> for k = 1:m\n          V[:,k] = Fun(S,[zeros(k-1);1]).(x,y)\n       end\n\n\njulia> f = Fun(S,V\\v);\n\njulia> f(0.1,0.2)\n1.1029700685084018\n\njulia> exp(0.1*cos(0.2))\n1.1029701284210731DocTestSetup = nothing"
},

{
    "location": "library.html#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#Library-1",
    "page": "Library",
    "title": "Library",
    "category": "section",
    "text": ""
},

{
    "location": "library.html#ApproxFun.Fun",
    "page": "Library",
    "title": "ApproxFun.Fun",
    "category": "Type",
    "text": "Fun(s::Space,coefficients::AbstractVector)\n\nreturns a Fun with the specified coefficients in the space s\n\n\n\nFun(f,s::Space)\n\nreturn a Fun representing the function, number, or vector f in the space s.  If f is vector-valued, it returns a vector-valued analogue of s.\n\n\n\nFun(f,d::Domain)\n\nreturns Fun(f,Space(d)), that is, it uses the default space for the specified domain.\n\n\n\nFun(s::Space)\n\nreturns Fun(identity,s)\n\n\n\nFun(f)\n\nreturns Fun(f,Chebyshev())\n\n\n\nFun()\n\nreturns Fun(identity,Chebyshev()).\n\n\n\n"
},

{
    "location": "library.html#Base.ones-Tuple{ApproxFun.Space}",
    "page": "Library",
    "title": "Base.ones",
    "category": "Method",
    "text": "ones(d::Space)\n\nReturn the Fun that represents the function one on the specified space.\n\n\n\n"
},

{
    "location": "library.html#Base.zeros-Tuple{ApproxFun.Space}",
    "page": "Library",
    "title": "Base.zeros",
    "category": "Method",
    "text": "zeros(d::Space)\n\nReturn the Fun that represents the function one on the specified space.\n\n\n\n"
},

{
    "location": "library.html#Constructing-a-Fun-1",
    "page": "Library",
    "title": "Constructing a Fun",
    "category": "section",
    "text": "Funones(::Space)zeros(::Space)"
},

{
    "location": "library.html#ApproxFun.Arc",
    "page": "Library",
    "title": "ApproxFun.Arc",
    "category": "Type",
    "text": "Arc(c,r,(θ₁,θ₂))\n\nrepresents the arc centred at c with radius r from angle θ₁ to θ₂.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Circle",
    "page": "Library",
    "title": "ApproxFun.Circle",
    "category": "Type",
    "text": "Circle(c,r,o)\n\nrepresents the circle centred at c with radius r which is positively (o=true) or negatively (o=false) oriented.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Curve",
    "page": "Library",
    "title": "ApproxFun.Curve",
    "category": "Constant",
    "text": "Curve Represents a domain defined by the image of a Fun.  Example usage would be\n\nx=Fun(1..2)\nCurve(exp(im*x))  # represents an arc\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Disk",
    "page": "Library",
    "title": "ApproxFun.Disk",
    "category": "Type",
    "text": "Disk(c,r)\n\nrepresents the disk centred at c with radius r.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Segment",
    "page": "Library",
    "title": "ApproxFun.Segment",
    "category": "Type",
    "text": "Segment(a,b)\n\nrepresents a line segment from a to b.  In the case where a and b are real and a < b, then this is is equivalent to an Interval(a,b).\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Interval",
    "page": "Library",
    "title": "ApproxFun.Interval",
    "category": "Function",
    "text": "Interval(a::Real,b::Real)\n\nrepresents the set {x : a ≤ x ≤ b}.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Line",
    "page": "Library",
    "title": "ApproxFun.Line",
    "category": "Type",
    "text": "Line{a}(c)\n\nrepresents the line at angle a in the complex plane, centred at c.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.PeriodicInterval",
    "page": "Library",
    "title": "ApproxFun.PeriodicInterval",
    "category": "Type",
    "text": "PeriodicInterval(a,b)\n\nrepresents a periodic interval from a to b, that is, the point b is identified with a.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Point",
    "page": "Library",
    "title": "ApproxFun.Point",
    "category": "Type",
    "text": "Point(x)\n\nrepresents a single point at x.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ProductDomain",
    "page": "Library",
    "title": "ApproxFun.ProductDomain",
    "category": "Type",
    "text": "ProductDomain((d1,d2))\n\nrepresents the product of two domains, the set {(x,y) : x ∈ d1 & y ∈ d2}.\n\nMultiplication of domains is overrident to return a ProductDomain. For example, the following represents the rectangle 1 ≤ x ≤ 2 & 3 ≤ y ≤ 4:\n\nInterval(1,2)*(3,4)\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Ray",
    "page": "Library",
    "title": "ApproxFun.Ray",
    "category": "Type",
    "text": "Ray{a}(c,o)\n\nrepresents a ray at angle a starting at c, with orientation out to infinity (o = true) or back from infinity (o = false).\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.UnionDomain",
    "page": "Library",
    "title": "ApproxFun.UnionDomain",
    "category": "Type",
    "text": "UnionDomain((d1,d2,…,dn))\n\nrepresents a union of multiple subdomains: {x : x ∈ d1 || … || x ∈ dn}.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.∂",
    "page": "Library",
    "title": "ApproxFun.∂",
    "category": "Function",
    "text": "∂(d::Domain)\n\nreturns the boundary of d.  For example, the boundary of a Disk() is a Circle(), and the boundary of Interval()^2 is a piecewise interval that sketches the boundary of a rectangle.\n\n\n\n"
},

{
    "location": "library.html#Domains-1",
    "page": "Library",
    "title": "Domains",
    "category": "section",
    "text": "ArcCircleCurveDiskSegmentIntervalLinePeriodicIntervalApproxFun.PointProductDomainRayUnionDomain∂"
},

{
    "location": "library.html#ApproxFun.canonicalspace",
    "page": "Library",
    "title": "ApproxFun.canonicalspace",
    "category": "Function",
    "text": "canonicalspace(s::Space)\n\nreturns a space that is used as a default to implement missing functionality, e.g., evaluation.  Implement a Conversion operator or override coefficients to support this.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.itransform",
    "page": "Library",
    "title": "ApproxFun.itransform",
    "category": "Function",
    "text": "itransform(s::Space,coefficients::AbstractVector)\n\nTransform coefficients back to values.  Defaults to using canonicalspace as in transform.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.transform",
    "page": "Library",
    "title": "ApproxFun.transform",
    "category": "Function",
    "text": "transform(s::Space,vals::Vector)\n\nTransform values on the grid specified by points(s,length(vals)) to coefficients in the space s. Defaults to coefficients(transform(canonicalspace(space),values),canonicalspace(space),space)\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.evaluate-Tuple{ApproxFun.Space,AbstractArray{T,1} where T,Any}",
    "page": "Library",
    "title": "ApproxFun.evaluate",
    "category": "Method",
    "text": "evaluate(sp::Space,coefficients::AbstractVector,x)\n\nEvaluates the expansion at a point x. If x is in the domain, then this should return zero.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.dimension-Tuple{ApproxFun.Space}",
    "page": "Library",
    "title": "ApproxFun.dimension",
    "category": "Method",
    "text": "dimension(s::Space)\n\nreturns the dimension of s, which is the maximum number of coefficients.\n\n\n\n"
},

{
    "location": "library.html#Accessing-information-about-a-spaces-1",
    "page": "Library",
    "title": "Accessing information about a spaces",
    "category": "section",
    "text": "ApproxFun.canonicalspaceitransformtransformevaluate(::Space,::AbstractVector,::)ApproxFun.dimension(::Space)"
},

{
    "location": "library.html#ApproxFun.SequenceSpace",
    "page": "Library",
    "title": "ApproxFun.SequenceSpace",
    "category": "Type",
    "text": "SequenceSpace is the space of all sequences, i.e., infinite vectors. Also denoted ℓ⁰.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ConstantSpace",
    "page": "Library",
    "title": "ApproxFun.ConstantSpace",
    "category": "Type",
    "text": "ConstantSpace is the 1-dimensional scalar space.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Chebyshev",
    "page": "Library",
    "title": "ApproxFun.Chebyshev",
    "category": "Type",
    "text": "Chebyshev() is the space spanned by the Chebyshev polynomials\n\n    T_0(x),T_1(x),T_2(x),…\n\nwhere T_k(x) = cos(k*acos(x)).  This is the default space as there exists a fast transform and general smooth functions on [-1,1] can be easily resolved.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Jacobi",
    "page": "Library",
    "title": "ApproxFun.Jacobi",
    "category": "Type",
    "text": "Jacobi(b,a) represents the space spanned by Jacobi polynomials P_k^{(a,b)}, which are orthogonal with respect to the weight (1+x)^β*(1-x)^α\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Ultraspherical",
    "page": "Library",
    "title": "ApproxFun.Ultraspherical",
    "category": "Type",
    "text": "Ultraspherical(λ) is the space spanned by the ultraspherical polynomials\n\n    C_0^{(λ)}(x),C_1^{(λ)}(x),C_2^{(λ)}(x),…\n\nNote that λ=1 this reduces to Chebyshev polynomials of the second kind: C_k^{(1)}(x) = U_k(x).\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Taylor",
    "page": "Library",
    "title": "ApproxFun.Taylor",
    "category": "Type",
    "text": "Taylor() is the space spanned by [1,z,z^2,...]. This is a type alias for Hardy{true}.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Hardy",
    "page": "Library",
    "title": "ApproxFun.Hardy",
    "category": "Type",
    "text": "Hardy{false}() is the space spanned by [1/z,1/z^2,...]. Hardy{true}() is the space spanned by [1,z,z^2,...].\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Fourier",
    "page": "Library",
    "title": "ApproxFun.Fourier",
    "category": "Type",
    "text": "Fourier() is the space spanned by the trigonemtric polynomials\n\n    1,sin(θ),cos(θ),sin(2θ),cos(2θ),…\n\nSee also Laurent.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Laurent",
    "page": "Library",
    "title": "ApproxFun.Laurent",
    "category": "Type",
    "text": "Laurent() is the space spanned by the complex exponentials\n\n    1,exp(-im*θ),exp(im*θ),exp(-2im*θ),…\n\nSee also Fourier.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.CosSpace",
    "page": "Library",
    "title": "ApproxFun.CosSpace",
    "category": "Type",
    "text": "CosSpace() is the space spanned by [1,cos θ,cos 2θ,...]\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.SinSpace",
    "page": "Library",
    "title": "ApproxFun.SinSpace",
    "category": "Type",
    "text": "SinSpace() is the space spanned by [sin θ,sin 2θ,...]\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.JacobiWeight",
    "page": "Library",
    "title": "ApproxFun.JacobiWeight",
    "category": "Type",
    "text": "JacobiWeight(β,α,s::Space)\n\nweights a space s by a Jacobi weight, which on -1..1 is (1+x)^β*(1-x)^α. For other domains, the weight is inferred by mapping to -1..1.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.LogWeight",
    "page": "Library",
    "title": "ApproxFun.LogWeight",
    "category": "Type",
    "text": "LogWeight(β,α,s::Space)\n\nrepresents a function on -1..1 weighted by log((1+x)^β*(1-x)^α). For other domains, the weight is inferred by mapping to -1..1.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ArraySpace",
    "page": "Library",
    "title": "ApproxFun.ArraySpace",
    "category": "Type",
    "text": "ArraySpace(s::Space,dims...)\n\nis used to represent array-valued expansions in a space s.  The coefficients are of each entry are interlaced.\n\nFor example,\n\nf = Fun(x->[exp(x),sin(x)],-1..1)\nspace(f) == ArraySpace(Chebyshev(),2)\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.TensorSpace",
    "page": "Library",
    "title": "ApproxFun.TensorSpace",
    "category": "Type",
    "text": "TensorSpace(a::Space,b::Space)\n\nrepresents a tensor product of two 1D spaces a and b. The coefficients are interlaced in lexigraphical order.\n\nFor example, consider\n\nFourier()*Chebyshev()  # returns TensorSpace(Fourier(),Chebyshev())\n\nThis represents functions on [-π,π) x [-1,1], using the Fourier basis for the first argument and Chebyshev basis for the second argument, that is, φ_k(x)T_j(y), where\n\nφ_0(x) = 1,\nφ_1(x) = sin x,\nφ_2(x) = cos x,\nφ_3(x) = sin 2x,\nφ_4(x) = cos 2x\n…\n\nBy Choosing (k,j) appropriately, we obtain a single basis:\n\nφ_0(x)T_0(y) (= 1),\nφ_0(x)T_1(y) (= y),\nφ_1(x)T_0(y) (= sin x),\nφ_0(x)T_2(y), …\n\n\n\n"
},

{
    "location": "library.html#Inbuilt-spaces-1",
    "page": "Library",
    "title": "Inbuilt spaces",
    "category": "section",
    "text": "SequenceSpaceConstantSpaceChebyshevJacobiUltrasphericalTaylorHardyFourierLaurentCosSpaceSinSpaceJacobiWeightApproxFun.LogWeightApproxFun.ArraySpaceTensorSpace"
},

{
    "location": "library.html#ApproxFun.domain",
    "page": "Library",
    "title": "ApproxFun.domain",
    "category": "Function",
    "text": "domain(f::Fun)\n\nreturns the domain that f is defined on\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.coefficients",
    "page": "Library",
    "title": "ApproxFun.coefficients",
    "category": "Function",
    "text": "coefficients(f::Fun) -> Vector\n\nreturns the coefficients of f, corresponding to the space space(f).\n\n\n\ncoefficients(f::Fun,s::Space) -> Vector\n\nreturns the coefficients of f in the space s, which may not be the same as space(f).\n\n\n\ncoefficients(cfs::AbstractVector,fromspace::Space,tospace::Space) -> Vector\n\nconverts coefficients in fromspace to coefficients in tospace\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.extrapolate",
    "page": "Library",
    "title": "ApproxFun.extrapolate",
    "category": "Function",
    "text": "extrapolate(f::Fun,x)\n\nreturns an extrapolation of f from its domain to x.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ncoefficients",
    "page": "Library",
    "title": "ApproxFun.ncoefficients",
    "category": "Function",
    "text": "ncoefficients(f::Fun) -> Integer\n\nreturns the number of coefficients of a fun\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.points",
    "page": "Library",
    "title": "ApproxFun.points",
    "category": "Function",
    "text": "points(f::Fun)\n\nreturns a grid of points that f can be transformed into values and back\n\n\n\npoints(s::Space,n::Integer)\n\nreturns a grid of approximately n points, for which a transform exists from values at the grid to coefficients in the space s.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.space",
    "page": "Library",
    "title": "ApproxFun.space",
    "category": "Function",
    "text": "space(f::Fun)\n\nreturns the space of f\n\n\n\n"
},

{
    "location": "library.html#Base.values",
    "page": "Library",
    "title": "Base.values",
    "category": "Function",
    "text": "values(f::Fun)\n\nreturns f evaluated at points(f)\n\n\n\n"
},

{
    "location": "library.html#Base.stride-Tuple{ApproxFun.Fun}",
    "page": "Library",
    "title": "Base.stride",
    "category": "Method",
    "text": "stride(f::Fun)\n\nreturns the stride of the coefficients, checked numerically\n\n\n\n"
},

{
    "location": "library.html#Accessing-information-about-a-Fun-1",
    "page": "Library",
    "title": "Accessing information about a Fun",
    "category": "section",
    "text": "domaincoefficientsextrapolatencoefficientspointsspaceApproxFun.valuesBase.stride(::Fun)"
},

{
    "location": "library.html#ApproxFun.reverseorientation",
    "page": "Library",
    "title": "ApproxFun.reverseorientation",
    "category": "Function",
    "text": "reverseorientation(f::Fun)\n\nreturn f on a reversed orientated contour.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.setdomain",
    "page": "Library",
    "title": "ApproxFun.setdomain",
    "category": "Function",
    "text": "setdomain(f::Fun,d::Domain)\n\nreturns f projected onto domain\n\n\n\n"
},

{
    "location": "library.html#Base.chop",
    "page": "Library",
    "title": "Base.chop",
    "category": "Function",
    "text": "chop(f::Fun,tol) -> Fun\n\nreduces the number of coefficients by dropping the tail that is below the specified tolerance.\n\n\n\n"
},

{
    "location": "library.html#Modify-a-Fun-1",
    "page": "Library",
    "title": "Modify a Fun",
    "category": "section",
    "text": "reverseorientationApproxFun.setdomainchop"
},

{
    "location": "library.html#ApproxFun.Operator",
    "page": "Library",
    "title": "ApproxFun.Operator",
    "category": "Type",
    "text": "Operator{T}\n\nis an abstract type to represent linear operators between spaces.\n\n\n\n"
},

{
    "location": "library.html#BandedMatrices.bandwidths-Tuple{ApproxFun.Operator}",
    "page": "Library",
    "title": "BandedMatrices.bandwidths",
    "category": "Method",
    "text": "bandwidths(op::Operator)\n\nreturns the bandwidth of op in the form (l,u), where l ≥ 0 represents the number of subdiagonals and u ≥ 0 represents the number of superdiagonals.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.domainspace",
    "page": "Library",
    "title": "ApproxFun.domainspace",
    "category": "Function",
    "text": "domainspace(op::Operator)\n\ngives the domain space of op.  That is, op*f will first convert f to a Fun in the space domainspace(op) before applying the operator.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.rangespace",
    "page": "Library",
    "title": "ApproxFun.rangespace",
    "category": "Function",
    "text": "rangespace(op::Operator)\n\ngives the range space of op.  That is, op*f will return a Fun in the space rangespace(op), provided f can be converted to a Fun in domainspace(op).\n\n\n\n"
},

{
    "location": "library.html#Base.getindex-Tuple{ApproxFun.Operator,Any,Any}",
    "page": "Library",
    "title": "Base.getindex",
    "category": "Method",
    "text": "op[k,j]\n\nreturns the kth coefficient of op*Fun([zeros(j-1);1],domainspace(op)).\n\n\n\n"
},

{
    "location": "library.html#Base.getindex-Tuple{ApproxFun.Operator,ApproxFun.Fun}",
    "page": "Library",
    "title": "Base.getindex",
    "category": "Method",
    "text": "op[f::Fun]\n\nconstructs the operator op*Multiplication(f), that is, it multiplies on the right by f first.  Note that op*f is different: it applies op to f.\n\n\n\n"
},

{
    "location": "library.html#Base.:\\-Tuple{ApproxFun.Operator,Any}",
    "page": "Library",
    "title": "Base.:\\",
    "category": "Method",
    "text": "\\(A::Operator,b;tolerance=tol,maxlength=n)\n\nsolves a linear equation, usually differential equation, where A is an operator or array of operators and b is a Fun or array of funs.  The result u will approximately satisfy A*u = b.\n\n\n\n"
},

{
    "location": "library.html#Base.LinAlg.qrfact-Tuple{ApproxFun.Operator}",
    "page": "Library",
    "title": "Base.LinAlg.qrfact",
    "category": "Method",
    "text": "qrfact(A::Operator)\n\nreturns a cached QR factorization of the Operator A.  The result QR enables solving of linear equations: if u=QR\\b, then u approximately satisfies A*u = b.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.cache-Tuple{ApproxFun.Operator}",
    "page": "Library",
    "title": "ApproxFun.cache",
    "category": "Method",
    "text": "cache(op::Operator)\n\nCaches the entries of an operator, to speed up multiplying a Fun by the operator.\n\n\n\n"
},

{
    "location": "library.html#Operators-1",
    "page": "Library",
    "title": "Operators",
    "category": "section",
    "text": "OperatorBandedMatrices.bandwidths(::Operator)domainspacerangespaceBase.getindex(::Operator,::,::)Base.getindex(::Operator,::Fun)\\(::Operator,::)qrfact(::Operator)cache(::Operator)"
},

{
    "location": "library.html#ApproxFun.Conversion",
    "page": "Library",
    "title": "ApproxFun.Conversion",
    "category": "Type",
    "text": "Conversion(fromspace::Space,tospace::Space)\n\nrepresents a conversion operator between fromspace and tospace, when available.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Derivative",
    "page": "Library",
    "title": "ApproxFun.Derivative",
    "category": "Type",
    "text": "Derivative(sp::Space,k::Int) represents the k-th derivative on sp.\n\n\n\nDerivative(sp::Space,k::Vector{Int}) represents a partial derivative on a multivariate space. For example,\n\nDx = Derivative(Chebyshev()^2,[1,0]) # ∂/∂x\nDy = Derivative(Chebyshev()^2,[0,1]) # ∂/∂y\n\n\n\nDerivative(sp::Space) represents the first derivative Derivative(sp,1).\n\n\n\nDerivative(k) represents the k-th derivative, acting on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\nDerivative() represents the first derivative on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Dirichlet",
    "page": "Library",
    "title": "ApproxFun.Dirichlet",
    "category": "Type",
    "text": "Dirichlet(sp,k) is the operator associated with restricting the k-th derivative on the boundary for the space sp.\n\n\n\nDirichlet(sp) is the operator associated with restricting the  the boundary for the space sp.\n\n\n\nDirichlet() is the operator associated with restricting on the  the boundary.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Evaluation",
    "page": "Library",
    "title": "ApproxFun.Evaluation",
    "category": "Type",
    "text": "Evaluation(sp,x,k) is the functional associated with evaluating the k-th derivative at a point x for the space sp.\n\n\n\nEvaluation(sp,x) is the functional associated with evaluating at a point x for the space sp.\n\n\n\nEvaluation(x) is the functional associated with evaluating at a point x.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Integral",
    "page": "Library",
    "title": "ApproxFun.Integral",
    "category": "Type",
    "text": "Integral(sp::Space,k::Int) represents a k-th integral on sp. There is no guarantee on normalization.\n\n\n\nIntegral(sp::Space) represents the first integral Integral(sp,1).\n\n\n\nIntegral(k)represents thek`-th integral, acting on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\nIntergral() represents the first integral on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Laplacian",
    "page": "Library",
    "title": "ApproxFun.Laplacian",
    "category": "Type",
    "text": "Laplacian(sp::Space) represents the laplacian on space sp.\n\n\n\nLaplacian() represents the laplacian on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Multiplication",
    "page": "Library",
    "title": "ApproxFun.Multiplication",
    "category": "Type",
    "text": "Multiplication(f::Fun,sp::Space) is the operator representing multiplication by f on functions in the space sp.\n\n\n\nMultiplication(f::Fun) is the operator representing multiplication by f on an unset space of functions.  Spaces will be inferred when applying or manipulating the operator.\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Neumann",
    "page": "Library",
    "title": "ApproxFun.Neumann",
    "category": "Function",
    "text": "Neumann(sp) is the operator associated with restricting the normal derivative on the boundary for the space sp. At the moment it is implemented as Dirichlet(sp,1).\n\n\n\n`Neumann( is the operator associated with restricting the normal derivative on the boundary.\n\n\n\n"
},

{
    "location": "library.html#Inbuilt-operators-1",
    "page": "Library",
    "title": "Inbuilt operators",
    "category": "section",
    "text": "ConversionDerivativeDirichletEvaluationIntegralLaplacianMultiplicationNeumann"
},

{
    "location": "library.html#ApproxFun.LowRankFun",
    "page": "Library",
    "title": "ApproxFun.LowRankFun",
    "category": "Type",
    "text": "LowRankFun gives an approximation to a bivariate function in low rank form.\n\n\n\n"
},

{
    "location": "library.html#Bivariate-1",
    "page": "Library",
    "title": "Bivariate",
    "category": "section",
    "text": "LowRankFun"
},

]}
