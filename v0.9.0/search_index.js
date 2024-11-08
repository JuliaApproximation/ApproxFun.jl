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
    "text": "A basic approach of computational mathematics that ApproxFun exploits is expansion in a basis $ f(x) \\approx \\sum{k=1}^n ck \\psik(x) $ Some traditional examples of bases \\psi1(x),\\psi_2(x),\\ldots$ areTaylor series: 1zz^2ldots\nFourier series (for periodic functions on 0..2π): 1sin x cos x sin 2 x ldots\nChebyshev series (for non-periodic functions on -1..1): 1xcos 2 hboxacos x cos 3 hboxacos x ldotsIn ApproxFun, functions are represented by a Fun with two components: space, which dictates the basis and coefficients which is a finite vector of coefficients.  Note that each Fun can have a different length vector of coefficients, allowing for approximation of many different functions to high accuracy.  The approximation by a Fun can be determined by a variety of methods:(1) Explicitly specifying the coefficients:julia> f = Fun(Taylor(),[1,2,3]) # Represents 1 + 2z + 3z^2\nFun(Taylor(🕒),[1.0,2.0,3.0])\n\njulia> f(1.0)\n6.0(2) Constructors take in a Function and adaptively determine the     number of coefficients.  For example,julia> Fun(exp)\nFun(Chebyshev(【-1.0,1.0】),[1.26607,1.13032,0.271495,0.0443368,0.00547424,0.000542926,4.49773e-5,3.19844e-6,1.99212e-7,1.10368e-8,5.5059e-10,2.49797e-11,1.03911e-12,3.99195e-14])determines that f can be approximated to roughly machine precision using 14 coefficients.  See Constructors for more information.(3) Manipulation of Funs give new Funs, where the number of coefficients is determined from the input.  The simplest example is addition, which for compatible bases is just padding the vectors to the same length and adding.  julia> a = Fun(cos,Chebyshev()); ncoefficients(a)\n13\n\njulia> b = Fun(x->cos(10cos(x^2)),Chebyshev()); ncoefficients(b)\n51\n\njulia> ncoefficients(a+b)\n51On the other hand, multiplication results in an approximation with more coefficients than either a or b, so that the result approximates the true a*b to roughly machine accuracy:julia> ncoefficients(a*b)\n63\n\njulia> a(0.1)*b(0.1) - (a*b)(0.1)\n1.1102230246251565e-16The example of multiplication highlights the importance of adaptivity: if with a fixed discretization size, operations like multiplication would lose accuracy when the true function is no longer resolved by the discretization.  More complicated examples are solving differential equations, where the coefficients of the solution can be determined adaptively, see Equations.ApproxFun supports a number of different spaces, as described in Spaces.  A key component of ApproxFun is support for interaction between different spaces.  This is crucial for efficient solution of differential equations, where linear operators are described as acting between different spaces, see Operators.  "
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
    "location": "usage/spaces.html#Classical-orthogonal-polynomial-spaces-1",
    "page": "Spaces",
    "title": "Classical orthogonal polynomial spaces",
    "category": "section",
    "text": "Chebyshev, Ultraspherical, Jacobi, Hermite, and Laguerre are spaces corresponding to expansion in classical orthogonal polynomials.Note that we always use the classical normalization: the basis are not orthonormal. This is because this normalization leads to rational recurrence relationships, which are more efficient than their normalized counterparts. See the Digital Library of Mathematical Functions for more information."
},

{
    "location": "usage/spaces.html#Chebyshev-space-1",
    "page": "Spaces",
    "title": "Chebyshev space",
    "category": "section",
    "text": "The default space in ApproxFun is Chebyshev, which represents expansions in Chebyshev polynomials:f(x) = sum_k=0^infty f_k T_k(x)where T_k(x) = cos k rm acos x, which are orthogonal polynomials with respect to the weight $ {1 \\over \\sqrt{1-x^2}} \\qquad\\hbox{for}\\qquad -1 \\leq x \\leq 1. $ Note that there is an intrinsic link between Chebyshev and CosSpace:  g(theta) = f(cos theta) = sum_k=0^infty f_k cos k thetaIn other words:DocTestSetup = quote\n    using ApproxFun\nendjulia> f=Fun(exp,Chebyshev());\n\njulia> g=Fun(CosSpace(),f.coefficients); # specify the coefficients directly\n\njulia> f(cos(0.1))\n2.70473560723178\n\njulia> g(0.1)\n2.7047356072317794"
},

{
    "location": "usage/spaces.html#Ultraspherical-spaces-1",
    "page": "Spaces",
    "title": "Ultraspherical spaces",
    "category": "section",
    "text": "A key tool for solving differential equations are the ultraspherical spaces, encoded as Ultraspherical(λ) for λ ≠ 0, which can be defined by the span of derivatives of Chebyshev polynomials, or alternatively as polynomials orthogonal with respect to the weight $ (1-x^2)^{\\lambda - 1/2} \\qquad\\hbox{for}\\qquad -1 \\leq x \\leq 1. $Note that Ultraspherical(1) corresponds to the Chebyshev basis of the second kind: U_k(x) = sin (k+1) rm acos x over sin rm acos x.   The relationship with Chebyshev polynomials follows from trigonemetric identities: T_k(x) = k U_k-1(x).  Converting between ultraspherical polynomials (with integer orders) is extremely efficient: it requires O(n) operations, where n is the number of coefficients."
},

{
    "location": "usage/spaces.html#Jacobi-spaces-1",
    "page": "Spaces",
    "title": "Jacobi spaces",
    "category": "section",
    "text": "Jacobi(b,a) represents the space spanned by the Jacobi polynomials, which are orthogonal polynomials with respect to the weight $ (1+x)^b(1-x)^a $ using the standard normalization."
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
    "text": "Some spaces are built out of other spaces:"
},

{
    "location": "usage/spaces.html#JacobiWeight-1",
    "page": "Spaces",
    "title": "JacobiWeight",
    "category": "section",
    "text": "JacobiWeight(β,α,space)  weights space, which is typically Chebyshev() or Jacobi(b,a),  by a Jacobi weight (1+x)^α*(1-x)^β: in other words, if the basis for space is psi_k(x) and the domain is the unit interval -1 .. 1, then the basis for JacobiWeight(β,α,space) is (1+x)^α(1-x)^β psi_k(x). If the domain is not the unit interval, then the basis is determined by mapping back to the unit interval: that is, if M(x) is the map dictated by tocanonical(space, x), where the canonical domain is the unit interval, then the basis is (1+M(x))^α(1-M(x))^β psi_k(x). For example, if the domain is another interval a .. b, then $ M(x) = {2x-b-a \\over b-a} $ and the basis becomes $ \\left({2 \\over (b-a)}\\right)^{\\alpha+\\eta}  (x-a)^α(b-x)^β \\psi_k(x) $"
},

{
    "location": "usage/spaces.html#SumSpace-1",
    "page": "Spaces",
    "title": "SumSpace",
    "category": "section",
    "text": "SumSpace((space_1,space_2,…,space_n)) represents the direct sum of the spaces, where evaluation is defined by adding up each component. A simple example is the following, showing that the coefficients are stored by interlacing:julia> x = Fun(identity,-1..1);\n\njulia> f = cos(x-0.1)*sqrt(1-x^2) + exp(x);\n\njulia> space(f) # isa SumSpace\n(1-x^2)^0.5[Chebyshev(【-1.0,1.0】)]⊕Chebyshev(【-1.0,1.0】)\n\njulia> a, b = components(f);\n\njulia> a(0.2) # cos(0.2-0.1)*sqrt(1-0.2^2)\n0.9749009987500246\n\njulia> b(0.2) # exp(0.2)\n1.2214027581601699\n\njulia> f(0.2) # a(0.2) + b(0.2)\n2.1963037569101944\n\njulia> norm(f.coefficients[1:2:end] - a.coefficients)\n0.0\n\njulia> norm(f.coefficients[2:2:end] - b.coefficients)\n0.0More complicated examples may interlace the coefficients using a different strategy. Note that it is difficult to represent the first component of function f by a Chebyshev series because the derivatives of f at its boundaries blow up, whereas the derivative of a polynomial is a polynomial.Note that Fourier and Laurent are currently implemented as SumSpace, but this may change in the future."
},

{
    "location": "usage/spaces.html#PiecewiseSpace-1",
    "page": "Spaces",
    "title": "PiecewiseSpace",
    "category": "section",
    "text": "PiecewiseSpace((space_1,space_2,…,space_n)) represents the direct sum of the spaces, where evaluation is defined in a piecewise way. A simple example is the following:julia> x = Fun(Domain(-1 .. 0) ∪ Domain( 1 .. 2));\n\njulia> f = exp(x);\n\njulia> a, b = components(f);\n\njulia> f(-0.5) - a(-0.5)\n0.0\n\njulia> f(1.5) - b(1.5)\n0.0\n\njulia> f(0.5) # outside domain components\n0.0\n\njulia> norm(f.coefficients[2:2:end] - b.coefficients)\n0.0\n\njulia> norm(f.coefficients[1:2:end] - a.coefficients)\n0.0More complicated examples may interlace the coefficients using a different strategy."
},

{
    "location": "usage/spaces.html#ArraySpace-1",
    "page": "Spaces",
    "title": "ArraySpace",
    "category": "section",
    "text": "ArraySpace(::Array{<:Space}) represents the direct sum of the spaces, where evaluation is defined in an array-wise manner. A simple example is the following:julia> x = Fun(identity,-1..1);\n\njulia> f = [exp(x); sqrt(1-x^2)*cos(x-0.1)];\n\njulia> space(f)\n2-element ArraySpace:\n Chebyshev(【-1.0,1.0】)             \n (1-x^2)^0.5[Chebyshev(【-1.0,1.0】)]\n\njulia> a, b = components(f);\n\njulia> norm(f(0.5) - [a(0.5); b(0.5)])\n0.0\n\njulia> norm(f.coefficients[1:2:end] - a.coefficients)\n0.0\n\njulia> norm(f.coefficients[2:2:end] - b.coefficients)\n0.0More complicated examples may interlace the coefficients using a different strategy."
},

{
    "location": "usage/spaces.html#TensorSpace-1",
    "page": "Spaces",
    "title": "TensorSpace",
    "category": "section",
    "text": "TensorSpace((space_1,space_2)) represents the tensor product of two spaces. See documentation of TensorSpace for more details on how the coefficients are interlaced. Note that more than two spaces is only partially supported."
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
    "text": "The default space is Chebyshev, which can represent non-periodic functions on intervals.  Each Space type has a default domain: for Chebyshev this is -1..1, for Fourier and Laurent this is -π..π.  Thus the following are synonyms:Fun(exp, Chebyshev(Interval(-1,1)))\nFun(exp, Chebyshev(Interval()))\nFun(exp, Chebyshev(-1..1))\nFun(exp, Chebyshev())\nFun(exp, -1..1)\nFun(exp, Interval())\nFun(exp, Interval(-1,1))\nFun(exp)If a function is not specified, then it is taken to be identity.  Thus we have the following synonyms:x = Fun(identity, -1..1)\nx = Fun(-1..1)\nx = Fun(identity)\nx = Fun()"
},

{
    "location": "usage/constructors.html#Specifying-coefficients-explicitly-1",
    "page": "Constructors",
    "title": "Specifying coefficients explicitly",
    "category": "section",
    "text": "It is sometimes necessary to specify coefficients explicitly.  This is possible via specifying the space followed by a vector of coefficients:julia> f = Fun(Taylor(), [1,2,3]);  # represents 1 + 2z + 3z^2\n\njulia> f(0.1)\n1.23\n\njulia> 1 + 2*0.1 + 3*0.1^2\n1.23In higher dimensions, ApproxFun will sum products of the 1D basis functions. So if T_i(x) is the ith basis function, then a 2D function can be approximated as the following: f(x  y) = sum_i j c_ij  T_i(x)  T_j(y)The products will be ordered lexicographically by the degree of the polynomial, i.e. in the order T_0(x)  T_0(y)  T_0(x)  T_1(y)   T_1(x)  T_0(y)   T_0(x)  T_2(y)   T_1(x)  T_1(y)   T_2(x)  T_0(y)    . For example, if we are in the two dimensional CosSpace space and we have coefficients c_1 c_2 c_3, then $ f(x, y) = c1 \\cos(0 x) \\cos(0 y) + c2 \\cos(0 x) \\cos(1 y) + c_3 \\cos(1 x) \\cos(0 y). $This is illustrated in the following code:julia> f = Fun(CosSpace()^2, [1,2,3])\nFun(CosSpace(【0.0,6.283185307179586❫)⊗CosSpace(【0.0,6.283185307179586❫),[1.0,2.0,3.0])\n\njulia> f(1,2)\n1.7886132445101346\n\njulia> 1cos(0*1)*cos(0*2) + 2cos(0*1)*cos(1*2) + 3cos(1*1)*cos(0*2)\n1.7886132445101346"
},

{
    "location": "usage/constructors.html#Using-ApproxFun-for-“manual”-interpolation-1",
    "page": "Constructors",
    "title": "Using ApproxFun for “manual” interpolation",
    "category": "section",
    "text": "The ApproxFun package for Julia implements all of the necessary operations for Chebyshev interpolation and operations (like differentiation or integration) on Chebyshev interpolants.Normally, you give it a function f and a domain d, and construct the Chebyshev interpolant by fc = Fun(f, d). The ApproxFun package figures out the necessary number of Chebyshev points (i.e., the polynomial order) required to interpolate f to nearly machine precision, so that subsequent operations on fc can be viewed as \"exact\".However, in cases where the function to be interpolated is extremely expensive, and possibly even is evaluated by an external program, it is convenient to be able to decide on the desired Chebyshev order in advance, evaluate the function at those points \"manually\", and then construct the Chebyshev interpolant. An example showing how to do this is given in the ApproxFun FAQ.DocTestSetup = nothing"
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
    "text": "Differential and integral operators are perhaps the most useful type of operators in mathematics.  Consider the derivative operator on CosSpace:julia> D = Derivative(CosSpace())\nConcreteDerivative:CosSpace(【0.0,6.283185307179586❫)→SinSpace(【0.0,6.283185307179586❫)\n 0.0  -1.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅    0.0  -2.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅     ⋅    0.0  -3.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅     ⋅     ⋅    0.0  -4.0    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅    0.0  -5.0    ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -6.0    ⋅     ⋅     ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -7.0    ⋅     ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -8.0    ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  -9.0  ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.0  ⋱\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱\n\njulia> f = Fun(θ->cos(cos(θ)), CosSpace());\n\njulia> fp = D*f;\n\njulia> fp(0.1) ≈ f\'(0.1) ≈ sin(cos(0.1))*sin(0.1)\ntrueHere, we specified the domain space for the derivative operator, and it automatically determined the range space:DocTestSetup = quote\n    using ApproxFun\n    D = Derivative(CosSpace())\n    f = Fun(θ->cos(cos(θ)),CosSpace())\n    fp = D*f\nendjulia> rangespace(D) == space(fp) == SinSpace()\ntrueOperators can be identified with infinite-dimensional matrices, whose entries are given by the canonical bases in the domain and range space.  In this case, the relevant formula is D cos k theta = -k sin k theta That is, the (k,k+1)th entry is as follows:julia> k,j = 5,6;\n\njulia> ej = Fun(domainspace(D),[zeros(j-1);1]);\n\njulia> D[k,j] ≈ (D*ej).coefficients[k] ≈ -k\ntrueThe Chebyshev space has the property that its derivatives are given by ultraspherical spaces:julia> Derivative(Chebyshev())\nConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n ⋅  1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅   2.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅   3.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅   4.0   ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅   5.0   ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅   6.0   ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   7.0   ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   8.0   ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   9.0  ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱"
},

{
    "location": "usage/operators.html#Functionals-1",
    "page": "Operators",
    "title": "Functionals",
    "category": "section",
    "text": "A particularly useful class of operators are functionals, which map from functions to scalar numbers.  These are represented by operators of size 1 × ∞: that is, infinite-dimensional analogues of row vectors.As an example, the evaluation functional f(0) on CosSpace has the form:julia> B = Evaluation(CosSpace(),0)\nConcreteEvaluation:CosSpace(【0.0,6.283185307179586❫)→ConstantSpace(Point(0))\n 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  ⋯\n\njulia> B*f ≈ f(0)\ntrueAs can be seen from the output, rangespace(B) is a ConstantSpace(Point(0)), a one-dimensional space used to represent scalars whose domain is a single point, 0.Closely related to functionals are operators with finite-dimensional range. For example, the Dirichlet operator represents the restriction of a space to its boundary. In the case, of Chebyshev(), this amounts to evaluation at the endpoints ±1:julia> B = Dirichlet(Chebyshev())\nConcreteDirichlet:Chebyshev(【-1.0,1.0】)→2-element ArraySpace:\n ConstantSpace(Point(-1.0))\n ConstantSpace(Point(1.0))\n 1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  ⋯\n 1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  ⋯\n\njulia> size(B)\n(2, ∞)\n\njulia> B*Fun(exp)\nFun(2-element ArraySpace:\n ConstantSpace(Point(-1.0))\n ConstantSpace(Point(1.0)) ,[0.367879, 2.71828])\n\njulia> B*Fun(exp) ≈ Fun([exp(-1),exp(1)])\ntrue"
},

{
    "location": "usage/operators.html#Multiplication-1",
    "page": "Operators",
    "title": "Multiplication",
    "category": "section",
    "text": "A Multiplication operator sends a Fun to a Fun in the corresponding space by multiplying a given function. The Multiplication operators are presented in matrix form in ApproxFun.julia> x = Fun();\n\njulia> M = Multiplication(1 + 2x + x^2, Chebyshev())\nConcreteMultiplication:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)\n 1.5  1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    ⋅\n 2.0  1.75  1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    ⋅\n 0.5  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅    ⋅\n  ⋅   0.25  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅     ⋅    ⋅\n  ⋅    ⋅    0.25  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅    ⋅\n  ⋅    ⋅     ⋅    0.25  1.0   1.5   1.0   0.25   ⋅     ⋅    ⋅\n  ⋅    ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   0.25   ⋅    ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   0.25  ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   ⋱\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   ⋱\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋱     ⋱    ⋱\n\njulia> (M * x).coefficients == ((1 + 2x + x^2) * x).coefficients == M[1:4,1:2] * x.coefficients\ntrueIt is possible for domain space and range space to be different under Mulitplication.julia> c = Fun(θ -> cos(θ), CosSpace());\n\njulia> Multiplication(c, SinSpace())\nConcreteMultiplication:SinSpace(【0.0,6.283185307179586❫)→SinSpace(【0.0,6.283185307179586❫)\n 8.974302782386682e-17  0.5                    …   ⋅                     ⋅\n 0.5                    8.974302782386682e-17      ⋅                     ⋅\n  ⋅                     0.5                        ⋅                     ⋅\n  ⋅                      ⋅                         ⋅                     ⋅\n  ⋅                      ⋅                         ⋅                     ⋅\n  ⋅                      ⋅                     …   ⋅                     ⋅\n  ⋅                      ⋅                         ⋅                     ⋅\n  ⋅                      ⋅                         ⋅                     ⋅\n  ⋅                      ⋅                        0.5                    ⋅\n  ⋅                      ⋅                        8.974302782386682e-17  ⋱\n  ⋅                      ⋅                     …   ⋱                     ⋱If a function is given by the expansion $ f(\\theta) = \\sum{n=1}^{\\infty}  {f}{n} * sin(n\\theta) $Then the matrix above can be easily derived from $ cos(\\theta) * f(\\theta) = cos(\\theta) \\cdot (\\sum{n=1}^{\\infty}  {f}{n} \\cdot sin(n\\theta) $ $ = \\sum{n=1}^{\\infty}  {f}{n} \\cdot cos(\\theta) \\cdot sin(n\\theta) $ $ = \\sum{n=1}^{\\infty}  {f}{n} \\cdot 0.5 \\cdot ((sin(n-1)\\theta) + (sin(n+1)\\theta) $ $ = \\sum{n=1}^{\\infty}  0.5 \\cdot ({f}{n-1} + {f}_{n+1}) \\cdot sin(n\\theta) $."
},

{
    "location": "usage/operators.html#Algebraic-manipulation-of-operators-1",
    "page": "Operators",
    "title": "Algebraic manipulation of operators",
    "category": "section",
    "text": "Operators can be algebraically manipulated, provided that the domain and range spaces are compatible, or can be made compatible.   As a simple example, we can add the second derivative of a Fourier space to the identity operator:julia> D2 = Derivative(Fourier(),2)\nDerivativeWrapper:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)\n 0.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅   -1.0    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅     ⋅   -1.0    ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅     ⋅     ⋅   -4.0    ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅   -4.0    ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅   -9.0    ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -9.0     ⋅      ⋅      ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -16.0     ⋅      ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅   -16.0     ⋅   ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅   -25.0  ⋅\n  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋱\n\njulia> D2 + I\nPlusOperator:Fourier(【0.0,6.283185307179586❫)→Fourier(【0.0,6.283185307179586❫)\n 1.0   ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅   0.0   ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅    ⋅   0.0    ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅    ⋅    ⋅   -3.0    ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅    ⋅    ⋅     ⋅   -3.0    ⋅     ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅    ⋅    ⋅     ⋅     ⋅   -8.0    ⋅      ⋅      ⋅      ⋅   ⋅\n  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅   -8.0     ⋅      ⋅      ⋅   ⋅\n  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅   -15.0     ⋅      ⋅   ⋅\n  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅   -15.0     ⋅   ⋅\n  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅   -24.0  ⋅\n  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋱When the domain and range space are not the same, the identity operator becomes a conversion operator.  That is, to represent D+I acting on the Chebyshev space, we would do the following:julia> D = Derivative(Chebyshev())\nConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n ⋅  1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅   2.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅   3.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅   4.0   ⋅    ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅   5.0   ⋅    ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅   6.0   ⋅    ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   7.0   ⋅    ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   8.0   ⋅   ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   9.0  ⋅\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱\n ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱\n\njulia> C = Conversion(Chebyshev(),Ultraspherical(1))\nConcreteConversion:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n 1.0  0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅   0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5  ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  ⋱\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5  ⋱\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱\n\n\njulia> D + C\nPlusOperator:Chebyshev(【-1.0,1.0】)→Ultraspherical(1,【-1.0,1.0】)\n 1.0  1.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅   0.5   2.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅    0.5   3.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅    0.5   4.0  -0.5    ⋅     ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅    0.5   5.0  -0.5    ⋅     ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅    0.5   6.0  -0.5    ⋅     ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.5   7.0  -0.5    ⋅   ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   8.0  -0.5  ⋅\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   9.0  ⋱\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5  ⋱\n  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱ApproxFun can automatically determine the spaces, so if one writes D + I it will translate it to D + C.  Now consider the Fredholm integral operator of the second kind:L u = u + rm e^x int_-1^1 u(x) rm dxWe can construct this usingjulia> x = Fun();\n\njulia> Σ = DefiniteIntegral(Chebyshev())\nConcreteDefiniteIntegral:Chebyshev(【-1.0,1.0】)→ConstantSpace\n 2.0  0.0  -0.6666666666666666  0.0  …  0.0  -0.031746031746031744  0.0  ⋯\n\njulia> L = I + exp(x)*Σ\nLowRankPertOperator:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)\n 3.5321317555040164     0.0  …  -0.0401925675476828      0.0  ⋯\n 2.260636415969941      1.0     -0.03588311771380859     0.0  ⋱\n 0.5429906790681531     0.0     -0.008618899667748462    0.0  ⋱\n 0.08867369969732766    0.0     -0.0014075190428147247   0.0  ⋱\n 0.010948480884187475   0.0     -0.00017378541086011864  0.0  ⋱\n 0.001085852623827888   0.0  …  -1.7235755933775998e-5   0.0  ⋱\n 8.995464590859038e-5   0.0     -1.4278515223585775e-6   0.0  ⋱\n 6.396872924803984e-6   0.0     -1.015376654730791e-7    0.0  ⋱\n 3.9842496133455937e-7  0.0      0.9999999936757943      0.0  ⋱\n 2.20735434510347e-8    0.0     -3.503737055719794e-10   1.0  ⋱\n  ⋮                      ⋱   …    ⋱                       ⋱   ⋱\n\njulia> u = cos(10x^2);\n\njulia> (L*u)(0.1)\n1.3777980523127336\n\njulia> u(0.1) + exp(0.1)*sum(u)\n1.3777980523127336Note that DefiniteIntegral is a functional, i.e., a 1 × ∞ operator.  when multiplied on the left by a function, it automatically constructs the operator rm e^x int_-1^1 f(x) dx viaDocTestSetup = quote\n    using ApproxFun\n    x = Fun()\n    Q = DefiniteIntegral(Chebyshev())\nendjulia> M = Multiplication(exp(x),ConstantSpace())\nConcreteMultiplication:ConstantSpace→Chebyshev(【-1.0,1.0】)\n 1.26607    \n 1.13032    \n 0.271495   \n 0.0443368  \n 0.00547424\n 0.000542926\n 4.49773e-5\n 3.19844e-6\n 1.99212e-7\n 1.10368e-8\n  ⋮         \n\njulia> M*Σ\nTimesOperator:Chebyshev(【-1.0,1.0】)→Chebyshev(【-1.0,1.0】)\n 2.53213     0.0  -0.844044     0.0  …  0.0  -0.0401926    0.0  ⋯\n 2.26064     0.0  -0.753545     0.0     0.0  -0.0358831    0.0  ⋱\n 0.542991    0.0  -0.180997     0.0     0.0  -0.0086189    0.0  ⋱\n 0.0886737   0.0  -0.0295579    0.0     0.0  -0.00140752   0.0  ⋱\n 0.0109485   0.0  -0.00364949   0.0     0.0  -0.000173785  0.0  ⋱\n 0.00108585  0.0  -0.000361951  0.0  …  0.0  -1.72358e-5   0.0  ⋱\n 8.99546e-5  0.0  -2.99849e-5   0.0     0.0  -1.42785e-6   0.0  ⋱\n 6.39687e-6  0.0  -2.13229e-6   0.0     0.0  -1.01538e-7   0.0  ⋱\n 3.98425e-7  0.0  -1.32808e-7   0.0     0.0  -6.32421e-9   0.0  ⋱\n 2.20735e-8  0.0  -7.35785e-9   0.0     0.0  -3.50374e-10  0.0  ⋱\n  ⋮           ⋱     ⋱            ⋱   …   ⋱     ⋱            ⋱   ⋱Note that Q*exp(x) applies the operator to a function.  To construct the operator that first multiplies by exp(x), use Q[exp(x)].  This is equivalent to Q*Multiplication(exp(x),Chebyshev())."
},

{
    "location": "usage/operators.html#Operators-and-space-promotion-1",
    "page": "Operators",
    "title": "Operators and space promotion",
    "category": "section",
    "text": "It is often more convenient to not specify a space explicitly, but rather infer it when the operator is used.  For example, we can construct Derivative(), which has the alias 𝒟, and represents the first derivative on any space:julia> f = Fun(cos,Chebyshev(0..1)); (𝒟*f)(0.1)\n-0.09983341664681705\n\njulia> f = Fun(cos,Fourier()); (𝒟*f)(0.1)\n-0.09983341664682804Behind the scenes, Derivative() is equivalent to Derivative(UnsetSpace(),1). When multiplying a function f, the domain space is promoted before multiplying, that is, Derivative()*f is equivalent to Derivative(space(f))*f.  This promotion of the domain space happens even when operators have spaces attached. This facilitates the following construction:julia> D = Derivative(Chebyshev());\n\njulia> D^2\nConcreteDerivative:Chebyshev(【-1.0,1.0】)→Ultraspherical(2,【-1.0,1.0】)\n ⋅  ⋅  4.0   ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n ⋅  ⋅   ⋅   6.0   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n ⋅  ⋅   ⋅    ⋅   8.0    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅\n ⋅  ⋅   ⋅    ⋅    ⋅   10.0    ⋅     ⋅     ⋅     ⋅   ⋅\n ⋅  ⋅   ⋅    ⋅    ⋅     ⋅   12.0    ⋅     ⋅     ⋅   ⋅\n ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅   14.0    ⋅     ⋅   ⋅\n ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅   16.0    ⋅   ⋅\n ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅   18.0  ⋅\n ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱\n ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱\n ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱Note that rangespace(D) ≠ Chebyshev(), hence the operators are not compatible. Therefore, it has thrown away its domain space, and thus this is equivalent to Derivative(rangespace(D))*D.DocTestSetup = nothing"
},

{
    "location": "usage/operators.html#Concatenating-operators-1",
    "page": "Operators",
    "title": "Concatenating operators",
    "category": "section",
    "text": "The concatenation functions vcat, hcat and hvcat are overriden for operators to represent the resulting combined operator, now with a rangespace or domainspace that is an ArraySpace."
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
    "text": "Linear equations such as ordinary and partial differential equations,  fractional differential equations and integral equations can be solved using ApproxFun. This is accomplished using A\\b where A is an Operator and b is a Fun.  As a simple example, consider the equationu(theta) + cu(theta) = costhetawhere we want a solution that is periodic on 02pi).  This can be solved succinctly as follows:DocTestSetup = quote\n    using ApproxFun\nendjulia> b = Fun(cos,Fourier());\n\njulia> c = 0.1; u = (𝒟+c*I)\\b;\n\njulia> u(0.6)\n0.64076835137228\n\njulia> (c*cos(0.6)+sin(0.6))/(1+c^2)  # exact solution\n0.6407683513722804Recall that 𝒟 is an alias to Derivative() == Derivative(UnsetSpace(),1).As another example, consider the Fredholm integral equationu + rm e^x int_-1^1 cos x  u(x) rm dx = cos rm e^xWe can solve this equation as follows:julia> Σ = DefiniteIntegral(Chebyshev()); x = Fun();\n\njulia> u = (I+exp(x)*Σ[cos(x)]) \\ cos(exp(x));\n\njulia> u(0.1)\n0.21864294855628802Note that we used the syntax op[f::Fun], which is a shorthand for op*Multiplication(f)."
},

{
    "location": "usage/equations.html#Boundary-conditions-1",
    "page": "Linear Equations",
    "title": "Boundary conditions",
    "category": "section",
    "text": "Incorporating boundary conditions into differential equations is important so that the equation is well-posed.  This is accomplished via combining operators and functionals (i.e., 1 × ∞ operators).  As a simple example, consider the first order initial value problemu = t u qquadhboxandqquad u(0) = 1To pose this in ApproxFun, we want to find a u such that Evaluation(0)*u == 1 and (𝒟 - t)*u == 0.  This is accomplished via:julia> t = Fun(0..1);\n\njulia> u = [Evaluation(0); 𝒟 - t]  \\ [1;0];\n\njulia> u(0)\n0.9999999999999996\n\njulia> norm(u\'-t*u)\n1.2016080299388273e-16Behind the scenes, the Vector{Operator{T}} representing the functionals and operators are combined into a single InterlaceOperator.A common usage is two-point boundary value problems. Consider the singularly perturbed boundary value problem:epsilon u-xu+u = u qquad u(-1) = 1quad u(1) = 2This can be solved in ApproxFun via:julia> x = Fun();\n\njulia> u = [Evaluation(-1);\n            Evaluation(1);\n            1/70*𝒟^2-x*𝒟+I] \\ [1,2,0];\n\njulia> u(0.1)\n0.049999999999960326Note in this case the space is inferred from the variable coefficient x.This ODE can also be solved using the Dirichlet operator:julia> x = Fun();\n\njulia> u = [Dirichlet();\n            1/70*𝒟^2-x*𝒟+I] \\ [[1,2],0];\n\njulia> u(0.1)\n0.04999999999996019"
},

{
    "location": "usage/equations.html#Eigenvalue-Problems-1",
    "page": "Linear Equations",
    "title": "Eigenvalue Problems",
    "category": "section",
    "text": "In analogy to linear algebra, many differential equations may be posed as eigenvalue problems. That is, for some differential operator L, there are a family of functions u_i(x) such that $ L~ui(x) = \\lambdai ui(x) $ where \\lambdai$ is the i^th eigenvalue of the L and has a corresponding eigenfunction u_i(x). A classic eigenvalue problem is known as the quantum harmonic oscillator where L = -frac12fracd^2dx^2 + frac12 x^2 and one demands that u(infty) = u(-infty) = 0. Because we expect the solutions to be exponentially suppressed for large x, we can approximate this with Dirichlet boundary conditions at a \'reasonably large\' x without much difference.We can express this in ApproxFun as the following:x = Fun(-8 .. 8)\nL = -𝒟^2/2 + x^2/2\nS = space(x)\nB = Dirichlet(S)\nλ, v = eigs(B, L, 500,tolerance=1E-10)note that boundary conditions must be specified in the call to eigs. Plotting the first 20 eigenfunctions offset vertically by their eigenvalue, we see(Image: harmonic_eigs)If the solutions are not relatively constant near the boundary then one should push the boundaries further out.For problems with different contraints or boundary conditions, B can be any zero functional constraint, eg. DefiniteIntegral()."
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
    "text": "Behind the scenes, A\\b where A is an Operator is implemented via an adaptive QR factorization.  That is, it is equivalent to qr(A)\\b.  (There is a subtly here in space inferring: A\\b can use     both A and b to determine the domain space, while qr(A) only     sees the operator A.)       Note that qr adaptively caches a partial QR Factorization as it is applied to different right-hand sides, so the same operator can be inverted much more efficiently in subsequent problems."
},

{
    "location": "usage/equations.html#Partial-differential-equations-1",
    "page": "Linear Equations",
    "title": "Partial differential equations",
    "category": "section",
    "text": "Partial differential operators are also supported.  Here\'s an example of solving the Poisson equation with zero boundary conditions:d = Domain(-1..1)^2\nx,y = Fun(d)\nf = exp.(-10(x+0.3)^2-20(y-0.2)^2)  # use broadcasting as exp(f) not implemented in 2D\nA = [Dirichlet(d);Δ]              # Δ is an alias for Laplacian()\n@time u = A \\ [zeros(∂(d));f]     #4s for ~3k coefficientsUsing a QR Factorization reduces the cost of subsequent calls substantially:QR = qr(A)\n@time QR \\ [zeros(∂(d));f]   # 4s\ng = exp.(-10(x+0.2)^2-20(y-0.1)^2)\n@time QR \\ [zeros(∂(d));g]  # 0.09sMany PDEs have weak singularities at the corners, in which case it is beneficial to specify a tolerance to reduce the time:\\(A, [zeros(∂(d));f]; tolerance=1E-6)"
},

{
    "location": "usage/equations.html#Nonlinear-equations-1",
    "page": "Linear Equations",
    "title": "Nonlinear equations",
    "category": "section",
    "text": "There is preliminary support for nonlinear equations, via Newton iteration in function space.  Here is a simple two-point boundary value problem:beginalign*\n    epsilon u + 6(1-x^2)u +u^2=1 cr\n    u(-1)=u(1)=0\nendalign*This can be solved usingx = Fun()\nN = u -> [u(-1.)-c; u(1.); ε*u\'\' + 6*(1-x^2)*u\' + u^2 - 1.0]\nu = newton(N,u0)"
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
    "text": "In the case where the grid is specified by points(space,n), you can apply the default transform to data:DocTestSetup = quote\n    using ApproxFun\nendjulia> S = Chebyshev(1..2);\n\njulia> p = points(S,20); # the default grid\n\njulia> v = exp.(p);      # values at the default grid\n\njulia> f = Fun(S,ApproxFun.transform(S,v));\n\njulia> f(1.1)\n3.0041660239464347\n\njulia> exp(1.1)\n3.0041660239464334ApproxFun has no inbuilt support for interpolating functions at other sets of points, but this can be accomplished manually by evaluating the basis at the set of points and using \\:julia> S = Chebyshev(1..2);\n\njulia> n = 50;\n\njulia> p = range(1,stop=2,length=n);   # a non-default grid\n\njulia> v = exp.(p);           # values at the non-default grid\n\njulia> V = Array(Float64,n,n); # Create a Vandermonde matrix by evaluating the basis at the grid\n\njulia> for k = 1:n\n           V[:,k] = Fun(S,[zeros(k-1);1]).(p)\n       end\n\njulia> f = Fun(S,V\\v);\n\njulia> f(1.1)\n3.0041660228311926\n\njulia> exp(1.1)\n3.0041660239464334Note that an evenly spaced grid suffers from instability for large n.  The easiest way around this is to use least squares with more points than coefficients, instead of interpolation:julia> S = Chebyshev(1..2);\n\njulia> n = 100; m = 50;\n\njulia> p = range(1,stop=2,length=n);   # a non-default grid\n\njulia> v = exp.(p);           # values at the non-default grid\n\njulia> V = Array(Float64,n,m); # Create a Vandermonde matrix by evaluating the basis at the grid\n\njulia> for k = 1:m\n           V[:,k] = Fun(S,[zeros(k-1);1]).(p)\n       end\n\njulia> f = Fun(S,V\\v);\n\njulia> f(1.1)\n3.004166023946434\n\njulia> exp(1.1)\n3.0041660239464334We can use this same approach for multivariate functions:julia> S = Chebyshev(0..1)^2;\n\njulia> n = 1000; m = 50;\n\njulia> Random.seed!(0); x = rand(n); y = rand(n);\n\njulia> v = exp.(x .* cos(y));  # values at the non-default grid\n\njulia> V = Array(Float64,n,m); # Create a Vandermonde matrix by evaluating the basis at the grid\n\njulia> for k = 1:m\n          V[:,k] = Fun(S,[zeros(k-1);1]).(x,y)\n       end\n\n\njulia> f = Fun(S,V\\v);\n\njulia> f(0.1,0.2)\n1.1029700685084018\n\njulia> exp(0.1*cos(0.2))\n1.1029701284210731DocTestSetup = nothing"
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
    "category": "type",
    "text": "GaussWeight(Hermite(L), L) is a space spanned by exp(-Lx²) * H_k(sqrt(L) * x) where H_k(x)\'s are Hermite polynomials.\n\nGaussWeight() is equivalent to GaussWeight(Hermite(), 1.0) by default.\n\n\n\n\n\nFun(s::Space,coefficients::AbstractVector)\n\nreturns a Fun with the specified coefficients in the space s\n\n\n\n\n\nFun(f,s::Space)\n\nreturn a Fun representing the function, number, or vector f in the space s.  If f is vector-valued, it returns a vector-valued analogue of s.\n\n\n\n\n\nFun(f,d::Domain)\n\nreturns Fun(f,Space(d)), that is, it uses the default space for the specified domain.\n\n\n\n\n\nFun(s::Space)\n\nreturns Fun(identity,s)\n\n\n\n\n\nFun(f)\n\nreturns Fun(f,Chebyshev())\n\n\n\n\n\nFun()\n\nreturns Fun(identity,Chebyshev()).\n\n\n\n\n\n"
},

{
    "location": "library.html#Base.ones-Tuple{Space}",
    "page": "Library",
    "title": "Base.ones",
    "category": "method",
    "text": "ones(d::Space)\n\nReturn the Fun that represents the function one on the specified space.\n\n\n\n\n\n"
},

{
    "location": "library.html#Base.zeros-Tuple{Space}",
    "page": "Library",
    "title": "Base.zeros",
    "category": "method",
    "text": "zeros(d::Space)\n\nReturn the Fun that represents the function one on the specified space.\n\n\n\n\n\n"
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
    "category": "type",
    "text": "Arc(c,r,(θ₁,θ₂))\n\nrepresents the arc centred at c with radius r from angle θ₁ to θ₂.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Circle",
    "page": "Library",
    "title": "ApproxFun.Circle",
    "category": "type",
    "text": "Circle(c,r,o)\n\nrepresents the circle centred at c with radius r which is positively (o=true) or negatively (o=false) oriented.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Curve",
    "page": "Library",
    "title": "ApproxFun.Curve",
    "category": "constant",
    "text": "Curve Represents a domain defined by the image of a Fun.  Example usage would be\n\nx=Fun(1..2)\nCurve(exp(im*x))  # represents an arc\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Disk",
    "page": "Library",
    "title": "ApproxFun.Disk",
    "category": "type",
    "text": "Disk(c,r)\n\nrepresents the disk centred at c with radius r.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Segment",
    "page": "Library",
    "title": "ApproxFun.Segment",
    "category": "type",
    "text": "Segment(a,b)\n\nrepresents a line segment from a to b.  In the case where a and b are real and a < b, then this is is equivalent to an Interval(a,b).\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Interval",
    "page": "Library",
    "title": "ApproxFun.Interval",
    "category": "function",
    "text": "Interval(a::Real,b::Real)\n\nrepresents the set {x : a ≤ x ≤ b}.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Line",
    "page": "Library",
    "title": "ApproxFun.Line",
    "category": "type",
    "text": "Line{a}(c)\n\nrepresents the line at angle a in the complex plane, centred at c.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.PeriodicInterval",
    "page": "Library",
    "title": "ApproxFun.PeriodicInterval",
    "category": "type",
    "text": "PeriodicInterval(a,b)\n\nrepresents a periodic interval from a to b, that is, the point b is identified with a.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Point",
    "page": "Library",
    "title": "ApproxFun.Point",
    "category": "type",
    "text": "Point(x)\n\nrepresents a single point at x.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ProductDomain",
    "page": "Library",
    "title": "ApproxFun.ProductDomain",
    "category": "type",
    "text": "ProductDomain((d1,d2))\n\nrepresents the product of two domains, the set {(x,y) : x ∈ d1 & y ∈ d2}.\n\nMultiplication of domains is overrident to return a ProductDomain. For example, the following represents the rectangle 1 ≤ x ≤ 2 & 3 ≤ y ≤ 4:\n\nInterval(1,2)*(3,4)\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Ray",
    "page": "Library",
    "title": "ApproxFun.Ray",
    "category": "type",
    "text": "Ray{a}(c,o)\n\nrepresents a ray at angle a starting at c, with orientation out to infinity (o = true) or back from infinity (o = false).\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.UnionDomain",
    "page": "Library",
    "title": "ApproxFun.UnionDomain",
    "category": "type",
    "text": "UnionDomain((d1,d2,…,dn))\n\nrepresents a union of multiple subdomains: {x : x ∈ d1 || … || x ∈ dn}.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.∂",
    "page": "Library",
    "title": "ApproxFun.∂",
    "category": "function",
    "text": "∂(d::Domain)\n\nreturns the boundary of d.  For example, the boundary of a Disk() is a Circle(), and the boundary of Interval()^2 is a piecewise interval that sketches the boundary of a rectangle.\n\n\n\n\n\n"
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
    "category": "function",
    "text": "canonicalspace(s::Space)\n\nreturns a space that is used as a default to implement missing functionality, e.g., evaluation.  Implement a Conversion operator or override coefficients to support this.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.itransform",
    "page": "Library",
    "title": "ApproxFun.itransform",
    "category": "function",
    "text": "itransform(s::Space,coefficients::AbstractVector)\n\nTransform coefficients back to values.  Defaults to using canonicalspace as in transform.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.transform",
    "page": "Library",
    "title": "ApproxFun.transform",
    "category": "function",
    "text": "transform(s::Space,vals::Vector)\n\nTransform values on the grid specified by points(s,length(vals)) to coefficients in the space s. Defaults to coefficients(transform(canonicalspace(space),values),canonicalspace(space),space)\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.dimension-Tuple{Space}",
    "page": "Library",
    "title": "ApproxFun.dimension",
    "category": "method",
    "text": "dimension(s::Space)\n\nreturns the dimension of s, which is the maximum number of coefficients.\n\n\n\n\n\n"
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
    "category": "type",
    "text": "SequenceSpace is the space of all sequences, i.e., infinite vectors. Also denoted ℓ⁰.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ConstantSpace",
    "page": "Library",
    "title": "ApproxFun.ConstantSpace",
    "category": "type",
    "text": "ConstantSpace is the 1-dimensional scalar space.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Chebyshev",
    "page": "Library",
    "title": "ApproxFun.Chebyshev",
    "category": "type",
    "text": "Chebyshev() is the space spanned by the Chebyshev polynomials\n\n    T_0(x),T_1(x),T_2(x),…\n\nwhere T_k(x) = cos(k*acos(x)).  This is the default space as there exists a fast transform and general smooth functions on [-1,1] can be easily resolved.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Jacobi",
    "page": "Library",
    "title": "ApproxFun.Jacobi",
    "category": "type",
    "text": "Jacobi(b,a) represents the space spanned by Jacobi polynomials P_k^{(a,b)}, which are orthogonal with respect to the weight (1+x)^β*(1-x)^α\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Ultraspherical",
    "page": "Library",
    "title": "ApproxFun.Ultraspherical",
    "category": "type",
    "text": "Ultraspherical(λ) is the space spanned by the ultraspherical polynomials\n\n    C_0^{(λ)}(x),C_1^{(λ)}(x),C_2^{(λ)}(x),…\n\nNote that λ=1 this reduces to Chebyshev polynomials of the second kind: C_k^{(1)}(x) = U_k(x).\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Taylor",
    "page": "Library",
    "title": "ApproxFun.Taylor",
    "category": "type",
    "text": "Taylor() is the space spanned by [1,z,z^2,...]. This is a type alias for Hardy{true}.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Hardy",
    "page": "Library",
    "title": "ApproxFun.Hardy",
    "category": "type",
    "text": "Hardy{false}() is the space spanned by [1/z,1/z^2,...]. Hardy{true}() is the space spanned by [1,z,z^2,...].\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Fourier",
    "page": "Library",
    "title": "ApproxFun.Fourier",
    "category": "type",
    "text": "Fourier() is the space spanned by the trigonemtric polynomials\n\n    1,sin(θ),cos(θ),sin(2θ),cos(2θ),…\n\nSee also Laurent.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Laurent",
    "page": "Library",
    "title": "ApproxFun.Laurent",
    "category": "type",
    "text": "Laurent() is the space spanned by the complex exponentials\n\n    1,exp(-im*θ),exp(im*θ),exp(-2im*θ),…\n\nSee also Fourier.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.CosSpace",
    "page": "Library",
    "title": "ApproxFun.CosSpace",
    "category": "type",
    "text": "CosSpace() is the space spanned by [1,cos θ,cos 2θ,...]\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.SinSpace",
    "page": "Library",
    "title": "ApproxFun.SinSpace",
    "category": "type",
    "text": "SinSpace() is the space spanned by [sin θ,sin 2θ,...]\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.JacobiWeight",
    "page": "Library",
    "title": "ApproxFun.JacobiWeight",
    "category": "type",
    "text": "JacobiWeight(β,α,s::Space)\n\nweights a space s by a Jacobi weight, which on -1..1 is (1+x)^β*(1-x)^α. For other domains, the weight is inferred by mapping to -1..1.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.LogWeight",
    "page": "Library",
    "title": "ApproxFun.LogWeight",
    "category": "type",
    "text": "LogWeight(β,α,s::Space)\n\nrepresents a function on -1..1 weighted by log((1+x)^β*(1-x)^α). For other domains, the weight is inferred by mapping to -1..1.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ArraySpace",
    "page": "Library",
    "title": "ApproxFun.ArraySpace",
    "category": "type",
    "text": "ArraySpace(s::Space,dims...)\n\nis used to represent array-valued expansions in a space s.  The coefficients are of each entry are interlaced.\n\nFor example,\n\nf = Fun(x->[exp(x),sin(x)],-1..1)\nspace(f) == ArraySpace(Chebyshev(),2)\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.TensorSpace",
    "page": "Library",
    "title": "ApproxFun.TensorSpace",
    "category": "type",
    "text": "TensorSpace(a::Space,b::Space)\n\nrepresents a tensor product of two 1D spaces a and b. The coefficients are interlaced in lexigraphical order.\n\nFor example, consider\n\nFourier()*Chebyshev()  # returns TensorSpace(Fourier(),Chebyshev())\n\nThis represents functions on [-π,π) x [-1,1], using the Fourier basis for the first argument and Chebyshev basis for the second argument, that is, φ_k(x)T_j(y), where\n\nφ_0(x) = 1,\nφ_1(x) = sin x,\nφ_2(x) = cos x,\nφ_3(x) = sin 2x,\nφ_4(x) = cos 2x\n…\n\nBy Choosing (k,j) appropriately, we obtain a single basis:\n\nφ_0(x)T_0(y) (= 1),\nφ_0(x)T_1(y) (= y),\nφ_1(x)T_0(y) (= sin x),\nφ_0(x)T_2(y), …\n\n\n\n\n\n"
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
    "category": "function",
    "text": "domain(f::Fun)\n\nreturns the domain that f is defined on\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.coefficients",
    "page": "Library",
    "title": "ApproxFun.coefficients",
    "category": "function",
    "text": "coefficients(f::Fun) -> Vector\n\nreturns the coefficients of f, corresponding to the space space(f).\n\n\n\n\n\ncoefficients(f::Fun,s::Space) -> Vector\n\nreturns the coefficients of f in the space s, which may not be the same as space(f).\n\n\n\n\n\ncoefficients(cfs::AbstractVector,fromspace::Space,tospace::Space) -> Vector\n\nconverts coefficients in fromspace to coefficients in tospace\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.extrapolate",
    "page": "Library",
    "title": "ApproxFun.extrapolate",
    "category": "function",
    "text": "extrapolate(f::Fun,x)\n\nreturns an extrapolation of f from its domain to x.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.ncoefficients",
    "page": "Library",
    "title": "ApproxFun.ncoefficients",
    "category": "function",
    "text": "ncoefficients(f::Fun) -> Integer\n\nreturns the number of coefficients of a fun\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.points",
    "page": "Library",
    "title": "ApproxFun.points",
    "category": "function",
    "text": "points(f::Fun)\n\nreturns a grid of points that f can be transformed into values and back\n\n\n\n\n\npoints(s::Space,n::Integer)\n\nreturns a grid of approximately n points, for which a transform exists from values at the grid to coefficients in the space s.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.space",
    "page": "Library",
    "title": "ApproxFun.space",
    "category": "function",
    "text": "space(f::Fun)\n\nreturns the space of f\n\n\n\n\n\n"
},

{
    "location": "library.html#Base.values",
    "page": "Library",
    "title": "Base.values",
    "category": "function",
    "text": "values(f::Fun)\n\nreturns f evaluated at points(f)\n\n\n\n\n\n"
},

{
    "location": "library.html#Base.stride-Tuple{Fun}",
    "page": "Library",
    "title": "Base.stride",
    "category": "method",
    "text": "stride(f::Fun)\n\nreturns the stride of the coefficients, checked numerically\n\n\n\n\n\n"
},

{
    "location": "library.html#Accessing-information-about-a-Fun-1",
    "page": "Library",
    "title": "Accessing information about a Fun",
    "category": "section",
    "text": "domaincoefficientsextrapolatencoefficientspointsspaceApproxFun.valuesstride(::Fun)"
},

{
    "location": "library.html#ApproxFun.reverseorientation",
    "page": "Library",
    "title": "ApproxFun.reverseorientation",
    "category": "function",
    "text": "reverseorientation(f::Fun)\n\nreturn f on a reversed orientated contour.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.setdomain",
    "page": "Library",
    "title": "ApproxFun.setdomain",
    "category": "function",
    "text": "setdomain(f::Fun,d::Domain)\n\nreturns f projected onto domain\n\n\n\n\n\n"
},

{
    "location": "library.html#Base.chop",
    "page": "Library",
    "title": "Base.chop",
    "category": "function",
    "text": "chop(f::Fun,tol) -> Fun\n\nreduces the number of coefficients by dropping the tail that is below the specified tolerance.\n\n\n\n\n\n"
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
    "category": "type",
    "text": "Operator{T}\n\nis an abstract type to represent linear operators between spaces.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.domainspace",
    "page": "Library",
    "title": "ApproxFun.domainspace",
    "category": "function",
    "text": "domainspace(op::Operator)\n\ngives the domain space of op.  That is, op*f will first convert f to a Fun in the space domainspace(op) before applying the operator.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.rangespace",
    "page": "Library",
    "title": "ApproxFun.rangespace",
    "category": "function",
    "text": "rangespace(op::Operator)\n\ngives the range space of op.  That is, op*f will return a Fun in the space rangespace(op), provided f can be converted to a Fun in domainspace(op).\n\n\n\n\n\n"
},

{
    "location": "library.html#Base.getindex-Tuple{Operator,Fun}",
    "page": "Library",
    "title": "Base.getindex",
    "category": "method",
    "text": "op[f::Fun]\n\nconstructs the operator op*Multiplication(f), that is, it multiplies on the right by f first.  Note that op*f is different: it applies op to f.\n\n\n\n\n\n"
},

{
    "location": "library.html#Base.qr-Tuple{Operator}",
    "page": "Library",
    "title": "Base.qr",
    "category": "method",
    "text": "qr(A::Operator)\n\nreturns a cached QR factorization of the Operator A.  The result QR enables solving of linear equations: if u=QR, then u approximately satisfies A*u = b.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.cache-Tuple{Operator}",
    "page": "Library",
    "title": "ApproxFun.cache",
    "category": "method",
    "text": "cache(op::Operator)\n\nCaches the entries of an operator, to speed up multiplying a Fun by the operator.\n\n\n\n\n\n"
},

{
    "location": "library.html#Operators-1",
    "page": "Library",
    "title": "Operators",
    "category": "section",
    "text": "OperatorBandedMatrices.bandwidths(::Operator)domainspacerangespacegetindex(::Operator,::,::)getindex(::Operator,::Fun)\\(::Operator,::)qr(::Operator)cache(::Operator)"
},

{
    "location": "library.html#ApproxFun.Conversion",
    "page": "Library",
    "title": "ApproxFun.Conversion",
    "category": "type",
    "text": "Conversion(fromspace::Space,tospace::Space)\n\nrepresents a conversion operator between fromspace and tospace, when available.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Derivative",
    "page": "Library",
    "title": "ApproxFun.Derivative",
    "category": "type",
    "text": "Derivative(sp::Space,k::Int) represents the k-th derivative on sp.\n\n\n\n\n\nDerivative(sp::Space,k::Vector{Int}) represents a partial derivative on a multivariate space. For example,\n\nDx = Derivative(Chebyshev()^2,[1,0]) # ∂/∂x\nDy = Derivative(Chebyshev()^2,[0,1]) # ∂/∂y\n\n\n\n\n\nDerivative(sp::Space) represents the first derivative Derivative(sp,1).\n\n\n\n\n\nDerivative(k) represents the k-th derivative, acting on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n\n\nDerivative() represents the first derivative on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Dirichlet",
    "page": "Library",
    "title": "ApproxFun.Dirichlet",
    "category": "type",
    "text": "Dirichlet(sp,k) is the operator associated with restricting the k-th derivative on the boundary for the space sp.\n\n\n\n\n\nDirichlet(sp) is the operator associated with restricting the  the boundary for the space sp.\n\n\n\n\n\nDirichlet() is the operator associated with restricting on the  the boundary.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Evaluation",
    "page": "Library",
    "title": "ApproxFun.Evaluation",
    "category": "type",
    "text": "Evaluation(sp,x,k) is the functional associated with evaluating the k-th derivative at a point x for the space sp.\n\n\n\n\n\nEvaluation(sp,x) is the functional associated with evaluating at a point x for the space sp.\n\n\n\n\n\nEvaluation(x) is the functional associated with evaluating at a point x.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Integral",
    "page": "Library",
    "title": "ApproxFun.Integral",
    "category": "type",
    "text": "Integral(sp::Space,k::Int) represents a k-th integral on sp. There is no guarantee on normalization.\n\n\n\n\n\nIntegral(sp::Space) represents the first integral Integral(sp,1).\n\n\n\n\n\nIntegral(k)represents thek`-th integral, acting on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n\n\nIntergral() represents the first integral on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Laplacian",
    "page": "Library",
    "title": "ApproxFun.Laplacian",
    "category": "type",
    "text": "Laplacian(sp::Space) represents the laplacian on space sp.\n\n\n\n\n\nLaplacian() represents the laplacian on an unset space. Spaces will be inferred when applying or manipulating the operator.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Multiplication",
    "page": "Library",
    "title": "ApproxFun.Multiplication",
    "category": "type",
    "text": "Multiplication(f::Fun,sp::Space) is the operator representing multiplication by f on functions in the space sp.\n\n\n\n\n\nMultiplication(f::Fun) is the operator representing multiplication by f on an unset space of functions.  Spaces will be inferred when applying or manipulating the operator.\n\n\n\n\n\n"
},

{
    "location": "library.html#ApproxFun.Neumann",
    "page": "Library",
    "title": "ApproxFun.Neumann",
    "category": "function",
    "text": "Neumann(sp) is the operator associated with restricting the normal derivative on the boundary for the space sp. At the moment it is implemented as Dirichlet(sp,1).\n\n\n\n\n\n`Neumann( is the operator associated with restricting the normal derivative on the boundary.\n\n\n\n\n\n"
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
    "category": "type",
    "text": "LowRankFun gives an approximation to a bivariate function in low rank form.\n\n\n\n\n\n"
},

{
    "location": "library.html#Bivariate-1",
    "page": "Library",
    "title": "Bivariate",
    "category": "section",
    "text": "LowRankFun"
},

]}
