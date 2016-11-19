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
    "text": ""
},

{
    "location": "index.html#ApproxFun.Interval",
    "page": "Home",
    "title": "ApproxFun.Interval",
    "category": "Type",
    "text": "Interval(a,b)\n\nrepresents an interval from a to b.  In the case where a and b are complex or 2-dimensional, it represents the line segment between a and b.\n\n\n\n"
},

{
    "location": "index.html#Domains-1",
    "page": "Home",
    "title": "Domains",
    "category": "section",
    "text": "Interval"
},

{
    "location": "index.html#ApproxFun.Fun",
    "page": "Home",
    "title": "ApproxFun.Fun",
    "category": "Type",
    "text": "Fun(coefficients,space)\n\nreturns a Fun with coefficients in the space\n\n\n\n"
},

{
    "location": "index.html#Constructing-a-Fun-1",
    "page": "Home",
    "title": "Constructing a Fun",
    "category": "section",
    "text": "Fun"
},

{
    "location": "index.html#ApproxFun.domain",
    "page": "Home",
    "title": "ApproxFun.domain",
    "category": "Function",
    "text": "domain(::Fun)\n\nreturns the domain that a Fun is defined on\n\n\n\n"
},

{
    "location": "index.html#ApproxFun.coefficients",
    "page": "Home",
    "title": "ApproxFun.coefficients",
    "category": "Function",
    "text": "coefficients(fun,space)\n\nreturns the coefficients of a fun in a possibly different space\n\n\n\n"
},

{
    "location": "index.html#ApproxFun.extrapolate",
    "page": "Home",
    "title": "ApproxFun.extrapolate",
    "category": "Function",
    "text": "extrapolate(fun,x)\n\nreturns an extrapolation of fun from its domain to x.\n\n\n\n"
},

{
    "location": "index.html#ApproxFun.ncoefficients",
    "page": "Home",
    "title": "ApproxFun.ncoefficients",
    "category": "Function",
    "text": "ncoefficients(fun) -> Integer\n\nreturns the number of coefficients of a fun\n\n\n\n"
},

{
    "location": "index.html#ApproxFun.points",
    "page": "Home",
    "title": "ApproxFun.points",
    "category": "Function",
    "text": "points(fun)\n\nreturns a grid of points that the fun can be transformed into values and back\n\n\n\n"
},

{
    "location": "index.html#ApproxFun.space",
    "page": "Home",
    "title": "ApproxFun.space",
    "category": "Function",
    "text": "space(fun)\n\nreturns the space of fun\n\n\n\n"
},

{
    "location": "index.html#Base.values",
    "page": "Home",
    "title": "Base.values",
    "category": "Function",
    "text": "values(fun)\n\nreturns fun evaluated at points(fun)\n\n\n\n"
},

{
    "location": "index.html#Accessing-information-about-a-Fun-1",
    "page": "Home",
    "title": "Accessing information about a Fun",
    "category": "section",
    "text": "domaincoefficientsextrapolatencoefficientspointsspaceApproxFun.values"
},

{
    "location": "index.html#ApproxFun.reverseorientation",
    "page": "Home",
    "title": "ApproxFun.reverseorientation",
    "category": "Function",
    "text": "reverseorientation(fun)\n\nreturn fun on a reversed orientated contour\n\n\n\n"
},

{
    "location": "index.html#ApproxFun.setdomain",
    "page": "Home",
    "title": "ApproxFun.setdomain",
    "category": "Function",
    "text": "setdomain(fun,domain)\n\nreturns fun projected onto domain\n\n\n\n"
},

{
    "location": "index.html#Modify-a-Fun-1",
    "page": "Home",
    "title": "Modify a Fun",
    "category": "section",
    "text": "reverseorientationApproxFun.setdomain"
},

{
    "location": "index.html#ApproxFun.Operator",
    "page": "Home",
    "title": "ApproxFun.Operator",
    "category": "Type",
    "text": "Operator{T}\n\nis an abstract type to represent linear operators between spaces.\n\n\n\n"
},

{
    "location": "index.html#ApproxFun.linsolve",
    "page": "Home",
    "title": "ApproxFun.linsolve",
    "category": "Function",
    "text": "linsolve(A,b;tolerance=tol,maxlength=n)\n\nsolves a linear equation, usually differential equation, where A is an operator or array of operators and b is a Fun or array of funs.  The result u will approximately satisfy A*u = b.\n\n\n\n"
},

{
    "location": "index.html#Base.LinAlg.qrfact-Tuple{ApproxFun.Operator}",
    "page": "Home",
    "title": "Base.LinAlg.qrfact",
    "category": "Method",
    "text": "qrfact(A::Operator)\n\nreturns a cached QR factorization of the Operator A.  The result QR enables solving of linear equations: if u=QR\\b, then u approximately satisfies A*u = b.\n\n\n\n"
},

{
    "location": "index.html#Operators-1",
    "page": "Home",
    "title": "Operators",
    "category": "section",
    "text": "Operatorlinsolveqrfact(::Operator)"
},

]}
