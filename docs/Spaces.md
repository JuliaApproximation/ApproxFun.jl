
The following commands must be overriden to provide features.

# Construction:

points(::FunctionSpace,n::Integer)

Return a grid 

transform(::FunctionSpace,values::Vector)

Transform values on the grid to coefficients.  

# Plotting:

itransform(::FunctionSpace,coefficients::Vector)

Transform coefficients back to values.  Should satisfy that

    transform(space,itransform(space,vals))==vals

# Evaluation:

evaluate(::FunctionSpace,coefficients::Vector,x::Number)

Evaluates the expansion at a point x.  There are no requirements for the default behaviour of x not in the domain.  For example, Chebyshev will give the analytic continuation while Piecewise will return Nothing.


# Conversion:

spacescompatible
canonicalspace 
Base.ones
spaceconversion



# Calculus

integrate




The following commands are standard but optional overrides


domain

    The default is to assume the space has .domain