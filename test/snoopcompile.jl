# this script generates the SnoopCompile precompile list


using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/approxfun_compiles.csv")

SnoopCompile.@snoop "/tmp/approxfun_compiles.csv" begin
    using ApproxFun
    a=Fun()
    b=Fun(exp)
    a+b
    a*b
    a/b
    a-b

    # following causes SEGFAULT
    # S = Chebyshev()
    # D = Derivative(S)
    # B = dirichlet(S)
    # A = [B;D^2+I]
    # u = A\[1.0,0.0,b]
end

### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts

# IMPORTANT: we must have the module(s) defined for the parcelation
# step, otherwise we will get no precompiles for the ApproxFun module
using ApproxFun

data = SnoopCompile.read("/tmp/approxfun_compiles.csv")


# Blacklist helps fix problems:
# - MIME uses type-parameters with symbols like :image/png, which is
#   not parseable
blacklist = ["MIME"]

# Use these two lines if you want to create precompile functions for
# individual packages
pc, discards = SnoopCompile.parcel(data[end:-1:1,2], blacklist=blacklist)
SnoopCompile.write("/tmp/precompile", pc)
