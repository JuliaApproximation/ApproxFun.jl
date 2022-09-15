# Random number sampling using [Olver & Townsend 2013](https://arxiv.org/abs/1307.1223).
# The following code samples 10,000 from a PDF given as the absolute value of the sine function on [-5,5]:
using ApproxFun

f = abs(Fun(sin,-5..5))
x = ApproxFun.sample(f, 10000);
