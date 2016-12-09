# Blocks

Every space is divided into blocks, which are used to indicate, for example,
polynomial degree.  The prototypical examples have trivial blocks,
for example `Taylor()` and `Chebyshev()` have blocks of length 1, as
each coefficient corresponds to a higher polynomial degree.   

Usually, non-trivial block lengths arise from modifications of the
spaces with trivial blocks.  For example, `Chebyshev(0..1) ∪ Chebyshev(2..3)`
has blocks of length 2, as the blocks of each component space are grouped
together to form a single block.  Another important example is
`Chebyshev() ⊗ Chebyshev()`, the tensor product space.  There are `d` polynomials
of degree `d`, thus the blocks of a tensor product space grow: that is, the first
block has length 1, then 2, and so on.

`blocklengths(::Space)` gives an iterator that encodes the lengths of the blocks.
For trivial blocks, this will return `ApproxFun.repeated(true)`.  For
`Chebyshev(0..1) ∪ Chebyshev(2..3)` it returns `ApproxFun.repeated(2)`.
For `Chebyshev() ⊗ Chebyshev()` it returns `1:∞`, which is equivalent to
`ApproxFun.UnitCount(1)`.  (Note that the `Base` iterators have been recreated in
ApproxFun to add missing features.)
