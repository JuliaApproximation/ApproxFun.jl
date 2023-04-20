using Documenter
using Literate
using ApproxFun, DomainSets, LinearAlgebra, BandedMatrices # these packages are needed for the docs library

# Generate examples using Literate
# See https://github.com/fredrikekre/Literate.jl/blob/master/docs/make.jl
const example_dir = joinpath(@__DIR__, "..", "examples")
const output_dir = joinpath(@__DIR__, "src/generated")

function replace_includes(str, included)
    for ex in included
        content = read(joinpath(example_dir, ex), String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end

function get_included_files(filename)
    v = [l for l in eachline(joinpath(example_dir, filename)) if contains(l, "include")]
    strip.(getindex.(split.(getindex.(split.(v, "include("), 2), ")"), 1), '\"')
end

file_with_includes(filename) = (filename, get_included_files(filename))

for (example, included) in map(file_with_includes,
                ["ODE.jl", "PDE.jl", "Sampling.jl", "Periodic.jl",
                "Eigenvalue.jl", "NonlinearBVP.jl", "system_of_eqn.jl"])
    filename = joinpath(example_dir, example)
    Literate.markdown(filename, output_dir, documenter=true,
        preprocess = str -> replace_includes(str, included))
end

makedocs(
            format = Documenter.HTML(),
            sitename = "ApproxFun.jl",
            authors = "Sheehan Olver, Jishnu Bhattacharya and contributors",
            pages = Any[
                    "Home" => "index.md",
                    "Usage" => Any[
                        "Domains" => "usage/domains.md",
                        "Spaces" => "usage/spaces.md",
                        "Constructors" => "usage/constructors.md",
                        "Operators" => "usage/operators.md",
                        "Linear Equations" => "usage/equations.md"
                    ],
                    "Examples" => [
                        "generated/ODE.md",
                        "generated/system_of_eqn.md",
                        "generated/PDE.md",
                        "generated/Sampling.md",
                        "generated/Periodic.md",
                        "generated/Eigenvalue.md",
                        "generated/NonlinearBVP.md",
                    ],
                    "Internals" => [
                        "internals/multivariate.md",
                        "internals/blocks.md",
                        "internals/spaces.md",
                    ],
                    "Frequently Asked Questions" => "faq.md",
                    "Library" => "library.md"
                ]
            )


deploydocs(
    repo   = "github.com/JuliaApproximation/ApproxFun.jl.git",
    push_preview = true,
    )
