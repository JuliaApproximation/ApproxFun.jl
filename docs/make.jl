using Documenter, ApproxFun, BandedMatrices, DomainSets, LinearAlgebra
using Literate

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

for (example, included) in [
            ("ODE.jl", ["ODE_BVP.jl", "ODE_increaseprec.jl"]),
            ("PDE.jl", ["PDE1.jl"]),
            ("Sampling.jl", String[]),
            ("Periodic.jl", ["Periodic1.jl"]),
            ("Eigenvalue.jl", String[]),
            ("NonlinearBVP.jl", ["NonlinearBVP1.jl"]),
            ("system_of_eqn.jl", String[])
            ]
    filename = joinpath(example_dir, example)
    Literate.markdown(filename, output_dir, documenter=true,
        preprocess = str -> replace_includes(str, included))
end

makedocs(
            doctest = false,
            clean = true,
            format = Documenter.HTML(),
            sitename = "ApproxFun.jl",
            authors = "Sheehan Olver",
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
                    "Frequently Asked Questions" => "faq.md",
                    "Library" => "library.md"
                ]
            )


deploydocs(
    repo   = "github.com/JuliaApproximation/ApproxFun.jl.git",
    push_preview = true,
    )
