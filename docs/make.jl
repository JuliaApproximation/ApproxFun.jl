using Documenter, ApproxFun, BandedMatrices, DomainSets, LinearAlgebra
using Literate

# Generate examples using Literate
# See https://github.com/fredrikekre/Literate.jl/blob/master/docs/make.jl
example_dir = joinpath(@__DIR__, "..", "examples")
output_dir = joinpath(@__DIR__, "src/generated")

for example in ["ODE.jl", "PDE.jl", "Sampling.jl", "Periodic.jl",
		"Eigenvalue.jl", "NonlinearBVP.jl"]
    filename = joinpath(example_dir, example)
    Literate.markdown(filename, output_dir, documenter=true)
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
    repo   = "github.com/JuliaApproximation/ApproxFun.jl.git"
    )
