using Documenter, ApproxFun, BandedMatrices, DomainSets, LinearAlgebra

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
					"Frequently Asked Questions" => "faq.md",
					"Library" => "library.md"
				]
			)


deploydocs(
    repo   = "github.com/JuliaApproximation/ApproxFun.jl.git",
    push_preview = true,
    )
