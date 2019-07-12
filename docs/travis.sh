#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia -e 'using Pkg; Pkg.add("Documenter")';
    julia -e 'using Pkg; cd(Pkg.dir("ApproxFun")); include(joinpath("docs", "make.jl"))';
fi
