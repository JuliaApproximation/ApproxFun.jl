#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia -e 'Pkg.add("ApproxFun")';
    julia -e 'Pkg.add("Documenter")';
    julia -e 'cd(Pkg.dir("ApproxFun")); include(joinpath("docs", "make.jl"))';
fi
