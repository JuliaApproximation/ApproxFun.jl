#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia -e 'using Pkg; Pkg.add("Documenter")';
    julia -e 'using ApproxFun; cd(joinpath(dirname(pathof(ApproxFun)), "..", "docs")); include("make.jl")';
fi
