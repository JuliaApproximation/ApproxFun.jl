#  Contains code that is based in part on Chebfun v5's chebfun/standardChop.m,
# which is distributed with the following license:

# Copyright (c) 2017, The Chancellor, Masters and Scholars of the University
# of Oxford, and the Chebfun Developers. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Oxford nor the names of its
#       contributors may be used to endorse or promote products derived from
#       this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Jared Aurentz and Nick Trefethen, July 2015.
#
# Copyright 2017 by The University of Oxford and The Chebfun Developers.
# See http://www.chebfun.org/ for Chebfun information.

function standardchoplength(coeffs, tol)
    # Check magnitude of TOL:
    if tol ≥ 1
        throw(ArgumentError("tolerance must be less than 1"))
    end

    # Make sure COEFFS has length at least 17:
    n = length(coeffs)

    if  n < 17
        return n
    end

    # Step 1: Convert COEFFS to a new monotonically nonincreasing
    #         vector ENVELOPE normalized to begin with the value 1.

    b = abs.(coeffs)
    m = b[end]*ones(n)
    for j = n-1:-1:1
        m[j] = max(b[j], m[j+1]);
    end
    if  m[1] == 0
        return 1
    end
    envelope = m/m[1];

    plateauPoint = 0

    for j = 2:n
        j2 = round(Int,1.25*j + 5);
        if  j2 > n
            # there is no plateau: exit
            return n
        end
        e1 = envelope[j]
        e2 = envelope[j2]
        r = 3*(1 - log(e1)/log(tol))
        plateau = (e1 == 0) || (e2/e1 > r)
        if plateau
            # a plateau has been found: go to Step 3
            plateauPoint = j - 1
            break
        end
    end

    # Step 3: fix CUTOFF at a point where ENVELOPE, plus a linear function
    # included to bias the result towards the left end, is minimal.
    #

    if  plateauPoint ≠ 0 && envelope[plateauPoint] == 0
        return plateauPoint
    end

    j3 = sum(envelope .≥ tol^(7/6))
    if j3 < j2
        j2 = j3 + 1
        envelope[j2] = tol^(7/6)
    end
    cc = log10.(envelope[1:j2])
    cc .+= linspace(0, (-1/3)*log10(tol), j2)
    d = indmin(cc)
    return max(d - 1, 1)
end
