######################################################
#
#   Library related to I2 Dilation Matrix
#   
#   Generic proc of 2I2 Dilation
#
#   M = <<2|0>,<0|2>>
#
#   Call  libFilter2D.mpl first by
#       read "libFilter2D.mpl";
#
#   Chenzhe
#   Jan, 2018
#

#with(LinearAlgebra):
#with(ArrayTools):
#with(combinat):
#with(PolynomialTools):

read "libFilter2D.mpl":

SR_I2 := proc(u)
    description "Sum Rule of a 2D lowpass filter on I2 dilation":
    local sr1, sr2, sr3:

    sr1 := simplify(eval(u, [z[1] = -z[1]])):
    sr1 := VM2D(sr1):
    sr2 := simplify(eval(u, [z[2] = -z[2]])):
    sr2 := VM2D(sr2):
    sr3 := simplify(eval(u, [z[1] = -z[1], z[2] = -z[2]])):
    sr3 := VM2D(sr3):
    return min(sr1, sr2, sr3):
end proc:

# generic I2
isCoset0_I2 := proc(t)
    description "Check if a term is in coset [0, 0] for I2 dilation.":

    local d1, d2, e1, e2:
    d1 := degree(t, z[1]):
    d2 := degree(t, z[2]):
    e1 := evalb(modp(d1, 2)=0):
    e2 := evalb(modp(d2, 2)=0):

    return (e1 and e2):
end proc:

# generic I2
getGammaCoset_I2 := proc(a, gammaj)
    description "Given a representative gammaj, compute this coset sequence of a, I2 dialtion":
    # gammaj is a vector or list, e.g., [0, 1] or <0, 1>

    local t, ashift, a0, M0:

    # shift gammaj coset to 0-coset
    t := z[1]^gammaj[1] * z[2]^gammaj[2]:
    ashift := expand(a/t):

    # delete other coset, only 0-coset is left
    t:= getAllTerms(ashift, [z[1], z[2]]):
    t:=select(isCoset0_I2, t):
    ashift:= add(t[j], j=1..nops(t)):

    # downsample
    M0 := <<2|0>,<0|2>>:
    a0 := downsampleM(ashift, M0):

    return a0:
end proc:

# generic I2
getCoset_I2 := proc(a)
    description "Get the coset sequenc of a given 2D Lpoly, with I2 Dilation matrix":
    # Dilation Matrix M = <<2|0>,<0|2>>
    # gamma0 = [0, 0], gamma1 = [1, 0], gamma2 = [0, 1], gamma3 = [1, 1]

    local a0new, a1new, a2new, a3new:

    a0new := getGammaCoset_I2(a, [0, 0]):
    a1new := getGammaCoset_I2(a, [1, 0]):
    a2new := getGammaCoset_I2(a, [0, 1]):
    a3new := getGammaCoset_I2(a, [1, 1]):

    return a0new, a1new, a2new, a3new:

end proc:


# currently only works for real input
getN_I2 := proc(a)
    description "Get matrix N for I2 Dilation":
    # gamma0 = [0, 0], gamma1 = [1,0], gamma2 = [0, 1], gamma3 = [1, 1]
    # Dilation Matrix M = <<2|0>,<0|2>>

    local v1, v2, a0, a1, a2, a3:

    a0, a1, a2, a3 := getCoset_I2(a):
    v1 := <a0| a1| a2| a3>:
    v2 := hcMatrix2D(v1):
    return v2.v1;
end proc:

# generic I2
upsample_I2 := proc(a0)
    description "Upsample 2D filter a with I2 Dilation Matrix":

    local a, M0:

    M0 := <<2|0>,<0|2>>:
    a := upsampleM(a0, M0):

    return a:
end proc:

downsample_I2 := proc(a)
    description "downsample with I2 dilation":

    local M0:

    M0 := <<2|0>,<0|2>>:
    return downsampleM(a, M0):
end proc:


# generic I2
getFilter_I2 := proc(a0, a1, a2, a3)
    description "Get back to filter from coset sequences.":
    # gamma0 = [0, 0], gamma1 = [1,0], gamma2 = [0, 1], gamma3 = [1, 1]

    local a:
    a := 0:

    a := a + upsample_I2~(a0):
    a := a + upsample_I2~(a1)*z[1]:
    a := a + upsample_I2~(a2)*z[2]:
    a := a + upsample_I2~(a3)*z[1]*z[2]:

    return simplify(a):
end proc:


# Defined as in Lemma 7
getEb0_I2 := proc(b)
    description "Generate matrix Eb0 for filter b.":
    # I2 Dilation, other dilations could be modified easily.

    local dm, Eb, gList, gamma0, gamma1, gamma2, gamma3, j, k:

    dm := 4:
    Eb := Matrix(dm, dm):

    gamma0:= <0, 0>:
    gamma1:= <1, 0>:
    gamma2:= <0, 1>:
    gamma3:= <1, 1>:
    gList := [gamma0, gamma1, gamma2, gamma3]:

    for j from 1 to dm do
        for k from 1 to dm do
            Eb[j,k] := getGammaCoset_I2(b, gList[k]-gList[j]):
        end do;
    end do;

    return Eb:
end proc:


# Defined as in Thm 8
getEmu_I2 := proc(mu)
    description "Generate Eb0 for divided difference delta":
    # mu is a list/vector with 2 component [mu1, mu2] or <mu1, mu2>

    local b, mu1, mu2:

    mu1:= mu[1]:
    mu2:= mu[2]:
    b := (z[1]-1)^mu1 * (z[2]-1)^mu2:
    b := expand(b):

    return getEb0_I2(b):
end proc:

CosetShiftMatrix_I2 := proc(fb)
    description "Coset Shift Matrix with I2 Dilation, used to verify PR":
    # input a column vector with n filters, output an n x 4 matrix.
    local f1, f2, f3, FB:

    f1 := simplify(eval(fb, [z[1] = -z[1]])):
    f2 := simplify(eval(fb, [z[2] = -z[2]])):
    f3 := simplify(eval(fb, [z[1] = -z[1], z[2] = -z[2]])):

    FB := <fb| f1| f2| f3>:
    return FB:
end proc:

SumRow := proc(A)
    local m, n, s, j, k:
    m, n := Dimension(A):
    s := [seq(add(A[j, k], k = 1..n), j = 1..m)]:
    return s:
end proc:

SumCol := proc(A)
    local m, n, s, j, k:
    m, n := Dimension(A):
    s := [seq(add(A[j, k], j = 1..m), k = 1..n)]:

    return s:
end proc:

SumDiag := proc(A)
    local m, n, s1, s2, s, j, k:
    m, n := Dimension(A):
    s1 := [seq(add(A[j, k-j+1], j = 1..min(m, k)), k = 1..n)]:
    s2 := [seq(add(A[j+k-1, n-k+1], k = 1..min(n, m+1-j)), j = 2..m)]:

    s := [op(s1), op(s2)]:

    return s:
end proc: