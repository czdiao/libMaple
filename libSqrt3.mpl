######################################################
#
#   Library related to Sqrt3 Dilation Matrix
#   
#   Generic proc of Sqrt3 Dilation
#
#   M = <<1|-2>,<2|-1>>
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

SR_Sqrt3 := proc(u)
    description "Sum Rule of a 2D lowpass filter on Sqrt3 dilation":
    local sr1, sr2:

    sr1 := simplify(eval(u, [z[1] = z[1]*exp(I*2*Pi/3), z[2] = z[2]*exp(-2*I*2*Pi/3)])):
    sr1 := VM2D(sr1):
    sr2 := simplify(eval(u, [z[1] = z[1]*exp(I*4*Pi/3), z[2] = z[2]*exp(-I*2*Pi/3)])):
    sr2 := VM2D(sr2):
    return min(sr1, sr2):
end proc:

# generic sqrt3
isCoset0_Sqrt3 := proc(t)
    description "Check if a term is in coset [0, 0] for Sqrt3 dilation.":

    local d1, d2:
    d1 := degree(t, z[1]):
    d2 := degree(t, z[2]):

    return evalb(modp(d1+d2, 3)=0):
end proc:

# generic sqrt3
getGammaCoset_Sqrt3 := proc(a, gammaj)
    description "Given a representative gammaj, compute this coset sequence of a, sqrt3 dialtion":
    # gammaj is a vector or list, e.g., [0, 1] or <0, 1>

    local t, ashift, a0, M0:

    # shift gammaj coset to 0-coset
    t := z[1]^gammaj[1] * z[2]^gammaj[2]:
    ashift := expand(a/t):

    # delete other coset, only 0-coset is left
    t:= getAllTerms(ashift, [z[1], z[2]]):
    t:=select(isCoset0_Sqrt3, t):
    ashift:= add(t[j], j=1..nops(t)):

    # downsample
    M0 := <<1|-2>,<2|-1>>:
    a0 := downsampleM(ashift, M0):

    return a0:
end proc:

# generic sqrt3
getCoset_Sqrt3 := proc(a)
    description "Get the coset sequenc of a given 2D Lpoly, with Sqrt3 Dilation matrix":
    # Dilation Matrix M = <<1|-2>,<2|-1>>
    # gamma0 = [0, 0], gamma1 = [1, 0], gamma2 = [1, 1]

    local a0new, a1new, a2new:

    a0new := getGammaCoset_Sqrt3(a, [0, 0]):
    a1new := getGammaCoset_Sqrt3(a, [1, 0]):
    a2new := getGammaCoset_Sqrt3(a, [1, 1]):

    return a0new, a1new, a2new:

end proc:


# currently only works for real input
getN_Sqrt3 := proc(a)
    description "Get matrix N for Sqrt3 Dilation":
    # gamma0 = [0, 0], gamma1 = [1,0], gamma2 = [1,1]
    # dilation matrix M = <<1|-2>,<2|-1>>

    local v1, v2, a0, a1, a2:

    a0, a1, a2 := getCoset_Sqrt3(a):
    v1 := <a0| a1| a2>:
    v2 := hcMatrix2D(v1):
    return v2.v1;
end proc:

# generic sqrt3
upsample_Sqrt3 := proc(a0)
    description "Upsample 2D filter a with Sqrt3 Dilation Matrix":

    local a, M0:

    M0 := <<1|-2>,<2|-1>>:
    a := upsampleM(a0, M0):

    return a:
end proc:

downsample_Sqrt3 := proc(a)
    description "downsample with Sqrt3 dilation":

    local M0:

    M0 := <<1|-2>,<2|-1>>:
    return downsampleM(a, M0):
end proc:


# generic sqrt3
getFilter_Sqrt3 := proc(a0, a1, a2)
    description "Get back to filter from coset sequences.":
    # gamma0 = [0, 0], gamma1 = [1,0], gamma2 = [1,1]

    local a:
    a := 0:

    a := a + upsample_Sqrt3~(a0):
    a := a + upsample_Sqrt3~(a1)*z[1]:
    a := a + upsample_Sqrt3~(a2)*z[1]*z[2]:

    return simplify(a):
end proc:


# Defined as in Lemma 7
getEb0_Sqrt3 := proc(b)
    description "Generate matrix Eb0 for filter b.":
    # Sqrt3 Dilation, other dilations could be modified easily.

    local dm, Eb, gList, gamma0, gamma1, gamma2, j, k:

    dm := 3:
    Eb := Matrix(dm, dm):

    gamma0:= <0, 0>:
    gamma1:= <1, 0>:
    gamma2:= <1, 1>:
    gList := [gamma0, gamma1, gamma2]:

    for j from 1 to dm do
        for k from 1 to dm do
            Eb[j,k] := getGammaCoset_Sqrt3(b, gList[k]-gList[j]):
        end do;
    end do;

    return Eb:
end proc:


# Defined as in Thm 8
getEmu_Sqrt3 := proc(mu)
    description "Generate Eb0 for divided difference delta":
    # mu is a list/vector with 2 component [mu1, mu2] or <mu1, mu2>

    local b, mu1, mu2:

    mu1:= mu[1]:
    mu2:= mu[2]:
    b := (z[1]-1)^mu1 * (z[2]-1)^mu2:
    b := expand(b):

    return getEb0_Sqrt3(b):
end proc:

CosetShiftMatrix_Sqrt3 := proc(fb)
    description "Coset Shift Matrix with Sqrt3 Dilation, used to verify PR":
    # input a column vector with n filters, output an n x 3 matrix.
    local f1, f2, FB:

    f1 := simplify(eval(fb, [z[1] = z[1]*exp(I*2*Pi/3), z[2] = z[2]*exp(2*I*Pi/3)])):
    f2 := simplify(eval(fb, [z[1] = z[1]*exp(I*4*Pi/3), z[2] = z[2]*exp(4*I*Pi/3)])):

    FB := <fb| f1| f2>:
    return FB:
end proc:

