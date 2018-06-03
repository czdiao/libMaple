######################################################
#
#   Library related to Quincunx Dilation Matrix
#   
#   Generic proc of Quincunx Dilation
#
#   M = Matrix([[1, 1], [1, -1]]):
#
#   Call  libFilter2D.mpl first by
#       read "libFilter2D.mpl";
#
#   Chenzhe
#   Jun, 2018
#

read "libFilter2D.mpl":

############ Filters with Quincunx Dilation Matrix #############

SR_QCX := proc(u)
    description "Sum Rule of a 2D lowpass filter on quincunx dilation":

    return VM2D(eval(u, [z[1] = -z[1], z[2]= -z[2]])):
end proc:

CosetShift_QCX := proc(a)
    description "Shift coset with Quincunx Dilation":
    return eval(a, {z[1]=-z[1], z[2]=-z[2]}):
end proc:


# generic qcx
isCoset0_QCX := proc(t)
    description "Check if a term is in coset [0, 0] for quincunx dilation.":

    local d1, d2:
    d1 := degree(t, z[1]):
    d2 := degree(t, z[2]):

    return evalb(modp(d1+d2, 2)=0):
end proc:

# generic qcx
getGammaCoset_QCX := proc(a, gammaj)
    description "Given a representative gammaj, compute this coset sequence of a, quincunx dialtion":
    # gammaj is a vector or list, e.g., [0, 1] or <0, 1>

    local t, ashift, a0, M0:

    # shift gammaj coset to 0-coset
    t := z[1]^gammaj[1] * z[2]^gammaj[2]:
    ashift := expand(a/t):

    # delete other coset, only 0-coset is left
    t:= getAllTerms(ashift, [z[1], z[2]]):
    t:=select(isCoset0_QCX, t):
    ashift:= add(t[j], j=1..nops(t)):

    # downsample
    M0 := Matrix([[1, 1], [1, -1]]):
    a0 := downsampleM(ashift, M0):

    return a0:
end proc:

# generic qcx
getCoset_QCX := proc(a)
    description "Get the coset sequenc of a given 2D Lpoly, with Quincunx Dilation matrix":
    # Dilation Matrix M = [1, 1; 1, -1]
    # gamma0 = [0, 0], gamma1 = [1, 0]

    local a0new, a1new:

    a0new := getGammaCoset_QCX(a, [0, 0]):
    a1new := getGammaCoset_QCX(a, [1, 0]):

    return a0new, a1new:

end proc:

# generic qcx
upsample_QCX := proc(a0)
    description "Upsample 2D filter a with Quincunx Dilation Matrix":

    local a, M0:

    M0 := Matrix([[1, 1], [1, -1]]):
    a := upsampleM(a0, M0):

    return a:
end proc:

downsample_QCX := proc(a)
    description "downsample with Quincunx dilation":

    local M0:

    M0 := Matrix([[1, 1], [1, -1]]):
    return downsampleM(a, M0):
end proc:

# generic qcx
getFilter_QCX := proc(a0, a1)
    description "Get back to filter from coset sequences.":

    local a:
    a := 0:

    a := a + upsample_QCX~(a0):
    a := a + upsample_QCX~(a1)*z[1]:

    return simplify(a):
end proc:


######## Matrix with Quincunx Dilation ####################

# currently only works for real input
getM_QCX := proc(a)
    description "Get matrix M for Quincunx Dilation":

    local a1, v, v1:

    a1:= eval(a, {z[1] = -z[1], z[2] = -z[2]}):
    v := <a| a1>:
    v1 := Transpose(eval(v, {z[1]=1/z[1], z[2]=1/z[2]})):

    return v1.v;
end proc:

# currently only works for real input
getN_QCX := proc(a)
    description "Get matrix N for Quincunx Dilation":
    # gamma0 = [0, 0], gamma1 = [1,0]
    # dilation matrix M = [1, 1; 1, -1]

    local v1, v2, a0, a1:

    a0, a1 := getCoset_QCX(a):
    v1 := <a0| a1>:
    v2 := Transpose(eval(v1, {z[1]=1/z[1], z[2]=1/z[2]})):
    return v2.v1;
end proc:


# Defined as in Lemma 7
getEb0_QCX := proc(b)
    description "Generate matrix Eb0 for filter b.":
    # Quincunx Dilation, other dilations could be modified easily.

    local dm, Eb, gList, gamma0, gamma1, j, k:

    dm := 2:
    Eb := Matrix(dm, dm):

    gamma0:= <0, 0>:
    gamma1:= <1, 0>:
    gList := [gamma0, gamma1]:

    for j from 1 to dm do
        for k from 1 to dm do
            Eb[j,k] := getGammaCoset_QCX(b, gList[k]-gList[j]):
        end do;
    end do;

    return Eb:
end proc:

# Defined as in Lemma 7
getEb1_QCX := proc(b)
    description "Generate matrix Eb1 for filter b.":
    # Quincunx Dilation, other dilations could be modified easily.

    local dm, Eb, j:

    dm := 2:
    Eb := getEb0_QCX(b):
    for j from 1 to dm do
        Eb[j,2] := -Eb[j,2]:
    end do;

    return Eb:
end proc:


# Defined as in Thm 8
getEmu_QCX := proc(mu)
    description "Generate Eb0 for divided difference delta":
    # mu is a list/vector with 2 component [mu1, mu2] or <mu1, mu2>

    local b, mu1, mu2:

    mu1:= mu[1]:
    mu2:= mu[2]:
    b := (z[1]-1)^mu1 * (z[2]-1)^mu2:
    b := expand(b):

    return getEb0_QCX(b):
end proc:
