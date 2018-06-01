######################################################
#
#   Library related to 2D Filters
#   With some functions for quincunx dilation.
#
#
#   Chenzhe
#   Jun, 2018
#


with(LinearAlgebra):
with(ArrayTools):
with(combinat):
with(PolynomialTools):

fsupp2D := proc(a)
    description "fsupp of 2D Lpoly":

    local fs1, fs2, a0:
    a0 := expand(a):
    fs1 := [ldegree(a0, z[1]), degree(a0, z[1])]:
    fs2 := [ldegree(a0, z[2]), degree(a0, z[2])]:

    return fs1, fs2:
end proc:

CosetShift := proc(a)
    description "Shift coset with Quincunx Dilation":
    return eval(a, {z[1]=-z[1], z[2]=-z[2]}):
end proc:

# currently only works for real input
hc2D := proc(a)
    description "Hermitian Conjugate of 2D filter":
    # only works for real case now.
    
    return eval(a, {z[1]=1/z[1], z[2]=1/z[2]}):
end proc:

Lpoly2Matrix := proc(f)
    description "Input a 2D filter in z-domain, and output its matrix form, help to debug";

    local a, ld1, ld2, d1, d2, M, fs1, fs2, t, di1, di2, c, m, n, di1m, di2m;

    a := expand(f):

    ld1 := ldegree(a, z[1]):
    d1 := degree(a, z[1]):
    ld2 := ldegree(a, z[2]):
    d2 := degree(a, z[2]):
    fs1 := [ld1, d1];
    fs2 := [ld2, d2];

    m := d1 - ld1 +1:
    n := d2 - ld2 +1:
    M := Matrix(n, m):

    for di1 from ld1 to d1 do
        for di2 from ld2 to d2 do
            c := coeff(coeff(a, z[1], di1), z[2], di2):
	        di1m := di1 - ld1 +1:
	        di2m := di2 - ld2 +1:
	        M[n-di2m+1, di1m] := c:
        end do:
    end do:

    return M, fs1, fs2;
end proc:

Matrix2Lpoly := proc(M, fs1, fs2)
    description "Input a 2D filter in matrix form, M is the matrix, fs1, fs2 are the fsupp of z[1] and z[2]":
    # fs1 and fs2 are lists with 2 elements. Matrix starts from lower-left corner, col for z[1], row for z[2].

    local m, n, a, t, c, i1, i2:

    a:= 0;
    n, m := Dimension(M):
    for i1 from 1 to m do
        for i2 from 1 to n do
            t := M[n-i2+1, i1] * z[1]^(fs1[1] + i1 - 1) * z[2]^(fs2[1] + i2 - 1):
            a := a + t:
        od:
    od:
    return a;
end proc:

# currently only works for real input
getM_QCX := proc(a)
    description "Get matrix M for Quincunx Dilation":

    local a1, v, v1:

    a1:= eval(a, {z[1] = -z[1], z[2] = -z[2]}):
    v := <a, a1>:
    v1 := Transpose(eval(v, {z[1]=1/z[1], z[2]=1/z[2]})):

    return v.v1;
end proc:


getCoset_QCX := proc(a)
    description "Get the coset sequenc of a given 2D Lpoly, with Quincunx Dilation matrix":
    # Dilation Matrix M = [1, 1; 1, -1]
    # gamma0 = [0, 0], gamma1 = [1, 0]

    local an, a0, a1, a0new, a1new, M, Minv, de, de1, de2, c, t:
    # find 2 coset sequences
    an := eval(a, {z[1]=-z[1], z[2]=-z[2]}):
    a0 := simplify(a+an)/2:
    a0 := expand(a0):
    a1 := (a - a0)/z[1]:
    a1 := expand(a1):

    a0new := 0:
    a1new := 0:

    # Quincunx Dilation
    M := Matrix([[1, 1], [1, -1]]):
    Minv := MatrixInverse(M):

    # downsample
    if type(a0, `+`) then
        for t in [op(a0)] do
            de1 := degree(t, z[1]):
            de2 := degree(t, z[2]):
            de := <de1, de2>:
            de := Minv.de:
            c := coeff(coeff(t, z[1], de1), z[2], de2):

            a0new := a0new + c * z[1]^de[1] * z[2]^de[2]:
        od:
    else
        de1 := degree(a0, z[1]):
        de2 := degree(a0, z[2]):
        de := <de1, de2>:
        de := Minv.de:
        c := coeff(coeff(a0, z[1], de1), z[2], de2):
        a0new := a0new + c * z[1]^de[1] * z[2]^de[2]:
    end if;
    
    if type(a1, `+`) then
        for t in [op(a1)] do
            de1 := degree(t, z[1]):
            de2 := degree(t, z[2]):
            de := <de1, de2>:
            de := Minv.de:
            c := coeff(coeff(t, z[1], de1), z[2], de2):

            a1new := a1new + c * z[1]^de[1] * z[2]^de[2]:
        od:
    else
        de1 := degree(a1, z[1]):
        de2 := degree(a1, z[2]):
        de := <de1, de2>:
        de := Minv.de:
        c := coeff(coeff(a1, z[1], de1), z[2], de2):
        a1new := a1new + c * z[1]^de[1] * z[2]^de[2]:
    end if;

    return a0new, a1new:

end proc:


# currently only works for real input
getN_QCX := proc(a)
    description "Get matrix N for Quincunx Dilation":
    # gamma0 = [0, 0], gamma1 = [1,0]
    # dilation matrix M = [1, 1; 1, -1]

    local v1, v2, a0, a1:

    a0, a1 := getCoset_QCX(a):
    v1 := <a0, a1>:
    v2 := Transpose(eval(v1, {z[1]=1/z[1], z[2]=1/z[2]})):
    return v1.v2;
end proc:

getFilter_QCX := proc(a0, a1)
    description "Get back to filter from coset sequences.":

    local a, M:
    a := 0:

    a := a + upsample_QCX(a0):
    a := a + upsample_QCX(a1)*z[1]:

    return simplify(a):
end proc:

upsample_QCX := proc(a0)
    description "Upsample 2D filter a with Quincunx Dilation Matrix":

    local a, M, ld1, d1, ld2, d2, de, di1, di2, c:
    a := 0:
    M := Matrix([[1, 1], [1, -1]]):

    ld1 := ldegree(a0, z[1]):
    d1 := degree(a0, z[1]):
    ld2 := ldegree(a0, z[2]):
    d2 := degree(a0, z[2]):
    for di1 from ld1 to d1 do
        for di2 from ld2 to d2 do
            c := coeff(coeff(a0, z[1], di1), z[2], di2):
            de := <di1, di2>:
            de := M.de:
            a := a + c * z[1]^de[1]*z[2]^de[2]:
        end do:
    end do:

    return a:
    
end proc:
