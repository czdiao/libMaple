######################################################
#
#   Library related to 2D Filters
#   Also includes some generic functions used for different dilations.
#
#
#   Chenzhe
#   Jun, 2018
#


with(LinearAlgebra):
with(ArrayTools):
with(combinat):
with(PolynomialTools):

VMnD := proc(u, zlist)
    description "Vanishing moments for multivariate Lpoly.":

    local ld, utmp, n, j, tmp:

    n := nops(zlist):
    ld := Vector(n):
    utmp := u:
    for j from 1 to n do
        ld[j] := ldegree(u, zlist[j]):
        utmp := utmp / zlist[j]^ld[j]:
    end do;
    utmp := expand(utmp):

    tmp := [seq(zlist[j] = zlist[j]+1, j = 1..nops(zlist))]:
    utmp := eval(utmp, tmp):
    utmp := expand(utmp):

    return ldegree(utmp):
end proc:

VM2D := proc(u)
    return VMnD(u, [z[1], z[2]]):
end proc:


fsupp2D := proc(a)
    description "fsupp of 2D Lpoly":

    local fs1, fs2, a0:
    a0 := expand(a):
    fs1 := [ldegree(a0, z[1]), degree(a0, z[1])]:
    fs2 := [ldegree(a0, z[2]), degree(a0, z[2])]:

    return fs1, fs2:
end proc:

hc2D := proc(a)
    description "Hermitian Conjugate of 2D filter":
    # only works for real case now.
    
    #return eval( a, {z[1]=1/z[1], z[2]=1/z[2], I = -I}):
    return conjugate( eval(a, [z[1] = 1/conjugate(z[1]), z[2] = 1/conjugate(z[2])]) ):
end proc:

hcMatrix2D := proc(A)
    description "Hermitian Conjugate of Matrix A":
    return Transpose(hc2D(A)):
end proc:

Lpoly2Matrix := proc(f)
    description "Input a 2D filter in z-domain, and output its matrix form, help to debug";

    local a, ld1, ld2, d1, d2, M, fs1, fs2, t, di1, di2, c, m, n, di1m, di2m;

    if evalb(f=0) then
        return Matrix(1,1,0), [-infinity, infinity], [-infinity, infinity]:
    end if;

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

############ Generic procedures ##################

# generic downsample with matrix M0
downsampleM := proc(a0, M0)
    description "Downsample of filter a0 with dilation matrix M0":
    # need to make sure all terms are in 0-coset first. Otherwise would get terms with fractional power.
    # e.g. for quincunx dilation, call 
    #   t:= getAllTerms(a0): t:=select(isCoset0_QCX, t):  a0:= add(t[j], j=1..nops(t)):
    # before this function.

    local a, t, de1, de2, de, c, a0new, Minv:

    Minv := MatrixInverse(M0):
    # downsample
    a := expand(a0):
    a0new := 0:
    if type(a, `+`) then
        for t in [op(a)] do
            de1 := degree(t, z[1]):
            de2 := degree(t, z[2]):
            de := <de1, de2>:
            de := Minv.de:
            c := coeff(coeff(t, z[1], de1), z[2], de2):

            a0new := a0new + c * z[1]^de[1] * z[2]^de[2]:
        od:
    else
        de1 := degree(a, z[1]):
        de2 := degree(a, z[2]):
        de := <de1, de2>:
        de := Minv.de:
        c := coeff(coeff(a, z[1], de1), z[2], de2):
        a0new := a0new + c * z[1]^de[1] * z[2]^de[2]:
    end if;

    return a0new:
end proc:

# generic upsample with matrix M0
upsampleM := proc(a0, M0)
    description "Upsample a 2D filter a0 with dilation matrix M0":

    local a, ld1, d1, ld2, d2, de, di1, di2, c:
    a := 0:

    ld1 := ldegree(a0, z[1]):
    d1 := degree(a0, z[1]):
    ld2 := ldegree(a0, z[2]):
    d2 := degree(a0, z[2]):
    for di1 from ld1 to d1 do
        for di2 from ld2 to d2 do
            c := coeff(coeff(a0, z[1], di1), z[2], di2):
            de := <di1, di2>:
            de := M0.de:
            a := a + c * z[1]^de[1]*z[2]^de[2]:
        end do:
    end do:

    return a:
    
end proc:

# Check whether a matrix is identically zero
is0Matrix := proc(A)
    local Atmp, j, k, m, n:
    m := RowDimension(A):
    n := ColumnDimension(A):
    for j from 1 to m do
        for k from 1 to n do
            if evalb(m=1) then  # row vector
                Atmp := simplify(A[k]):
            elif evalb(n=1) then  # col vector
                Atmp := simplify(A[j]):
            else   # matrix
                Atmp := simplify(A[j,k]):
            end if;

            if not evalb(Atmp = 0) then
                return false:
            end if;
        end do;
    end do;

    return true:
end proc:

######## only for debugging, not used
# "Compute U_1uv as in eq (3.17)"
getU1uv := proc(a0, m)
    description "Compute U_1uv as in eq (3.17)":
    # output Uuv[u, v] with  range: v <= u
    # Uuv[u, v] with nonzero:
    #   u = 0..m, v = 0;
    #   u = m, v = 1..m:

    local Uuv, mu0, u, v, U0, d1, d2, tmp:

    U0 := DividedDiffSp(a0, 2*m):
    Uuv := table():
    for mu0 from 0 to m do
        u := mu0:
        v := 0:
        Uuv[u, v] := U0[mu0, 2*m-mu0] * (-z[1])^u * (-z[2])^(m-u):
    end do;
    for mu0 from m+1 to 2*m do
        u := m:
        v := mu0 - u:
        Uuv[u, v] := U0[mu0, 2*m-mu0] * (-z[1])^u * (-z[2])^(m-u):
    end do;

    # add 0 to other v <= u
    for u from 1 to (m-1) do
        for v from 1 to u do
            Uuv[u, v] := 0:
        end do;
    end do;

    # verify
    tmp := 0:
    for u from 0 to m do
        for v from 0 to u do
            d1 := (z[1]-1)^u * (z[2]-1)^(m-u):
            d1 := hc2D(d1):
            d2 := (z[1]-1)^v * (z[2]-1)^(m-v):
            tmp := tmp + d1*d2*Uuv[u,v]:
        end do;
    end do;
    tmp := simplify(tmp-a0):
    if not evalb(tmp=0) then
        error("Error in getU1uv!"):
    end if;

    return eval(Uuv):
end proc:

# "Compute U_2uv as in eq (3.17)"
getU2uv := proc(a, m)
    description "Compute U_2uv as in eq (3.17)":

    local a1, u2uv, Usp, u, v, tmp, tmp1, tmp2:

    u2uv:= table():
    a1 := CosetShift_QCX(a):
    a1 := hc2D(a1):
    Usp := DividedDiffSp(a1, m):
    for u from 0 to m do
        for v from 0 to m do
            tmp := hcMatrix2D(CosetShift_QCX(Usp[v, m-v])):
            u2uv[u,v] := -(-z[1])^u * (-z[2])^(m-u) /(z[1])^v/(z[2])^(m-v) * Usp[u, m-u] * tmp:
        end do;
    end do;

    # verify
    tmp1 := 0:
    for u from 0 to m do
        for v from 0 to m do
            tmp2 := (z[1]-1)^u * (z[2]-1)^(m-u):
            tmp2 := hc2D(tmp2) * CosetShift_QCX((z[1]-1)^v * (z[2]-1)^(m-v)):
            tmp1 := tmp1 +  tmp2 * u2uv[u,v]:
        end do;
    end do;
    tmp1 := simplify(tmp1 + hc2D(a)*CosetShift_QCX(a) ):
    if not evalb(tmp1 = 0) then
        error("getU2uv not correct!"):
    end if;
    return eval(u2uv):
end proc:
