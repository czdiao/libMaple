###########################################################
#
#   Library related to Matrix Split and construction of 
#   framelet filter banks.
#
#   This library relies on:
#       1. libLpoly.mpl
#       2. libFejerRiesz.mpl
#
#
#   Chenzhe
#   Oct, 2016
#

with(LinearAlgebra):
with(ArrayTools):
with(combinat):
with(PolynomialTools):



isVectorOrMatrix := proc(A)
    not type(A, scalar):
end:

negVec := proc(a)
    description "Generate matrix [a(z); a(-z)]";
    Matrix([ [a], [eval(a, z = -z)] ]):
end proc:


m2prod := proc(M1, M2)
    description " M1*hermConj(M2) , M1 and M2 can be either matrices of scalars";
    local B:

    if ( isVectorOrMatrix(M2) ) then
	B := simplify( M1.hermConj(M2) ):
    else ### assume total_elements == 1, scalar case
	B := M1*hermconj(M2):
    end if:
    B:
end:

m2sqr := proc(M)
    description " M*hermConj(M) ";
    m2prod(M, M):
end:

neg2prod := proc(a1, a2)
    m2prod( negVec(a1), negVec(a2) ):
end:

neg2sqr := proc(a)
    m2sqr(negVec(a)):
end:

getMTheta := proc(a1, a2, theta)
## get 2x2 matrix theta, input is Laurant polys a1, a2, theta
## output is a 2x2 matrix theta
    local a1_neg, a1_hermite, a1_neg_herm, a2_neg, a2_hermite,
    a2_neg_herm, theta_neg, theta2, M:

    theta_neg := eval(theta, z = -z):
    theta2 := eval(theta, z = z^2):
    M := Matrix([ [theta, 0], [0, theta_neg] ])
    -theta2*neg2prod(a1, a2):
    M := simplify(M):
end proc:

seq2f := proc(a)
# convert a seqence to a Laurent poly, ldeg is the lowest degree.

    local f, n, ldeg, hdeg:

    ldeg := a[1]: hdeg := a[2]: f := 0:
    for n from ldeg by 1 to hdeg do
	f := f + a[n - ldeg + 3]*z^n:
    od:
    f:
end proc:

f2seq := proc(f)
    local vec, i, ldeg, hdeg:

    ldeg := ldegree(f, z): hdeg := degree(f, z):
    if (f = 0) then
	ldeg := 0: hdeg := 0:
    end if:
    vec := Vector[row](hdeg - ldeg + 1 + 2):
    vec[1] := ldeg: vec[2] := hdeg:
    for i from 0 to hdeg - ldeg do
	vec[3 + i] := coeff(f, z, ldeg + i):
    od:
    vec:
end proc:


f2seq_reverse := proc(f)
    local vec, i, ldeg, hdeg:

    ldeg := ldegree(f, z): hdeg := degree(f, z):
    if (f = 0) then
	ldeg := 0: hdeg := 0:
    end if:
    vec := Vector[row](hdeg - ldeg + 1 + 2):
    vec[1] := ldeg: vec[2] := hdeg:
    for i from 0 to hdeg - ldeg do
	vec[3 + i] := coeff(f, z, hdeg - i):
    od:
    vec:
end proc:



getM2 := proc(a1, a2, theta, n1, n2)
    ### get simplified MTheta, i.e., factored-outed M.

    local M, d, d2_hermite, M2:

    M := getMTheta(a1, a2, theta):
    d := Matrix( [ [1/(1 - z)^n1, 0], [0, 1/(1 + z)^n1] ] ):
    d2_hermite := Matrix( [ [1/(1 - 1/z)^n2, 0], [0, 1/(1 + 1/z)^n2] ] ):
    M2 := simplify(d.M.d2_hermite):
end proc:

getP1 := proc()
    ## get the 2x2 matrix [ [1, 1], [z, -z] ]
    local M:

    M := Matrix( [ [1, 1], [z, -z] ]):
end proc:

getP1Herm := proc()
    ## get the 2x2 matrix [ [1, 1/z], [1, -1/z] ]
    local M:

    M := Matrix( [ [1, 1/z], [1, -1/z] ] ):
end proc:

P1Prod := proc(M)
    ## get the product of P1.M.P1*
    local M2:

    M2 := simplify( getP1().M.getP1Herm() ):
    return M2:
end proc:

getP2 := proc()
## get the 2x2 matrix [ [1 + z, 1 - z], [1 - z, 1 + z] ]
	local M:

	M := Matrix( [ [1 + z, 1 - z], [1 - z, 1 + z] ]):
end proc:

getP2Herm := proc()
## get the 2x2 matrix [ [1 + 1/z, 1 - 1/z], [1 - 1/z, 1 + 1/z] ]
	local M:

	M := Matrix( [ [1 + 1/z, 1 - 1/z], [1 - 1/z , 1 + 1/z] ] ):
end proc:

P2Prod := proc(M)
## get the product of P2.M.P2*
	local M2:

	M2 := simplify( getP2().M.getP2Herm() ):
end proc:



############# Generate Dual lowpass filter for biorthogonal filter bank #########

DualLowpass := proc(a)
    description "Find the dual lowpass filter";
    local m, ashift, ainv, u, v, adual:

    m := ldegree(a, z):
    ashift:= simplify(a/z^m):
    ainv:= eval(ashift, z=-z):
    gcdex(ashift, ainv, 1, z, 'u', 'v'):
    u := hc(collect(u, z)):
    v := hc(collect(v, z)):
    v := eval(v, z=-z):
    adual:= (u + v)*z^m/2:

    return simplify(adual):

end proc:

DualLowpassComplexSym := proc(a)
    description "For a with complex sym, find dual lowpass with same sym type";

    local SS, adual:

    adual:= DualLowpass(a):
    if isComplexSym(a) then
        SS:= Lquo(a, hc(a)):
        adual := (adual + SS*hc(adual))/2:
    else
        error("Input lowpass is not complex sym"):
    end if:

    return simplify(adual):

end proc:

CheckDual := proc(a, aa)
    description "Check if a and aa form pair of dual lowpass filters, output should be 1":
    local ainv, aainv, tmp:

    ainv:= eval(a, z=-z):
    aainv:= eval(aa, z=-z):

    tmp:= a*hc(aa) + ainv*hc(aainv):
    return simplify(tmp):
end proc:


BsplineLowpass := proc(m)
    description "Generate lowpass of B-spline filter bank of sum rule m":
    local a:

    a:=collect((1+z)^m/2^m,z);
    return a:

end proc:

# Construction of a_{2m, 2n}
Pmn := proc(m, n, t)
    local f, x:
    f:= MTM[taylor](1/(1-x)^m, n):
    f:= eval(f, x=t):

    return f:
end proc:

LowpassSrLpm := proc(M, N)
    description "Construction of lowpass filters with M sum rule, N lpm, N has to be even";
    local p1, p2, p, m, n:

    if modp(N, 2)<>0 then
        error("Cannot construct odd lpm"):
    end if:

    n:= N/2:
    p1:= (z+1/z+2)/4:   # cos^2(xi/2)
    p2:= (2-z-1/z)/4:   # sin^2(xi/2)

    if modp(M, 2)=0 then
        m:= M/2:
        p:= p1^m * Pmn(m, n, p2):

        return collect(p,z):
    else
        m:= (M+1)/2:
        p:= (1+z) * p1^(m-1) * Pmn(m-1/2, n, p2)/2:

        return collect(p,z):
    end if:

end proc:

LowpassSrLpmEven := proc(m, n)
    description "Lowpass with 2m sum rule, 2n lpm";
    local p1, p2, p:

    p1:= (z+1/z+2)/4:   # cos^2(xi/2)
    p2:= (2-z-1/z)/4:   # sin^2(xi/2)

    p:= p1^m * Pmn(m, n, p2):

    return collect(p,z):
end proc:

LowpassSrLpmOdd := proc(m, n)
    description "Lowpass with 2m-1 sum rule, 2n lpm";
    local p1, p2, p:

    p1:= (z+1/z+2)/4:   # cos^2(xi/2)
    p2:= (2-z-1/z)/4:   # sin^2(xi/2)

    p:= (1+z) * p1^(m-1) * Pmn(m-1/2, n, p2)/2:

    return collect(p,z):
end proc:


################ Compute minimum nb ###############

ZLpoly := proc(p, z, z0)
    description "Compute the multiplicity of zeroes at z0, p(z) is a Lpoly":

    local k, len:
    len:= lenLpoly(p):
    for k from 1 to len+1 do
        if not Ldivide(p, (z-z0)^k) then
            break:
        end if:
    od:

    return k-1:
end proc:

SumRule := proc(a)
    description "Compute the sum rule of a":
    local k, len:

    len:= lenLpoly(a):
    for k from 1 to len+1 do
        if not Ldivide(a, (1+z)^k) then
            break:
        end if:
    od:

    return k-1:
end proc:

VanMom := proc(b)
    description "Compute the vanishing moment of b";
    local k, len:

    len:= lenLpoly(b):
    for k from 1 to len+1 do
        if not Ldivide(b, (z-1)^k) then
            break:
        end if:
    od:

    return k-1:
end proc:

MinVM := proc(a, Theta)
    description "Compute minimum possible vanishing moment nb";

    local sr, vm, aa, t:

    sr:= SumRule(a):
    aa:= simplify(Theta - eval(Theta, z=z^2)*a*hc(a)):
    vm:= VanMom(aa):

    t:= min(sr, vm/2):

    return t:


end proc:


#### TODO list:
#3. real/complex SOS/DOS


############### For Indefinite Matrix Split #####################
getN2 := proc(M2)
    description "Given M, calculate N";
    local N2, T, T_herm, T_inverse;

    T:= Matrix([[1, z], [1, -z]]):
    T_herm:=hermConj(T):
    T_inverse:= 1/2 * T_herm:
    N2:= simplify(T_inverse.M2.T/2):

    N2:= collect(N2, z):
    N2:= eval(N2, z=z^(1/2)):

    return N2:

end proc:

MatrixGCD := proc(N)
    description "Input 2x2 laurent matrix N, calculate gcd of 4 elements";
    local deg, ldeg, p, tmp, ii, jj;

    p:=0:
    for ii from 1 to 2 do
		for jj from 1 to 2 do
			tmp:=N[ii,jj]:
			ldeg:=-ldegree(tmp,z):
			tmp:=collect(tmp*z^ldeg, z):

			p:=gcd(p, tmp):
		od:
    od:

    ## Still need to normalize it to be real on T
	if modp(lenLpoly(p),2)=0 and isComplexSym(p) then
		# now p could be normalized to be real on T
		p:=NormRealonT(p):
	else
		print("gcd of all elements of the matrix is not real on T"):
	fi:
	return p:
end proc:

MatrixRowGCD := proc(N)
    description "Return a list of 2 Matrices, [Q, N0], N=Q.N0.Q*";
    local Q, N0, p1, p2, q1, q2, N11, N12, N21, N22;

    q1 := RowGCD(N[1,1], N[1,2]):
    p1 := q1*hc(q1):

    N11 := Lquo(N[1,1], p1):
    N12 := Lquo(N[1,2], q1):
    N21 := Lquo(N[2,1], hc(q1)):

    q2 := RowGCD(N[2,2], N21):
    p2 := q2*hc(q2):
    N22 := Lquo(N[2,2], p2):
    N21 := Lquo(N21, q2):
    N12 := Lquo(N12, hc(q2)):

    N0:= Matrix([[N11, N12], [N21, N22]]):
    Q := Matrix([[q1, 0], [0, q2]]):


    return Q, N0:

end proc:


RowGCD := proc(N11, N12)
    description "See Algo 3.2.1.";

    local tmp, tmp2, p, q1, q2, ptilde, q;

    tmp := collect(N12*hc(N12), z):
    tmp := Lpoly2poly(tmp):
    tmp2 := Lpoly2poly(N11):
    p := gcd(tmp, tmp2):
    q1 := gcd(p, Lpoly2poly(N12)):
    ptilde := Lquo(q1*hc(q1), p):
    ptilde := NormRealonT(ptilde):
    q2 := SqrtFejerRiesz(ptilde):
    q := Lquo(q1, q2):

    return q:

end proc:


solve_new := proc(eqns, x::evaln)
    # x is the name of the free variable
    local res, S, R, f, sol, df:
    res:= solve(eqns):
    S, R := selectremove(evalb, res):   # S is free var, R is the rest var

    f:= {seq(lhs(S[j])=x[j], j=1..nops(S) )}:
    sol := eval(R, f) union f:
    df := nops(S):  # degree of freedom

    return sol, df:

end proc:



DOS_compare:= proc(f1, f2, fv)
    description "Compute |f1|^2 - |f2|^2, return col vector of <f1, hc(f2), DOS>";
    # f1 and f2 are Lpoly with unknown parameters, fv are the values of the unknown params

    local g1, g2, v:

    g1 := eval(f1, fv):
    g2 := collect(hc(eval(f2, fv)),z):
    v := simplify(g1*hc(g1) - g2*hc(g2)):

    return <g1, g2, v>:
end proc:


################# Indefinite Split for 2x2 case ###################

SplitConstDet := proc(A)
    description "Spectral Decompose a 2x2 matrix of Lpoly, if the det is constant":
    # input should be a Hermitian 2x2 matrix of Lpoly, this is not checked in the code
    # also the input is assumed to be real
    # return UU, DD, such that:   A = U.DD.hermConj(U)

    local det, len11, len22, UU, DD, tmp, ii:

    det := simplify(Determinant(A)):
    if is(det=0) then
        error("Determinant is zero, singular case"):
    end if:

    if is(lenLpoly(det)>0) then
        error("Determinant is not constant!"):
    end if:

    len11:=lenLpoly(A[1,1]):
    len22:=lenLpoly(A[2,2]):

    # initialization
    if is((len11-len22)>0) then
        UU:= Matrix([[0, 1],[1, 0]]):
        DD:= simplify(UU.A.UU):
    else
        UU:= Matrix([[1,0], [0, 1]]):
        DD:= A:
    end if:

    # loop while the of diagonal element is not zero
    for ii from 1 while not is(DD[2,1] = 0) do
        if not is(DD[1,1]=0) then      # DD[1,1] element is not zero
            len11:=lenLpoly(DD[1,1]):
            len22:=lenLpoly(DD[2,2]):
            if is((len11-len22)>0) then     # swich the two diagonal elements
				tmp := Matrix([[0, 1],[1, 0]]):
                UU:= simplify(UU. tmp):
                DD:= simplify(tmp.DD.tmp):
                next:   # similar to "continue"
            end if:

            if is(Lquo(DD[2,1], DD[1,1]) = 0) then
                error("Should not happen for constant det matrix!"):
            end if:

            tmp:= Matrix([[1, 0], [-Lquo(DD[2,1], DD[1,1]), 1]]):
            DD := simplify(tmp.DD.hermConj(tmp)):
            UU := simplify(UU.MatrixInverse(tmp)):

        else    # D[1,1]=0, then off-diag elements must be const
            if is(lenLpoly(DD[2,1])>0) then
                error("Should not happen for constant det matrix!"):
            end if:

            tmp:= Matrix([[1, 0], [ -Lquo(DD[2,2], DD[1,2])/2, 1]]):
            tmp := Matrix([[1, 1], [1, -1]]).tmp:   # assume the matrix is real, otherwise this is matrix cannot diagonalize DD

            DD := simplify(tmp.DD.hermConj(tmp)):
            UU := simplify(UU.MatrixInverse(tmp)):

            break:

        end if:
    od:

    tmp:= Matrix([[sqrt(abs(DD[1,1])), 0], [0, sqrt(abs(DD[2,2]))]]):
    UU:= simplify(UU.tmp):
    tmp := MatrixInverse(tmp):
    DD:= simplify(tmp.DD.hermConj(tmp)):

    return UU, DD:    
end proc:

SplitConstDet_RealSym := proc(A)
    description "Spectral Decompose a 2x2 matrix of Lpoly, if the det is constant":
    # input should be a Hermitian 2x2 matrix of Lpoly, this is not checked in the code
    # also the input is assumed to be real
    # return UU, DD, such that:   A = U.DD.hermConj(U)

    local det, len11, len22, UU, DD, tmp, ii:

    det := simplify(Determinant(A)):
    if is(det=0) then
        error("Determinant is zero, singular case"):
    end if:

    if is(lenLpoly(det)>0) then
        error("Determinant is not constant!"):
    end if:

    len11:=lenLpoly(A[1,1]):
    len22:=lenLpoly(A[2,2]):

    # initialization
    if is((len11-len22)>0) then
        UU:= Matrix([[0, 1],[1, 0]]):
        DD:= simplify(UU.A.UU):
    else
        UU:= Matrix([[1,0], [0, 1]]):
        DD:= A:
    end if:

    # loop while the of diagonal element is not zero
    for ii from 1 while not is(DD[2,1] = 0) do
        if not is(DD[1,1]=0) then      # DD[1,1] element is not zero
            len11:=lenLpoly(DD[1,1]):
            len22:=lenLpoly(DD[2,2]):
            if is((len11-len22)>0) then     # swich the two diagonal elements
				tmp := Matrix([[0, 1],[1, 0]]):
                UU:= simplify(UU. tmp):
                DD:= simplify(tmp.DD.tmp):
                next:   # similar to "continue"
            end if:

            if is(Lquo(DD[2,1], DD[1,1]) = 0) then
                error("Should not happen for constant det matrix!"):
            end if:

            tmp:= Matrix([[1, 0], [-LquoSym(DD[2,1], DD[1,1]), 1]]):
            DD := simplify(tmp.DD.hermConj(tmp)):
            UU := simplify(UU.MatrixInverse(tmp)):

        else    # D[1,1]=0, then off-diag elements must be const
            if is(lenLpoly(DD[2,1])>0) then
                error("Should not happen for constant det matrix!"):
            end if:

            tmp:= Matrix([[1, 0], [ -LquoSym(DD[2,2], DD[1,2])/2, 1]]):
            tmp := Matrix([[1, 1], [1, -1]]).tmp:   # assume the matrix is real, otherwise this is matrix cannot diagonalize DD

            DD := simplify(tmp.DD.hermConj(tmp)):
            UU := simplify(UU.MatrixInverse(tmp)):

            break:

        end if:
    od:

    tmp:= Matrix([[sqrt(abs(DD[1,1])), 0], [0, sqrt(abs(DD[2,2]))]]):
    UU:= simplify(UU.tmp):
    tmp := MatrixInverse(tmp):
    DD:= simplify(tmp.DD.hermConj(tmp)):

    return UU, DD:    
end proc:
