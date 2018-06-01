######################################################
#
#   Library related to Laurent Polynomials
#
#   This should be the lowest level library, other
#   libraries will depend on the functions here.
#
#   Chenzhe
#   Jun, 2018
#


with(LinearAlgebra):
with(ArrayTools):
with(combinat):
with(PolynomialTools):


############ Basic Tools ###############
#   Hermitian Conjugate, length, etc

img := proc(f)
    coeff(simplify(f), I):
end proc:

reaL := proc(f)
    simplify(f - img(f)*I ):
end proc:


hc := proc(f)
    ## returns \bar f
    ## here we assume 'z' is the variable and parameters are real
	# sometimes wrong if not calling collect(f,z) first, Chenzhe

    simplify( eval( f - 2*img(f)*I, z = 1/z) ):
end proc:

hermconj := proc(f)
    #hc(f):		# commented by Chenzhe, sometimes wrong if not calling collect(f,z) first
    #eval(f, {z = 1/z, I = -I}):	# Chenzhe
    conjugate( eval(f, z = 1/conjugate(z)) ):
end proc:



hermConj := proc(A)
    local i, j, total_elements, m, n, lst, B:

    lst := [op(1, A)];
    if ( isVectorOrMatrix(A) ) then
	m := lst[1]: n := lst[2]:
	B := Transpose(A):
	for i from 1 to m do
	    for j from 1 to n do
		B[j, i] := hermconj(A[i, j]):
	    od:
	od:
    else ### assume total_elements == 1, scalar case
	B := hermconj(A):
    end if:
    B:
end proc:


Lpoly2poly:= proc(f)
    description "Shift Laurent poly to polynomial";
	# since quo() only accept polynomial input

    local ldeg, s;
    ldeg:= -ldegree(f, z):
    s:= collect(f*z^ldeg,z):
	
	# if called as n,a:=Lpoly2poly(f), return ldegree and s separately
	# otherwise, only return s
	if _nresults = 1 or _nresults = undefined then
		return s:
	else
		return -ldeg, s:
	fi:
end proc:

lenLpoly := proc(f)
	description "length of Laurent polynomial";
	local ldeg, deg, len;

    if is(simplify(f)=0) then
        return -infinity:
    end if:

	deg:= degree(f):
	ldeg:= ldegree(f):
	len:= deg-ldeg:
	
	return len:
end proc:


fsupp := proc(A)
    description "A is a Matrix/scalar of Lpoly, return the fsupp of each elements also as a matrix/scalar":
    local m, n, k, l, tmp, out:

    if not type(A, scalar) then
        m, n := Dimension(A):
        out:= Matrix(m, n):

        for k from 1 to m do
            for l from 1 to n do
                tmp:=fsupp(A[k, l]):
                out[k, l] := tmp:
            od:
        od:
    else    # scalar case
        out:= {ldegree(A,z), degree(A, z)}:
    end if:

    return out:
end proc:

printfilter := proc(a)
    local v, supp, ldeg:

    interface(displayprecision = 8):
    v:= RoundingLpoly(a):
    supp := fsupp(v):
    ldeg:= ldegree(v, z):

    v:= CoefficientList(simplify(v/z^ldeg), z):

    print(v, supp)
end proc:

########## Division of Lpoly #################
# Lpoly version of quo, rem, divide, gcd, etc.
#

Lquo:= proc(a, b)
	description "quo of Laurent Polynomial";

	local az, na, bz, nb, q:
	na, az := Lpoly2poly(a):
	nb, bz := Lpoly2poly(b):
	q:= quo(az, bz, z):
	q:= collect(q*z^(na-nb),z):
	
	return q:

end proc:

LquoSym := proc(a, b)
	description "For Lpoly a, b both with (real) symmetry, output q=Lquo(a, b) would also be Lpoly with symmetry":
	# a = b*q + r, symmetry type: S(a) = S(b)*S(q) = S(r)
	
	local Sa, Sb, Sq, ii, az, bz, f, q:
	if not isRealSym(a) then
		error("Input a without (real) symmetry"):
	end if:
	if not isRealSym(b) then
		error("Input b without (real) symmetry"):
	end if:
		
	Sa := RealSymType(a):
	Sb := RealSymType(b):
	Sq := simplify(Sa/Sb):
	az := collect(a, z):
	bz := collect(b, z):
	q := 0:
	for ii from 1 while lenLpoly(az)>=lenLpoly(bz) do
		f:= lcoeff(az, z)/lcoeff(bz, z) * z^(degree(az, z) - degree(bz, z)):
		f:= f + Sq*eval(f, z=1/z):
		az := collect(az - f*bz, z ):
		q:= q + f:
	od:
	
	return collect(q, z):
end proc:

Lrem:= proc(a, b)
    description "rem of Laurent Polynomial";

    local q, r:
    q:= Lquo(a, b):
    r:= a - b*q:

    return simplify(r):
end proc:


Ldivide := proc(a, b)
    description "Check if b|a, where a and b are Laurent poly of z. Return is true/false";
	# similar to divide() function, which sometimes gives wrong results for complex Laurent poly.
	
	local apoly, bpoly, r:
	apoly:= Lpoly2poly(a):
	bpoly:= Lpoly2poly(b):
	r:= rem(apoly, bpoly, z):
	
	return evalb(r=0):
	
end proc:

Lgcd := proc(a, b)
    description "gcd of Laurent poly";

    local pa, pb, d:
    pa := Lpoly2poly(a):
    pb := Lpoly2poly(b):

    d := gcd(pa, pb):
    return d:

end proc:


################ Symmetry and normalization of Lpoly ##########

isRealSym := proc(f)
    description "Check if a Lpoly is real symmetric, return true/false";

    local f1, finv:

    f1:= collect(f, z):
    finv := eval(f, z=1/z):

    return Ldivide(f, finv):

end proc:


isComplexSym := proc(f)
	description "Check if a Lpoly is (essential) complex symmetric, return true/false";
	
	local f1, fstar;
	f1:= collect(f, z):
	fstar:=hc(f):
	
	# return is a logical true/false
	return Ldivide(f, fstar):
end proc:


isRealonT := proc(p)
	description "Check p is real on T";
    # p is real on T <==> complex sym factor SS=1
	
	local r1, r2:
	r1:= isComplexSym(p):
	r2:= evalb(Lquo(p, hc(p))=1 ):
	return evalb(r1 and r2):
	
end proc:

NormComplexSym := proc(f)
    description "If f has essential complex sym, normalize it to make the complex sym factor 1 or z";

    local fz, fz_star, epsilon, c, SS, tmp, p, ctr:

    if is(f=0) then
        return 0:
    end if:

    if not isComplexSym(f) then
        error("The input is not essential complex symmetric!"):
    end if:

    fz:= collect(f, z):
    fz_star:= collect(hc(f),z):

    SS:= Lquo(fz, fz_star, z):
	c:=degree(SS):
	epsilon:= coeff(SS, z, c):

    ctr:= floor(c/2);
	tmp:= z^(-ctr)/sqrt(epsilon):

	p:=collect(simplify(f*tmp),z):

	return p:

end proc:

NormRealonT := proc(f)
    description "For f with complex sym factor: epsilon*z^c, c is even. We normalize f to be real on T";
    # If it not possible to normalize f to be real on T, generate error

    local fz, fz_star, epsilon, c, SS, tmp, p:
	
	if modp(lenLpoly(f),2)<>0 or (not isComplexSym(f)) then
		error("Cannot be normalized to real valued function on T!"):
	end if:
	
    # p is real on T <==> complex sym factor SS=1
    p:= NormComplexSym(f):

	return p:

end proc:


RealSymType := proc(A)
    description "A is a matrix/scalar of Lpoly with (real) sym, return the sym type of each ele also as a matrix/scalar";
    local m, n, k, l, tmp, out:

    if not type(A, scalar) then
        m, n := Dimension(A):
        out:= Matrix(m, n):

        for k from 1 to m do
            for l from 1 to n do
                tmp:=RealSymType(A[k, l]):
                out[k, l] := tmp:
            od:
        od:

    else    # scalar case
        if not isRealSym(A) then
            error("The Laurent polynomial does not have (real) symmetry"):
        end if:

        out:= Lquo(A, eval(A, z=1/z)):
    end if:

    return out:
end proc:

ComplexSymType := proc(A)
    description "A is a matrix/scalar of Lpoly with complex sym, return the sym type of each ele also as a matrix/scalar";
    local m, n, k, l, tmp, out, fz, fz_star:


    if not type(A, scalar) then
        m, n := Dimension(A):
        out:= Matrix(m, n):

        for k from 1 to m do
            for l from 1 to n do
                tmp:=ComplexSymType(A[k, l]):
                out[k, l] := tmp:
            od:
        od:
    else    # scalar case
        if not isComplexSym(A) then
            error("The Laurent polynomial does not have complex symmetry"):
        end if:

        fz:= collect(A, z):
        fz_star:= collect(hc(fz),z):

        out:= Lquo(fz, fz_star, z):

    end if:


    return out:
end proc:


################ Smith Normal Form ########################

FactorLMatrix := proc(M)    # no use: we could just use factor~(M) instead of this function
    description "Factor each elements of matrix of Lpoly":
    local N, ii,jj, m, n:

    m, n := Dimension(M):
    N := Matrix(m, n):
    for ii from 1 to m do
        for jj from 1 to n do
            N[ii, jj] := factor(M[ii, jj]):
        od:
    od:

    return N:
end proc:

(*
LSmithForm_old := proc(M)   # old implementation, not using the generic SmithForm
    description "Smith Normal Form of a matrix of Laurent polynomials. Return SS, UU, VV: SS=UU.M.VV":
    # notice: avoid S, U, V in the workspace when calling this function, the names are used for SmithForm()
    # Could we find a way to get rid of this trouble?

    local m, n, ii, jj, ld, Mnew, SS, UU, VV, NV, N:

    m, n := Dimension(M):
    ld:=ldegree(M[1,1],z):

    for ii from 1 to m do
        for jj from 1 to n do
            ld := min(ldegree(M[ii,jj]), ld):
        od:
    od:
    Mnew:=M/z^ld:
    Mnew:=FactorLMatrix(Mnew):   
    # SmithForm() only works well after factorization, 
    # otherwise will get wrong answer if the Lpoly is complex

    SS, UU, VV:= SmithForm(Mnew,z,output=['S', 'U', 'V']):
    
    UU:=UU/z^ld:
    NV := Vector(m):
    for ii from 1 to m do
        ld:=ldegree(SS[ii,ii], z):
        NV[ii] := 1/z^ld/lcoeff(SS[ii, ii], z):
    od:
    N := DiagonalMatrix(NV):
    UU:= N.UU:
    SS:= N.SS:

    SS:=FactorLMatrix(SS):

    return SS, UU, VV: # S=U.M.V

end proc:
*)

LSmithForm:= proc(M)
    description "Smith Normal Form of a matrix of Laurent polynomials. Return SS, UU, VV: SS=UU.M.VV":
    # notice: avoid S, U, V in the workspace when calling this function, the names are used for SmithForm()
    # Could we find a way to get rid of this trouble?

    local m, n, ii, jj, ld, Mnew, SS, UU, VV, NV, N, Cz:

    Cz[`0`]:=0: 
    Cz[`1`]:=1: 
    Cz[`+`]:=`+`: 
    Cz[`-`]:=`-`: 
    Cz[`*`]:=`*`: 
    Cz[`=`]:= `=`:
    Cz[Quo] := proc(a,b,r) if nargs=3 then quo(a,b,z,r) else quo(a,b,z) fi end:
    Cz[Rem] := (a,b,q) -> rem(a,b,z,q):
    Cz[EuclideanNorm] := a -> degree(a,z):
    Cz[Gcdex] := (a,b,s,t)->gcdex(a,b,z,s,t):
    Cz[UnitPart] := sign:

    m, n := Dimension(M):
    ld:=ldegree(M[1,1],z):

    for ii from 1 to m do
        for jj from 1 to n do
            ld := min(ldegree(M[ii,jj]), ld):
        od:
    od:
    Mnew:=simplify(M/z^ld):

    SS, UU, VV := LinearAlgebra:-Generic:-SmithForm[Cz](Mnew, output=['S','U','V']);


    UU:=simplify(UU/z^ld):
    SS:= simplify(SS):
    NV := Vector(m):
    for ii from 1 to m do
        ld:=ldegree(SS[ii,ii], z):
        NV[ii] := 1/z^ld/lcoeff(SS[ii, ii], z):
    od:
    N := DiagonalMatrix(NV):
    UU:= simplify(N.UU):
    SS:= N.SS:

    SS:=factor~(SS):
    return SS, UU, VV: # S=U.M.V

end proc:

Matrixdiff := proc(M, x)
    local m, n, ii, jj, out:
    m, n := Dimension(M):
    out := Matrix(m, n):

    for ii from 1 to m do
        for jj from 1 to n do
            out[ii,jj] := diff(M[ii,jj], x):
        od:
    od:

    return out:

end proc:

Lsimplify := proc(p)
	description "Simplify a Laurant polynomial":
	# simplify each coeff in the Lpoly, somehow simplify() function cannot work this way by itself.
	local q, qout, ldeg, deg, iter, c:
	
	q:= collect(p, z):
	qout:= 0:
	ldeg:= ldegree(p, z):
	deg:= degree(p, z):
	
	for iter from ldeg to deg do
		c:= coeff(q, z, iter):
		c:= simplify(c):
        c:= simplify(convert(c, radical)):
		qout:= qout + c*z^iter:
	od:
	
	return qout:
	
end proc:

RoundingLpoly := proc(p)
    description "Return evalf(p), dropping coeffs smaller than 10e-(Digits -2 )":

    local iter, deg, ldeg, q, qout, CompC, tol:

    deg := degree(p, z):
    ldeg := ldegree(p, z):

    q:= evalf(collect(p,z)):
    qout:=0:

    tol:= 10^(-Digits+2):
    for iter from ldeg to deg do
        CompC:= coeff(q, z, iter):
        if abs(Re(CompC)) < tol then
            CompC:= I*Im(CompC):
        fi:
        if abs(Im(CompC)) < tol then
            CompC:= Re(CompC):
        fi:

        qout := qout + CompC * z^iter:
    od:

    return qout:

end proc:


SmoothExp2 := proc(a)
    description "Smooth Exponent of lowpass filter a":

    local sm, v, w, rho, sr, K, j, k, W, eig, indrule:
    sr := SumRule(a):

    v:= Lquo(a, (1+z)^sr):
    v:= RoundingLpoly(v):
    w:= RoundingLpoly(collect(v*hc(v), z)):

    K := degree(w, z):

    indrule := (j, k, w) -> coeff(w, z, (2*j-k) ):

    W := Matrix( 2*K+1, (i,j)-> indrule(i-K-1, j-K-1, w) ):

    eig := abs~(Eigenvalues(W)):
    rho := max(eig):

    sm := -0.5 - log[2](sqrt(rho)):

    return sm:

end proc:

Lpoly2List := proc(a)
    description "Print the Lpoly a as a list, print its fsupp":
    local aList, fs, deg, ldeg, j:

    deg := degree(a, z):
    ldeg := ldegree(a, z):
    aList := []:
    fs := [ldeg, deg]:

    for j from ldeg to deg do
        aList := [op(aList), coeff(a, z, j)]:
    end do;

    return aList, fs:
end proc:


CheckPR := proc(fb, Sig)
    description "Check PR of 1D dilation 2 filter bank":
    # fb is a column vector, Sig is the diagonal signature matrix with +/- 1
    # if PR holds, should return 2x2 identity matrix
    local FB:

    FB := <fb| eval(fb, z=-z)>:

    return simplify~(hermConj(FB).Sig.FB):
end proc:


LpolySym2Poly := proc(Lpoly)
    description "Given a Lpoly with real sym, and sym center is 0, write it as a poly of zeta = (z+1/z)";
    local Poly, c, deg, k, Lpoly1:

    #evaln(z):
    deg := degree(Lpoly, z):
    Lpoly1 := Lpoly:
    Poly := 0:

    for k from deg by (-1) to 1 do
        c:= coeff(Lpoly1, z, k):
        Poly := Poly + c * zeta^k:
        Lpoly1 := Lpoly1 - c * (z+1/z)^k:
    end do;

    Poly := Poly + simplify(Lpoly1):

    if evalb(simplify(eval(Poly, zeta = (z+1/z)) - Lpoly) <> 0 ) then
        error("Error in LpolySym2Poly, no real sym?");
    end if;

    return Poly:
end proc:



##################### for high dimension, maybe move to somewhere else?

getAllTerms := proc(Lpoly, vlist)
    description "Get all terms in a Lpoly as a list, works for multivariate case":

    local a, termList:

    a := collect(Lpoly, vlist, distributed):
    if type(a, `+`) then
        termList := [op(a)]:
    else
        termList := [a]:
    end if;

    return termList:
end proc:

degreeVec := proc(t, vlist)
    description "Input a multivariate polynomial with variables in vlist, return its degree as vector":

    local n, degV, j:

    n := nops(vlist):
    degV := Vector(n):
    for j from 1 to n do
        degV[j] := degree(t, vlist[j]):
    end do;

    return degV:

end proc:

isConstTerm := proc(t, vlist)
    description "Input a multivariate monomial, determine if it is a contant term":

    local degV, deg:

    if evalb(simplify(t) = 0) then
        return true:
    fi:

    degV := degreeVec(t, vlist):
    deg := add(abs(degV[j]), j = 1..nops(vlist)):

    return evalb(deg = 0):
end proc:

hcND := proc(a, vlist)
    local s, j:
    s := {}:
    for j from 1 to nops(vlist) do
        s := s union {vlist[j] = 1/vlist[j]}:
    end do;

    return eval(a, s):
end proc: