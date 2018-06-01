############# 3 Types of Fejer Riesz ##############
# 1. No symmetry
# 2. Complex symmetry
# 3. Real symmetry
#
#   This library relies on libLpoly.mpl, libMatrixSplit.mpl
#
#   Chenzhe
#   Jun, 2018
#


with(LinearAlgebra):
with(ArrayTools):
with(combinat):
with(PolynomialTools):


SqrtFejerRiesz := proc(p)
	description "Square Root of a Laurent poly p, one solution of F-R lemma. p=(+/-)q*hc(q)":
	# This proc is symbolic, which means the coeffs of input p has to be fractions, rather than floating point numbers.
	# Otherwise, the PolynomialTools[Splits] will fail.
    # This proc works for input p always positive/negative on T, otherwise, it will generate error message.
    #   If p >= 0 on T, p=q*hc(q)
    #   If p <= 0 on T, p=-q*hc(q)

    # We only choose the roots that are:
    # 1. inside T
    # 2. half of the ones on the unit circle T.
	
    local q, c, f, fr, nr, f0, ii, r, mr, lambda, tol:


    if is(p=0) then
        return 0:
    end if:

    if not isRealonT(p) then
        error("Cannot be factorized by F-R lemma, the Lpoly is not real on T!"):
    end if:

	f:= Splits(p, z):
	c:= f[1]:	# leading coeff of factorization
	fr:= f[2]:	# all roots and their multiplicity
	nr:= nops(fr):	# number of distinct roots
	
	tol:=10^(-Digits+2):
	q:=1:
	for ii from 1 to nr do
		f0:= fr[ii]:	# factor and its multiplicity
		r:= eval(-f0[1], z=0):		# the root
		mr:= f0[2]:				# its multiplicity
		if (1-evalf(abs(r)))> tol   then	# inner roots
			q:= q*f0[1]^mr:
		elif abs(evalf(abs(r))-1)<=tol    then	# roots on T
			if modp(mr, 2)=0    then
				q:=q*f0[1]^(mr/2):
			else
				error("Cannot be factorized by F-R Lemma, the multiplicity of roots on T is not even!"):
			end if:
		end if:
		
	od:

    lambda:= eval(q, z=0):
    lambda:= sqrt(lcoeff(p)/lambda);
    q:=simplify(q*lambda):
	
    return q:
	
end proc:


SqrtComplexSym := proc(p)
	description "Special case of F-R lemma, p=q*hc(q), where q has complex symmetry";
	# This proc only checks that the input p is real on T.
    #
    #   Input p has to be always positive on T
    #   Otherwise, it will give wrong answer
    #
    # The output has complex sym factor 1 or z

    local N, c, q, t, j, k, n, tmp:

    if not isRealonT(p) then
        error("Cannot be factorized by Complex Sym F-R lemma, the Lpoly is not real on T!"):
    end if:

    N:= degree(p, z):
    if modp(N, 2)=0 then
        # Output has complex sym factor 1
        c:=sqrt(lcoeff(p, z)): # c=sqrt(abs(p_N)) * exp(i*alpha)
        n:= N/2:
        t[n]:= 1:
        for j from 1 by 1 to n do
            t[n-j] := coeff(p, z, N-j)/coeff(p, z, N):
            for k from (n-j+1) by 1 to (n-1) do
                t[n-j] := t[n-j] - t[k]*t[2*n-j-k]:
            od:

            t[n-j] := t[n-j]/2:
        od:

        q:= c * add(t[k]*z^k, k=1..n):
        q:= q + hc(q):
        t[0] := c*t[0]:
        t[0] := (t[0] + hc(t[0]))/2:
        q:= q + t[0]:
        
    else
        # Output has complex sym factor 1/z
        c:=sqrt(lcoeff(p, z)): # c=sqrt(abs(p_N)) * exp(i*alpha)
        n:= floor(N/2):
        t[n]:= 1:
        for j from 1 by 1 to n do
            t[n-j] := coeff(p, z, N-j)/coeff(p, z, N):
            for k from (n-j+1) by 1 to (n-1) do
                t[n-j] := t[n-j] - t[k]*t[2*n-j-k]:
            od:

            t[n-j] := t[n-j]/2:
        od:

        q:= c * add(t[k]*z^k, k=0..n):
        tmp:= hc(q)/z:
        q:= q+tmp:

    end if:

    if type(lcoeff(p,z), realcons) and (not type(c, realcons)) then
        q := q*I:
    end if:


    # Check the result
    if simplify(q*hc(q)-p)<>0 then
        print("The input  cannot be splitted as square of complex sym Laurent Poly"):
    end if:

    return simplify(q):



end proc:


SqrtRealSym_param := proc(p)
	description "Special case of F-R lemma, p=q*hc(q), where q has real symmetry, parameterized case";
	# This proc only checks nothing, p can be parametrized.
    #
    #   Find the solution using formula.
    #
    #   Output 2 solutions for the cases p_N positive or negative.
    #
    # The output has sym factor 1 or z

    local N, c, q1, q2, t, j, k, n, tmp:

    t := table():
    N:= degree(p, z):
    if modp(N, 2)=0 then
        # Output has real sym factor 1
        n:= N/2:
        t[n]:= 1:
        for j from 1 by 1 to n do
            t[n-j] := coeff(p, z, N-j)/coeff(p, z, N):
            for k from (n-j+1) by 1 to (n-1) do
                t[n-j] := t[n-j] - t[k]*t[2*n-j-k]:
            od:

            t[n-j] := t[n-j]/2:
        od:

        c:=sqrt(lcoeff(p, z)):  # p_N > 0
        q1:= c * add(t[k]*(z^k+1/z^k), k=1..n):
        t[0] := c*t[0]:
        q1:= q1 + t[0]:

        c:=sqrt(-lcoeff(p, z)): # p_N < 0
        q2 := c * add(t[k]*(z^k-1/z^k), k=1..n):
        
    else
        # Output has real sym factor 1/z
        n:= floor(N/2):
        t[n]:= 1:
        for j from 1 by 1 to n do
            t[n-j] := coeff(p, z, N-j)/coeff(p, z, N):
            for k from (n-j+1) by 1 to (n-1) do
                t[n-j] := t[n-j] - t[k]*t[2*n-j-k]:
            od:

            t[n-j] := t[n-j]/2:
        od:

        c:=sqrt(lcoeff(p, z)):  # p_N > 0
        q1:= c * add(t[k] * (z^k + 1/z^k/z), k=0..n):

        c:=sqrt(-lcoeff(p, z)):  # p_N < 0
        q2:= c * add(t[k] * (z^k - 1/z^k/z), k=0..n):

    end if:

    return simplify(q1), simplify(q2):

end proc:


(*
SqrtRealSym := proc(p)
	description "Special case of F-R lemma, p=(+/-)q*hc(q), where q has real symmetry";
	# This proc is symbolic, which means the coeffs of p has to be fractions, rather than floating point numbers.
	# Otherwise, the PolynomialTools[Splits] will fail.
    #
    # This proc checks:
    #   1. input p has to be always positive/negative on T
    #       If p >= 0 on T, p=q*hc(q)
    #       If p <= 0 on T, p=-q*hc(q)
    #   2. Necessary root condition for real sym split
    # Otherwise, it will generate error message
    #
    # We only choose the roots that are 
    #   1. inner T and above real axis
    #   2. outter T and below real axis
    #   3. Half of the ones on the unit circle T
    #   4. Half of the ones on real axis
	#
	#
	#
	#	Needs to be rewritten, with explicit formula
	#
	#	Chenzhe
	#	Nov, 2016
	#
	
    local q, c, f, fr, nr, f0, ii, r, mr, lambda, tol, len:


    if (not isRealonT(p)) or (not isRealSym(p)) or (Lquo(p, eval(p, z=1/z))<>1) then
        error("Cannot be factorized by Real Sym F-R lemma, the Lpoly is not real/complex sym!"):
        # this checks real sym with factor 1 and p has real coeffs
    end if:

	f:= Splits(p, z):
	c:= f[1]:	# leading coeff of factorization
	fr:= f[2]:	# all roots and their multiplicity
	nr:= nops(fr):	# number of distinct roots
	
	q:=1:
	tol:=10^(-Digits+2):
    for ii from 1 to nr do
        f0:= fr[ii]:	# factor and its multiplicity
		r:= eval(-f0[1], z=0):		# the root
        mr:= f0[2]:				# its multiplicity

        if evalf(Im(r))>=tol and (1-evalf(abs(r)))>=tol then
            # Case 1: upper & inner
            q:= q*f0[1]^mr:
        elif evalf(Im(r))<=-tol and (evalf(abs(r))-1) >= tol then
            # Case 2: lower & outter
            q:= q*f0[1]^mr:
        elif abs(evalf(abs(r))-1)<tol then
            # Case 3: on T
            if modp(mr, 2)=0    then
                q:=q*f0[1]^(mr/2):
            else
                error("Cannot be factorized by Real Sym F-R Lemma, the multiplicity of roots on T is not even!"):
            end if:
        elif evalf(abs(Im(r)))<tol then
            # Case 4: on real axis
            if modp(mr, 2)=0    then
                q:=q*f0[1]^(mr/2):
            else
                error("Cannot be factorized by Real Sym F-R Lemma, the multiplicity of roots on real line is not even!"):
            end if:

        end if:

    od:

    lambda:= eval(q, z=0):
    lambda:= sqrt(lcoeff(p)/lambda);
    len:=floor(lenLpoly(q)/2):
    q:=simplify(q*lambda/z^len):
	
    return q:

end proc:
*)

SqrtDaubechiesP0 := proc(P)
    description "Implementation of Algorithm 2.2.2";
    # P is a polynomial in x satisfiying P(x)>=0 for x on [0,1], and P(0) = 1.
    # P(x) should have real coefficients, output Q(z) also has real coeffs.
    #
    # The algorithm is similar to Fejer Riesz:
    #   Output Q(z) is a Laurent polynomial satisfying |Q(z)|^2 = P(sin^2(\xi/2)), Q(z=1)=1.
    #
    # If we take P(x) = P_{m,m}(x), it is used to construct Daubechies filters:
    #       |a(z)|^2 = cos^{2m}(\xi/2) P(sin^2(\xi/2))
    #       a(z) = 2^{-m} (z+1)^m Q(z)
    # Solving in this way has 2 benefits than using FR directly:
    #   1. we don't need to take care of the sum rule m, this reduces the degree of the polynomial
    #   2. polynomial P(x) has lower degree than P(sin^2(\xi/2)) in z domain (half). 
    #      Since both methods require us to factorize and find the roots, lower degree would be 
    #      much easier in symbolic case.



end proc:

SqrtDaubechiesP1 := proc(P)
    description "SqrtDaubechiesP() function with FR":

    local P1, Q, tmp, psin:

    psin:= (2-z-1/z)/4:
    P1:= eval(P, x=psin):
    Q := SqrtFejerRiesz(P1):
    tmp := eval(Q, z=1):
    Q := simplify(Q/tmp):

    Q := collect(Q, z):

    return Q:

end proc:

SqrtDaubechiesP_num := proc(P)
    description "SqrtDaubechiesP() function with numerical output":

    local P1, Q, tmp, psin:

    psin:= (2-z-1/z)/4:
    P1:= eval(P, x=psin):
    Q := SqrtFejerRiesz(P1):
    tmp := eval(Q, z=1):
    Q := simplify(Q/tmp):

    Q := evalf(collect(Q, z)):
    Q := reaL(Q):

    return Q:

end proc:

DaubechiesLowpass := proc(m)
    description "construct Daubechies Lowpass filter of order m, using SqrtDaubechiesP() function above":
    local P, Q, a:

    P := Pmn(m, m, x):
    Q := SqrtDaubechiesP_num(P):
    Q := collect(Q, z):
    a := 2^(-m) * z^(1-m) * (1+z)^m * Q;

    return a:

end proc:





