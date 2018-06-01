##### Copyright and Notice #####
#Copyright by Bin Han at University of Alberta, since 2018
#All Rights Reserved
#Report errors, mistakes, and bugs to Bin Han at bhan@ualberta.ca.
#Send comments/suggestions to Bin Han at bhan@ualberta.ca

#This maple routines are based and for the book:
#Bin Han, Framelets and Wavelets: Algorithms, Analysis,
#    and Applications, in Applied and Numerical Harmonic Analysis,
#    Birkhauser/Springer, Cham, (2017), 724 pages.

#The use of all the following software is permitted for
#non-commercial, educational, and research use only. The software
#and/or related materials are provided "as-is" without warranty of any
#kind including any warranties of performance or merchantability or
#fitness for a particular use or purpose or for any purpose
#whatsoever, for the licensed product, however used. In no event shall
#University of Alberta and/or Bin Han be liable for any damages and/or
#costs, including but not limited to incidental or consequential
#damages of any kind, including economic damage or injury to property
#and lost profits, regardless of whether University of Alberta shall
#be advised, have reason to know, or in fact shall know of the
#possibility. User bears all risk relating to quality and performance
#of the software and/or related materials.

#Any use other than non-commercial, educational, or research, or any
#redistribution in original or modified form requires prior written
#authorization from the copyright holder.
##################################################################

lprint("0. Copyrighted and developed by Bin Han at University of Alberta, Canada.  Send questions to bhan@ualberta.ca"):
lprint("1. maple routines can handle any dilation, any multiplicity of real/complex-valued filters/masks w/o (complex) symmetry."):
lprint("2. Do not assign values to protected variables: x,y,z,xi,zeta,eta."):
lprint("3. A mask/filter $a={a(k)}_{k in Z}:Z->C^{r*r}$ is stored by a matrix $sum_{k in Z}a(k)z^k$ of Laurent polynomials. 1 must be an eigenvalue of $sum_{k in Z} a(k)$ for low-pass filters."):
lprint("4. Routines ending with C for filters with complex-valued parameters. For complex filters without unknowns/parameters, use other routines as usual."):
lprint("5. For related theory, see book: Bin Han, Framelets and Wavelets, Birkhauser/Springer, (2017), 724 pages."):

protect(x,y,z,xi,zeta,eta):
with(LinearAlgebra):
Digits:=128:

##Generalize command subs for matrices of Laurent polynomials
##subsset: the set used for substitution
##mpoly: Matrix of (Laurent) polynomials or any thing
Subs:=proc(subsset,mpoly)
    local j,k,rowdim,coldim,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    result:=Matrix(rowdim,coldim,0):
    for j from 1 to rowdim do
    for k from 1 to coldim do
        result[j,k]:=expand(subs(subsset,mpoly[j,k])):
    end do: end do:
    result:=result:
end proc:


##Find support of coefficients of a matrix of Laurent polynomials
##Output: [low,high], a Vector for support of mask/filter.
##Input: a matrix of (Laurent) polynomials in terms of variable var.
D1SupportPoly:=proc(mpoly,var)
    local j,k,rowdim,coldim,tmppoly,result:
    result:=Vector[row](2,[10000000,-10000000]):
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    for j from 1 to rowdim do
    for k from 1 to coldim do
        tmppoly:=collect(expand(mpoly[j,k]),var):
        result[1]:=min(result[1],-degree(subs({var=1/var},tmppoly),var)):
        result[2]:=max(result[2],degree(tmppoly,var)):
    end do: end do:
    result:=result:
end proc:

##Write a fraction of complex numbers into real+I*imaginary
##maple does not automatically do this.
##This routine is useful for a fraction of two complex numbers
SimplifyComplex:=proc(cnum)
    local ncf,dcf,result:
    ncf:=numer(cnum):
    dcf:=denom(cnum):
    if ((dcf<>1)and(coeff(dcf,I,1)<>0)) then
        result:=expand(ncf*(coeff(dcf,I,0)-I*coeff(dcf,I,1))):
        result:=expand(result/simplify((coeff(dcf,I,0))^2+(coeff(dcf,I,1))^2)):
    else
        result:=cnum:
    end if:
    result:=result:
end proc:

##Simplify Poly with complicated complex coefficients
##so that each coefficient is expressed as real+I*imaginary
D1SimplifyComplexPoly:=proc(mpoly,var)
    local j,k,n,supp,rowdim,coldim,tmp,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    supp:=D1SupportPoly(mpoly,var):
    result:=Matrix(rowdim,coldim,0):
    for j from 1 to rowdim do
    for k from 1 to coldim do
    for n from supp[1] to supp[2] do
        tmp:=simplify(coeff(expand(mpoly[j,k]),var,n)):
        tmp:=SimplifyComplex(tmp):
        result[j,k]:=result[j,k]+tmp*var^n:
    end do:
        result[j,k]:=collect(result[j,k],var):
    end do: end do:
    result:=result:
end proc:


##Simplify/expand/factorize/evaluate coefficients of a Matrix poly
##choice=0 for simplify; choice=-1 for expansion, choice=-2 for factorization; otherwise, evaluate it to digits given by choice.
D1SimplifyPoly:=proc(mpoly,choice)
    local j,k,n,rowdim,coldim,tmppoly,ncf,dcf,result:
#    lprint("choice=0 for simplification; -2 for factorization; -1 for expansion;  otherwise for evaluating to choice digits."):
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    result:=Matrix(rowdim,coldim,0):
    for j from 1 to rowdim do
    for k from 1 to coldim do
        tmppoly:=expand(mpoly[j,k]):
        if (choice=0) then
            result[j,k]:=simplify(tmppoly):
        elif (choice=-1) then
            result[j,k]:=expand(tmppoly):
        elif (choice=-2) then
            result[j,k]:=factor(tmppoly):
        else
            result[j,k]:=evalf(tmppoly,choice):
        end if:
    end do: end do:
    result:=result:
end proc:


##solve a set of equations, if no solutions, give warning and return the same filters;
##if there are solutions, substitute the first solution to a poly (Matrix)
SolveEqs:=proc(seteqns,mpoly)
    local solt,result:
    solt:=solve(seteqns):
    if (solt=NULL) then
        print("In SolveEqs(seteqns,mpoly): No solutions to the given set of equations in seteqns!"):
        result:=mpoly:
    else
        result:=D1SimplifyPoly(Subs(solt,mpoly),-1):
    end if:
    result:=result:
end proc:


##Output: Matrix u(k) of a matrix filter at position k.
D1MaskElement:=proc(mpoly,xposition,var)
    local j,k,rowdim,coldim,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    result:=Matrix(rowdim,coldim,0):
    for j from 1 to rowdim do
    for k from 1 to coldim do
        result[j,k]:=simplify(coeff(mpoly[j,k],var,xposition)):
    end do: end do:
    result:=result:
end proc:


##Find the real part and imaginary part of a complex filter
##result[1]:=real part, result[2]:=imaginary.
D1RealImagPoly:=proc(mpoly)
    local result:
    result:=Vector(2):
    result[1]:=D1MaskElement(mpoly,0,I):
    result[2]:=D1MaskElement(mpoly,1,I):
    result:=result:
end proc:


##Check whether a filter is complex-valued or real-valued.
##This will affect the parametrization of its dual filters,
##\hat{\phi}, \vgu, etc. If a filter has real coefficients,
##then use real numbers as well
##Ouput: true for real-valued filters, false for complex-valued
D1IsRealFilter:=proc(mpoly)
    local j,k,rowdim,coldim,val,tmp,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    tmp:=D1MaskElement(mpoly,1,I): #imaginary parts
    val:=0:
    for j from 1 to rowdim do
    for k from 1 to coldim do
        val:=val+abs(tmp[j,k]):
    end do: end do:
    if (val=0) then
        result:=true:
    else
        result:=false:
    end if:
    result:=result:
end proc:


##Output: set of all coefficients in a matrix of Laurent polynomials in variable var
##If the coefficients are complex numbers, then
##instead of output the complex coefficients, we output
##its real part and imaginary part instead.
D1PolyToCoeffSet:=proc(mpoly,var)
    local j,k,n,rowdim,coldim,supp,tmppoly,cf, result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    supp:=D1SupportPoly(mpoly,var):
    result:={}:
    for j from 1 to rowdim do
    for k from 1 to coldim do
        tmppoly:=expand(mpoly[j,k]):
        for n from supp[1] to supp[2] do
            cf:=expand(coeff(tmppoly,var,n)):
            result:=result union {coeff(cf,I,0),coeff(cf,I,1)}:
        end do:
    end do: end do:
    result:=result:
end proc:


##Output: set of equations for Big O of a matrix of polynomials
##        up to the degree degord
##Explanation: mpoly(var)=O(|var|^degord) as var->0.
##For choice<=0, find set of all coefficients up to \xi^(degord-1).
##if choice>=0, find set of all coefficients with \xi^choice.
##If the coefficients are complex numbers, then
##instead of output the complex coefficients, we output
##its real part and imaginary part instead.
D1BigOEqs:=proc(mpoly,var,degord,choice)
    local j,k,n,rowdim,coldim,supp,cf,tmppoly,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    supp:=D1SupportPoly(mpoly,var):
    result:={}:
    for j from 1 to rowdim do
    for k from 1 to coldim do
        tmppoly:=expand(mpoly[j,k]):
        if (choice<0) then
            for n from max(0,supp[1]) to min(degord-1,supp[2]) do
                cf:=expand(coeff(tmppoly,var,n)):
                result:=result union {coeff(cf,I,0),coeff(cf,I,1)}:
            end do:
        else
            cf:=expand(coeff(tmppoly,var,choice)):
            result:=result union {coeff(cf,I,0),coeff(cf,I,1)}:
        end if:
    end do: end do:
    result:=result:
end proc:


##Convert Matrix of Laurent polynomials into Vector of Matrix elements
##Output: Display mask/filter in Vector/Array of matrix elements u(k)
##Input: mpoly, matrix of Laurent polynomials for a mask/filter
D1PrintMask:=proc(mpoly)
    local j,k,n,rowdim,coldim,supp,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    supp:=D1SupportPoly(mpoly,z):
    result:=Vector[row](supp[2]-supp[1]+1,0):
    for j from supp[1] to supp[2] do
        result[j-supp[1]+1]:=Matrix(rowdim,coldim,0):
    end do:
    for j from 1 to rowdim do
    for k from 1 to coldim do
        for n from supp[1] to supp[2] do
            result[n-supp[1]+1][j,k]:=simplify(coeff(mpoly[j,k],z,n)):
        end do:
    end do: end do:
    print(convert(result,matrix)):
    print("The matrix-valued mask/filter is supported on",supp):
    result:=result:
end proc:


##Generalize conjugate for complex numbers to complex expressions
##except I, all other involved parameters are real valued
ConjComplex:=proc(cnum)
    local tmp,result:
    tmp:=SimplifyComplex(cnum): #convert it into real+I*imaginar
    result:=coeff(expand(tmp),I,0)-I*coeff(expand(tmp),I,1):
    result:=result:
end proc:

##conjugate&transpose (choice<>1,2,3), conjugate coeff only (choice=1), conjugate and z->1/z (choice=2), conjugate and transpose coeff only (choice=2) no z->1/z for coefficients of a matrix polynomial
##Explanation: For mpoly:=\sum_{k\in\Z} u(k)var^k, then output
##choice=1,\sum_{k\in\Z}\overline{u(k)} var^k
##choice=2,\sum_{k\in\Z}\overline{u(k)} var^(-k)
##choice=3,\sum_{k\in\Z}\overline{u(k)}^T var^k
##otherwise,mpoly^\star(var):=\sum_{k\in\Z}\overline{u(k)}^T var^(-k)
D1ConjTranPoly:=proc(mpoly,var,choice)
    local j,k,n,supp,rowdim,coldim,isrealcf,tmppoly,tmp,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    supp:=D1SupportPoly(mpoly,var):
    isrealcf:=D1IsRealFilter(mpoly):
    if D1IsRealFilter(mpoly) then
        tmp:=mpoly:
    else  #take complex conjugate on coefficients
        tmp:=Matrix(rowdim,coldim,0):
        for j from 1 to rowdim do
        for k from 1 to coldim do
            tmppoly:=expand(mpoly[j,k]):
            for n from supp[1] to supp[2] do #conjugate coeff only
                tmp[j,k]:=tmp[j,k]+ConjComplex(coeff(tmppoly,var,n))*var^(n):
            end do:
        end do: end do:
    end if:
    if (choice=1) then #conjugate coeffs only
        result:=tmp:
    elif (choice=2) then #conjugate both coeffs and var:z->1/z.
        result:=Subs({var=1/var},tmp):
    elif (choice=3) then #conjugate and transpose coeffs.
        result:=Transpose(tmp):
    else  #conjugate, transpose both coeffs and var:(z->1/z)
        result:=Transpose(Subs({var=1/var},tmp)):
    end if:
    result:=result:
end proc:


##Calculate the coset sequence of a matrix poly
##choice=0 for \sum_{k\in\Z} a(\gamma+dilation*k)z^k for standard
##coset at gamma: a^{[gamma]}(k):= a(gamma+dilation*k), k\in Z
##choice=otherwise, output
##\sum_{k\in\Z} a(\gamma+dilation*k)z^(gamma+dilation*k)
##Remark: the choice<>0 sometimes is convenient.
D1CosetSeq:=proc(mpoly,dilation,gamma,choice)
    local j,k,n,nn,rowdim,coldim,supp,result:
#    lprint("choice=0 for sum_{k in Z}a(gamma+dilation*k)z^k; choice=otherwise for sum_{k in Z} a(gamma+dilation*k)z^(gamma+dilation*k)."):
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    supp:=D1SupportPoly(mpoly,z):
    result:=Matrix(rowdim,coldim,0):
    for n from supp[1] to supp[2] do
        nn:=trunc((n-gamma)/dilation):
        if (n=(nn*dilation+gamma)) then
            if (choice<>0) then nn:=n: end if:
            for j from 1 to rowdim do
            for k from 1 to coldim do
                result[j,k]:=result[j,k]+coeff(mpoly[j,k],z,n)*z^nn:
            end do: end do:
        end if:
    end do:
    result:=result:
end proc:


##Calculate Taylor polynomials to degree degord of a matrix poly
##mpoly(exp(-I*(xi+pt)) at the base point \xi=0
##output a Taylor polynomials of power xi.
D1TaylorPoly:=proc(mpoly,pt,degord)
    local j,k,n,rowdim,coldim,tmppoly,result:
    rowdim:=RowDimension(mpoly):
    coldim:=ColumnDimension(mpoly):
    result:=Matrix(rowdim,coldim,0):
    for j from 1 to rowdim do
    for k from 1 to coldim do
        tmppoly:=expand(subs({z=exp(-I*(xi+pt))},mpoly[j,k])):
        result[j,k]:=convert(taylor(tmppoly,xi=0,degord),polynom):
    end do: end do:
    result:=result:
end proc:


##Output: set of equations for Hermite interpolatory filters
##See (6.2.5) in page of the book.
D1HermiteCond:=proc(mpoly,dilation)
    local j,multi,tmppoly,result:
    tmppoly:=D1CosetSeq(mpoly,dilation,0,0):
    multi:=RowDimension(mpoly):
    for j from 1 to multi do
          tmppoly[j,j]:=tmppoly[j,j]-dilation^(-j):
    end do:
    result:=D1PolyToCoeffSet(tmppoly,z):
    result:=result:
end proc:


##Output: vgu (matching filter) for Hermite interpolatory filters
##See (6.2.6) in page 503 and (6.2.8) in page 505 of the book.
D1Hermitevgu:=proc(hmord)
    local j,result:
    result:=Vector(hmord,0):
    result[1]:=1:
    for j from 2 to hmord do
        result[j]:=(I*xi)^(j-1):
    end do:
    result:=result:
end proc:


##Parametrize a Matrix of Laurent polynomials supported on [low,high] with real coefficients and variable z
##ouput: result[1]=matrix polynomials
##       result[2]=set of parameters/unknowns in result[1]
D1ParaMaskVar:=proc(rowdim,coldim,low,high,para)
    local j,k,n,ind,tmp,result:
    result:=Vector(2):
    tmp:=Matrix(rowdim,coldim,0):
    result[2]:={}:
    ind:=1:
    for j from 1 to rowdim do
    for k from 1 to coldim do
        for n from low to high do
            tmp[j,k]:=tmp[j,k]+evaln(para||ind)*z^n:
            result[2]:=result[2] union {evaln(para||ind)}:
            ind:=ind+1:
        end do:
    end do: end do:
    result[1]:=tmp:
    result:=result:
end proc:


##Parametrize a Matrix of Laurent polynomials supported on [low,high] with real coefficients and variable z
##ouput: matrix polynomials
D1ParaMask:=proc(rowdim,coldim,low,high,para)
    local tmp,result:
    tmp:=D1ParaMaskVar(rowdim,coldim,low,high,para):
    result:=tmp[1]:
    result:=result:
end proc:


##Parametrize a Matrix of polynomials in variable xi with  O(|\xi|^degord)
##ouput: result[1]=matrix polynomials
##       result[2]=set of parameters/unknowns in result[1]
D1ParaPolyVar:=proc(rowdim,coldim,degord,para)
    local j,k,n,ind,result:
    result:=Vector(2):
    result[1]:=Matrix(rowdim,coldim,0):
    result[2]:={}:
    ind:=1:
    for j from 1 to rowdim do
    for k from 1 to coldim do
        for n from 0 to (degord-1) do
            result[1][j,k]:=result[1][j,k]+evaln(para||ind)*(I*xi)^n:
            result[2]:=result[2] union {evaln(para||ind)}:
            ind:=ind+1:
        end do:
    end do: end do:
    result:=result:
end proc:


##Parametrize a Matrix of polynomials in variable xi with  O(|\xi|^degord)
D1ParaPoly:=proc(rowdim,coldim,degord,para)
    local tmp,result:
    tmp:=D1ParaPolyVar(rowdim,coldim,degord,para):
    result:=tmp[1]:
    result:=result:
end proc:


##Output: set of equations to determine \hat{\phi}(\xi)+O(|\xi|^degord) at \xi=0 such that
##\hat{\phi}(dilation\xi)=apoly(\exp(-I*\xi)\hat{\phi}(\xi)+O(|\xi|^degord) at \xi=0
D1phiMomentEqs:=proc(apoly,hatphimt,dilation,degord)
    local tmp,result:
    tmp:=Subs({xi=dilation*xi},hatphimt):
    tmp:=tmp-D1TaylorPoly(apoly,0,degord).hatphimt:
    result:=D1BigOEqs(tmp,xi,degord,-1):
    result:=result:
end proc:


##Input: mpoly, a Matrix(column Vector) of polynomials in variable xi
##Output: normalize mpoly so that its first nonzero entry is 1.
##If the first column is 0, do nothing and return mpoly
D1NormalizePoly:=proc(mpoly)
    local j,rowdim,tmp,result:
    rowdim:=RowDimension(mpoly):
    tmp:=Subs({xi=0},mpoly):
    for j from 1 to rowdim do
        if (simplify(tmp[j,1])<>0) then
            result:=SolveEqs({simplify(tmp[j,1])-1},mpoly):
        end if:
    end do:
    result:=result:
end proc:


##Output: set of equations to determine \hat{\vgo}(\xi)+O(|\xi|^degord) at \xi=0 such that
#\hat{\vgu}(dilation\xi)apoly(\exp(-I*\xi)=\hat{\vgu}(\xi)+O(|\xi|^degord) at \xi=0
##For explanation, see (5.6.5) in page 418 of the book
D1vguMomentEqs:=proc(apoly,dilation,vgu,degord)
    local tmp,result:
    tmp:=Subs({xi=dilation*xi},vgu):
    tmp:=tmp.D1TaylorPoly(apoly,0,degord)-vgu:
    result:=D1BigOEqs(tmp,xi,degord,-1):
    result:=result:
end proc:


##Output:set of equations for prescribed symmetry of matrix filters
##Matrix-valued mpoly can be any size of matrix of Laurent polys.
##Symmetry:\hat{u}(\xi)=symb(dilation*xi)\hat{u}(-\xi)syma(\xi)^{-1}
#for simplicity, syma/symb are Vectors (instead diagonal matrices)
##For symmetry of filters, see Exercise 6.42(items(c), (g)) on pages 575-576 of the book.
##Remark: only for symmetry, cannot handle complex symmetry
##Use the routine D1SymMaskEqsC(apoly,dilation,syma,symb)
##instead for complex symmetry
D1SymMaskEqs:=proc(apoly,dilation,syma,symb)
    local j,rowdim,coldim,symleft,symright,tmp,result:
    rowdim:=RowDimension(apoly):
    coldim:=ColumnDimension(apoly):
    symleft:=Matrix(rowdim,rowdim,0):
    symright:=Matrix(coldim,coldim,0):
    for j from 1 to rowdim do
        symleft[j,j]:=subs({z=z^dilation},symb[j]):
    end do:
    for j from 1 to coldim do
        symright[j,j]:=subs({z=1/z},syma[j]):
    end do:
    tmp:=Subs({z=1/z},apoly):
    tmp:=apoly-symleft.tmp.symright:
    result:=D1PolyToCoeffSet(tmp,z):
    result:=result:
end proc:


##Output: Vector of sets of equations for sum rules of a matrix-valued filter: result[0] consists of all terms from 0 to srord-1.
##result[j] consists of only terms for \xi^j.
##Matrix-valued poly must be a square matrix of Laurent polys.
##For scalar filters, row poly vector vgu of xi is not needed.
##srord is the order of sum rules.
##For definition of sum rules, see (5.6.28) in page 425 and (5.5.13) in page 407 of the book.
##Remark: much faster by using cosets of filters/masks
D1SumRuleEqs:=proc(apoly,dilation,vgu,srord)
    local j,k,multi,tmp,mvgu,mvgu2,cf,result:
    multi:=RowDimension(apoly):
    result:=Array(0..srord):
    if (multi=1) then #special case for a scalar filter
        cf:=expand(subs({z=1},apoly[1,1])-1):
        result[0]:={coeff(cf,I,0),coeff(cf,I,1)}:
        result[1]:=result[0]:
        for k from 2 to srord do:
            result[k]:={}:
        end do:
        mvgu:=D1CosetSeq(apoly,dilation,0,1): #0-coset
        for j from 1 to (abs(dilation)-1) do
            tmp:=D1TaylorPoly(D1CosetSeq(apoly,dilation,j,1)-mvgu,0,srord):
            result[0]:=result[0] union D1BigOEqs(tmp,xi,srord,-1):
            for k from 1 to srord do
                result[k]:=result[k] union D1BigOEqs(tmp,xi,srord,k-1):
            end do:
        end do:
    else
        mvgu:=Matrix(1,multi,0):
        mvgu[1,1..multi]:=vgu:
        mvgu2:=Subs({xi=dilation*xi},mvgu):
        for k from 0 to srord do
            result[k]:={}:
        end do:
        for j from 0 to (abs(dilation)-1) do
            tmp:=abs(dilation)*mvgu2.D1TaylorPoly(D1CosetSeq(apoly,dilation,j,1),0,srord)-mvgu:
            result[0]:=result[0] union D1BigOEqs(tmp,xi,srord,-1):
            for k from 1 to srord do
                result[k]:=result[k] union D1BigOEqs(tmp,xi,srord,k-1):
            end do:
        end do:
    end if:
    result:=result:
end proc:

##Output: a matrix-valued filter satisfying sum rules
D1SumRulePoly:=proc(apoly,dilation,vgu,srord)
    local j,mvgu2,tmp,result:
    tmp:=D1SumRuleEqs(apoly,dilation,vgu,srord):
    result:=SolveEqs(tmp[0],apoly):
    result:=result:
end proc:

##Find the sum rules order of a given filter/mask
##maxsr is a guessed the largest possible sum rules
##Output: largest sum rule order(<=maxsr) satisfied by mpoly
##For definition of sum rules, see (5.6.28) in page 425 and (5.5.13) in page 407 of the book.
D1FindMaxSumRule:=proc(apoly,dilation,maxsr)
    local j,EQ,tmp,solt,vgu,result:
    vgu:=D1vguMoment(apoly,dilation,maxsr):
    tmp:=D1SumRuleEqs(apoly,dilation,vgu,maxsr):
    ##check when we have a solution
    EQ:={}:
    result:=0:
    for j from 1 to maxsr do
        EQ:=EQ union tmp[j]:
        solt:=solve(EQ):
        if (solt=NULL) then
            break;
        else
            result:=j:
        end if:
    end do:
    result:=result:
end proc:


##Output: set of equations for biorthogonal or perpendicular filters
##choice=0, biorthogonal filters with real coefficients:
##\sum_{j=0}^{dilation-1}\hat{\tilde{a}}(\xi+2\pi j/dilation)\overline{\hat{a}(\xi+2\pi j/dilation)=I_multi
##choice=otherwise, perpendicular filters with real coefficients: ##\sum_{j=0}^{dilation-1}\hat{\tilde{a}}(\xi+2\pi j/dilation)\overline{\hat{a}(\xi+2\pi j/dilation)=0
##Note: dpoly=\tilde{a} which is biorthogonal to the filter a
##For definition of biorthogonal/perpendicular, see (4.5.5) in page 309, (6.4.14) in page 527, or (6.4.21) in page 531 of the book.
D1BiorthEqs:=proc(apoly,dpoly,dilation,choice)
    local j,multi,tmppoly,result:
    multi:=RowDimension(apoly):
    tmppoly:=Matrix(multi,multi,0):
    if (choice=0) then
        for j from 1 to multi do #for biorthogonal filters
            tmppoly[j,j]:=1/abs(dilation):
        end do:
    end if:
    tmppoly:=tmppoly-dpoly.D1ConjTranPoly(apoly,z,0):
    tmppoly:=D1CosetSeq(tmppoly,dilation,0,0):
    result:=D1PolyToCoeffSet(tmppoly,z):
    result:=result:
end proc:



##Output: a dual filter for a given filter with prescribed sum rules
##and with prescribed symmetry
##This routine only handles real-valued filters
D1DualMaskSym:=proc(apoly,dilation,low,high,dsrord,syma,para)
    local j,multi,EQ,dualvgu,tmppoly,tmp,result:
    multi:=RowDimension(apoly):
    tmp:=D1ParaMaskSym(multi,multi,low,high,para,syma,syma):
    dualvgu:=D1Dualvgu(apoly,dilation,dsrord):
    tmppoly:=D1SumRulePoly(tmp,dilation,dualvgu,dsrord):
    EQ:=D1BiorthEqs(apoly,tmppoly,dilation,0):
    result:=SolveEqs(EQ,tmppoly):
    result:=result:
end proc:



##Plot a wavelet vector function \psi satisfying
##\psi(x)=|dilation| \sum_{k\in \Z} b(k) \phi(dilation*x-k)
##Or plot a refinable vector function \phi satisfying
##\phi(x)=|dilation| \sum_{k\in\Z} a(k) \phi(dilation*x-k)
##with the normalization condition: the first nonzero entry in \wh{\phi}(0) is 1. This is the special case for psi with b=a
##reslevel, the data is on lattice: dilation^(-reslevel)*Z
##Ouput: multi+1 Vectors of data xval, psi/phi_1,...,psi/phi_multi
##We use abs(dilation)^n*[b(2^{n-1}\cdot)*a_{n-1}*\wh{\phi}(0)](k) located at dilation^{-n}k, n=reslevel
##where the Fourier series of b(2^{n-1}\cdot)*a_{n-1} is ##\hat{b}(2^{n-1}\xi)\hat{a_{n-1}}(\xi)=\hat{b}(2^{n-1}\xi)
##\hat{a}(dilation^{n-2}\xi)\cdot\hat{a}(dilation*\xi)\hat{a}(\xi)
##For boundary, we pad zero at both endpoints for better effect.
##Note that supp(phi)=fsupp(a)/(dilation-1)
##and supp(psi)=[supp(phi)+supp(b)]/dilation
##This routine handles complex-valued filters/phi as well
##Use bpoly=apoly to plot a refinable vector function
##see Subsection 6.2.3 (with m=0) in page 511 of the book for detail
D1Plotphi:=proc(apoly,bpoly,dilation,reslevel)
    local j,k,multi,colb,xval,absxval,suppa,suppb,supp,mpolyn,hatphi0,tmp,result:
    multi:=RowDimension(apoly): #multi*multi
    colb:=RowDimension(bpoly):
    print("Column Vectors of xval, psi/phi_1,...psi/phi_s for plotting psi/phi normalizing the first nonzero entry of \hat{\phi}(0) to be 1"):
    suppa:=D1SupportPoly(apoly,z):
    suppb:=D1SupportPoly(bpoly,z):
    mpolyn:=bpoly:
    hatphi0:=Subs({xi=0},D1phiMoment(apoly,dilation,1)):
    for j from 1 to (reslevel-1) do
        mpolyn:=Subs({z=z^dilation},mpolyn).apoly:
    end do:
    xval:=dilation^(-reslevel):
    absxval:=1/abs(xval):
    mpolyn:=D1SimplifyPoly(absxval*mpolyn.hatphi0,16):
    supp:=D1SupportPoly(mpolyn,z):
    result:=Array(0..colb):
    for j from 0 to colb do
        result[j]:=Vector(supp[2]-supp[1]+3,0):
    end do:
    if (dilation>0) then #modify boundary for better looking.
        result[0][1]:=evalf((suppa[1]/(dilation-1)+suppb[1])/dilation,16):
        result[0][supp[2]-supp[1]+3]:=evalf((suppa[2]/(dilation-1)+suppb[2])/dilation,16):
    else
        result[0][1]:=evalf((suppa[2]/(dilation-1)+suppb[2])/dilation,16):
        result[0][supp[2]-supp[1]+3]:=evalf((suppa[1]/(dilation-1)+suppb[1])/dilation,16):
    end if:
    for j from 1 to multi do
        result[j][1]:=0:  result[j][supp[2]-supp[1]+3]:=0:
    end do:
    for k from supp[1] to supp[2] do
        result[0][k-supp[1]+2]:=evalf(xval*k,16):
        tmp:=D1MaskElement(mpolyn,k,z):
        for j from 1 to multi do
            result[j][k-supp[1]+2]:=tmp[j,1]:
        end do:
    end do:
    result:=result:
end proc:


##plot orthogonal phi with normalization \|\hat{\phi}(0)\|_2=1
##This routine handles complex-valued filters/phi as well
D1PlotOrthphi:=proc(apoly,dilation,reslevel)
    local j,multi,hatphi0,tmp,result:
    multi:=RowDimension(apoly):
    print("Column Vectors of xval, phi_1,...phi_r for plotting phi with normalization \|\hat{\phi}(0)\|_2=1."):
    hatphi0:=Subs({xi=0},D1phiMoment(mpoly,dilation,1)):
    result:=D1Plotphi(apoly,apoly,dilation,reslevel):
    tmp:=0:
    for j from 1 to multi do
        tmp:=tmp+(abs(hatphi0[j]))^2:
    end do:
    for j from 0 to multi do
        result[j]:=evalf(1/sqrt(tmp),16)*result[j]:
    end do:
    result:=result:
end proc:

##Ouput: L2 smoothness exponent sm_2(a,dilation), the L2 of filter mpoly
##Output: Array of All information: result[0]=sm_2(a,dilation)
##result[1]=Maximum sum rules, which is used for computing
##result[2]=Vector of the eigenvalues of \hat{a}(0) in modulus
##result[3]=Vector of the two different eigenvalues in modulus
##result[4]=Vector of the special eigenvalues in modulus
##result[5]=Vector of all the eigenvalues Tmatrix in modulus
##result[6]=Matrices TMatrix to generate all the eigenvalues
##For algorithm, see theorem 5.8.4 in page 461 of the book
D1L2SmoothnessPoly:=proc(apoly,dilation,maxsr,errtol)
    local j,k,J,K,multi,sr,supp,len,lenc,convpoly,clist,Tmatr,eigTmatr,eigA,eigsp,result:
    result:=Array(0..6):
    multi:=RowDimension(apoly):
    sr:=D1FindMaxSumRule(apoly,dilation,maxsr):
    eigsp:=Vector(2*sr*multi,0):
    ##generate the Vector of special eigenvalues
    eigA:=sort(abs(Eigenvalues(evalf(Subs({z=1},apoly),16))),`>`):
    result[1]:=sr:
    result[2]:=Transpose(eigA):
    #print("The first element of sorted eigenvalues of hat{a}(0) must be 1!",eigA):
    for j from 1 to (2*sr) do
        eigsp[j]:=dilation^(1-j):
    end do:
    J:=2*sr:
    for k from 2 to multi do
    for j from 0 to (sr-1) do
        eigsp[J]:=dilation^(-j)*eigA[k]:
        eigsp[J+1]:=dilation^(-j)*eigA[k]:
        J:=J+2:
    end do: end do:
    sr:=J-1:
    eigsp:=sort(evalf(eigsp,16),`>`):
    result[4]:=Transpose(eigsp):
    supp:=D1SupportPoly(apoly,z):
    len:=supp[2]-supp[1]:
    convpoly:=D1ConjTranPoly(apoly,z,2): #conjugate both coeff and variable
    convpoly:=KroneckerProduct(convpoly,apoly):
    supp:=D1SupportPoly(convpoly,z):
    clist:=Array(supp[1]..supp[2]):
    for j from supp[1] to supp[2] do
        clist[j]:=Transpose(D1MaskElement(convpoly,j,z)):
    end do:
    multi:=multi^2:
    Tmatr:=Matrix(multi*(2*len+1),multi*(2*len+1),0):
    for j from -len to len do
    for k from -len to len do
        if (abs(j-dilation*k)<=supp[2]) then
            J:=multi*(j+len): K:=multi*(k+len):
            Tmatr[(J+1)..(J+multi),(K+1)..(K+multi)]:=abs(dilation)*clist[j-dilation*k]:
        end if:
    end do: end do:
    eigTmatr:=sort(abs(Eigenvalues(evalf(Tmatr,16))),`>`):
    result[5]:=Transpose(eigTmatr):
    result[6]:=evalf(Tmatr,16):
    ##Now find the first different element between two Vectors
    result[0]:=-1:
    for j from 1 to sr do
        if (evalf(abs(eigTmatr[j]-eigsp[j]),32)>errtol) then
            result[0]:=eigTmatr[j]:
            result[3]:=Vector[row](2, [eigsp[j],result[0]]):
            break;
        end if:
    end do:
    #To handle: all the special eigenvalues belong to eigTmarix
    if (j=(sr+1)) then
        result[0]:=eigTmatr[j]:
        result[3]:=Vector[row](2, [infinity,result[0]]):
    end if: 
    if (result[0]<=0) then
        print("Something is wrong when calculating sm_2(a,dilation)."):
    else
        result[0]:=evalf(-0.5*ln(result[0])/ln(abs(dilation)),16):
    end if:
    print("the L2 smoothness exponent sm_2(filter,dilation) is",result[0]):
    result:=result:
end proc:


##Output: a matrix of of product of each pair of filters
##that is, for a pair of given filter tupoly and upoly with
##multiplicity multi, output a (2multi)*(2multi) matrix:
##[\overline{\hat{tupoly^{[0]}}(\xi)}]
##[\overline{\hat{tupoly^{[1]}}(\xi)}]
##...
##[\overline{\hat{tupoly^{[abs(dilation)-1]}}(\xi)}]
##multiply [\hat{upoly^{[0]}}(\xi),...,\hat{upoly^{[abs(dilatioon)-1]}}(\xi)]
##For information, see (6.4.5) in page 523 of the book (need to rewrite it into coset form
D1FrameletEqs:=proc(mpoly,dmpoly,dilation)
    local gamma,multi,tmp,tmppoly,tdmpoly,result:
    multi:=RowDimension(mpoly):
    tmppoly:=Matrix(multi,multi*abs(dilation),0):
    tdmpoly:=Matrix(multi,multi*abs(dilation),0):
    for gamma from 0 to (abs(dilation)-1) do
        tmp:=D1CosetSeq(mpoly,dilation,gamma,0):
        tmppoly[1..multi,(gamma*multi+1)..(gamma+1)*multi]:=tmp:
        tmp:=D1CosetSeq(dmpoly,dilation,gamma,0):
        tdmpoly[1..multi,(gamma*multi+1)..(gamma+1)*multi]:=tmp:
    end do:
    result:=D1ConjTranPoly(tdmpoly,z,0).tmppoly:
    result:=result:
end proc:


##Generate Matrix representation for gamma-shifted transition operator
D1TransitOp:=proc(mpoly,dilation,gamma)
    local j,k,J,K,multi,supp,len,clist,result:
    multi:=RowDimension(mpoly):
    supp:=D1SupportPoly(mpoly,z):
    len:=supp[2]-supp[1]+1:
    clist:=Array(supp[1]..supp[2]):
    for j from supp[1] to supp[2] do
        clist[j]:=D1ConjTranPoly(D1MaskElement(mpoly,j,z),z,0):
    end do:
    result:=Matrix(multi*len,multi*len,0):
    for j from -supp[2] to -supp[1] do
    for k from -supp[2] to -supp[1] do
        if (((gamma+j-dilation*k)>=supp[1]) and ((gamma+j-dilation*k)<=supp[2])) then
            J:=multi*(j+supp[2]): K:=multi*(k+supp[2]):
            result[(J+1)..(J+multi),(K+1)..(K+multi)]:=abs(dilation)*clist[gamma+j-dilation*k]:
        end if:
    end do: end do:
    result:=result:
end proc:


##Find out all the free parameter in a solution set
##The criterion used here is to find left=right, that is the left side
##and the right side are the same
##choice=0, show free parameters, otherwise not
##otherwise, replace the free parameters with the ordered parameters:
##c.startvalue,c.(startvalue+1),...
##So in set solt, the parameter should NOT use c.number
##Output: a list of all variables in the set such that c=c.
##        that is, there are free parameters.
##Output all the free parameters
FreeParameter:=proc(solt,newpara,startvalue,choice)
    local j,ind,cfpara,tmp,var,result:
    result:={}: tmp:={}: var:={}: ind:=0:
    if (solt={}) then
        printf("It is the empty set!\n"):
    else
        for j from 1 to nops(solt) do
            cfpara:=op(1,solt[j]):
            if(cfpara=op(2,solt[j])) then #free parameter
                tmp:=tmp union {cfpara=evaln(newpara||(ind+startvalue))}:
                var:=var union {evaln(newpara||(ind+startvalue))}:
                ind:=ind+1:
            else
                result:=result union {solt[j]}:
            end if:
        end do:
    end if:
    if (choice=0) then
        printf("\n Total %d free parameters given by:\n", ind):
        print(var):
    end if:
    result:=subs(tmp,result):
    result:=result union tmp:
    result:=result:
end proc:


##solve a set of equations, if no solutions, give warning and return the same filters;
##if there are solutions, substitute the first solution to a poly (Matrix)
##choice=0, don't show free parameters, otherwise show it
SolveEqsVar:=proc(seteqns,var,newvar,mpoly,choice)
    local solt,result:
    solt:=solve(seteqns,var):
    if (solt=NULL) then
        print("In SolveEqsVar(seteqns,var,newvar,mpoly,choice): No solutions to the given set of equations in seteqns!"):
        result:=mpoly:
    else
        solt:=FreeParameter(solt,newvar,1,choice):
        result:=D1SimplifyPoly(Subs(solt,mpoly),-1):
    end if:
    result:=result:
end proc:



##Output: Matrix(column Vector) of Taylor polynomials of \wh{\phi}(\xi)+O(|\xi|^degord) at \xi=0
##such that the first nonzero entry of \hat{\phi}(0) is 1 and
##\hat{\phi}(dilation\xi)=apoly(\exp(-I*\xi)\hat{\phi}(\xi)+O(|\xi|^{degord}) at \xi=0
##if mpoly has real coefficients, then output real phimoments
D1phiMoment:=proc(apoly,dilation,degord)
    local coldim,hatphi,tmp,result:
    coldim:=ColumnDimension(apoly):
    if D1IsRealFilter(apoly) then
        hatphi:=D1ParaPolyVar(coldim,1,degord,hphi):
    else
        hatphi:=D1ParaPolyVarC(coldim,1,degord,hphi):
    end if:
    tmp:=D1phiMomentEqs(apoly,hatphi[1],dilation,degord):
    result:=D1NormalizePoly(SolveEqsVar(tmp,hatphi[2],phivar,hatphi[1],1)):
    print("With normalization \hat{\phi}(0)=",D1SimplifyPoly(Subs({xi=0},result),0)):
    result:=D1SimplifyComplexPoly(D1SimplifyPoly(result,0),xi):
end proc:


##Output: Matrix (row vector) Taylor polynomials of \hat{\vgu}(\xi)+O(|\xi|^\degord) at \xi=0 such that
##the first nonzero entry of \hat{\vgu}(0) is 1 and
##\hat{\vgu}(dilation\xi)apoly(\exp(-I*\xi)=\hat{\vgu}(\xi)+O(|\xi|^degord) at 0
##For explanation, see (5.6.5) in page 418 of the book
D1vguMoment:=proc(apoly,dilation,degord)
    local rowdim,vgu,tmp,result:
    rowdim:=RowDimension(apoly):
    if D1IsRealFilter(apoly) then
        vgu:=D1ParaPolyVar(1,rowdim,degord,vvar):
    else
        vgu:=D1ParaPolyVarC(1,rowdim,degord,vvar):
    end if:
    tmp:=D1vguMomentEqs(apoly,dilation,vgu[1],degord):
    result:=Transpose(D1NormalizePoly(Transpose(SolveEqsVar(tmp,vgu[2],vguvar,vgu[1],1)))):
#    print("With normalization hat{vgu}(0)=",Subs({xi=0},result)):
    result:=result:
end proc:


##Output: dualvgu(matching filter) for a dual/biorthogonal filter of a given filter
##such that the first nonzero entry in dualvgu(\xi=0) is 1 and satisfies
##\hat{dualvgu}(\xi)=\overline{\hat{a}(\xi)}^T=\hat{\dualvgu}(dilation*\xi)+O(|\xi|^dualvguord),\xi\to 0
##See (6.5.1) in page 537 of the book.
D1Dualvgu:=proc(apoly,dilation,dualvguord)
    local tmp,result:
    tmp:=D1phiMoment(apoly,dilation,dualvguord):
    result:=D1ConjTranPoly(tmp,xi,3): #conjugate and transpose coeffs
#    print("The dualvgu is normalized to be \overline{\hat{\phi}(0)}=",Subs({xi=0},result)):
    result:=result:
end proc:


##Output: symmetric vgu for apoly with symmetry
##Symmetry:\hat{\vgu}(\xi)=\hat{\vgu}(-xi)syma(\xi)^{-1}+O(|xi|^vguord)
##See Exercise 6.43 on page 576.
D1vguSym:=proc(apoly,syma,vguord,vguvar)
    local j,multi,symright,vgu,tmp,result:
    multi:=RowDimension(apoly):
    vgu:=D1ParaPolyVar(1,multi,vguord,tvguvar):
    symright:=Matrix(multi,multi,0):
    for j from 1 to multi do
        symright[j,j]:=syma[j]:
    end do:
    symright:=D1TaylorPoly(Subs({z=1/z},symright),0,vguord):
    tmp:=Subs({xi=-xi},vgu[1]):
    tmp:=vgu[1]-tmp.symright:
    result:=SolveEqsVar(D1BigOEqs(tmp,xi,vguord,-1),vgu[2],vguvar,vgu[1],0):
    result:=result:
end proc:



##Output: symmetric vgu for apoly with symmetry
##Symmetry:\hat{\vgu}(\xi)=\hat{\vgu}(-xi)syma(\xi)^{-1}+O(|xi|^vguord)
##Complex symmetry: \hat{\vgu}(\xi)=\overline{\hat{\vgu}(xi)}syma(\xi)^{-1}+O(|xi|^vguord)
##See Exercise 6.43 on page 576.
##choice<0 for complex symmetry; and choice>=0 for symmetry
D1vguSymC:=proc(apoly,syma,vguord,vguvar,choice)
    local varr,vari,apolyri,result;
    apolyri:=D1RealImagPoly(apoly):
    varr:=evaln(vguvar||r):
    vari:=evaln(vguvar||i):
    result:=D1vguSym(apolyri[1],syma,vguord,varr):
    if (choice<0) then
        result:=result+I*D1vguSym(apolyri[2],syma,vguord,vari):
    else
        result:=result+I*D1vguSym(apolyri[2],-syma,vguord,vari):
    end if:
    result:=result:
end proc:


##Parametrize a Matrix of Laurent polynomials in z with symmetry
D1ParaMaskSym:=proc(rowdim,coldim,low,high,para,syma,symb)
    local tmp,result:
    tmp:=D1ParaMaskVar(rowdim,coldim,low,high,tpara):
    result:=SolveEqsVar(D1SymMaskEqs(tmp[1],dilation,syma,symb),tmp[2],para,tmp[1],1):
    result:=result:
end proc:


##Parametrize a matrix-valued Hermite interpolatory filters
##supported on [low,high] with variable z
D1ParaHermiteMask:=proc(hmord,dilation,low,high,para)
    local tmp,tpara,result;
    tpara:=evaln(t||para):
    tmp:=D1ParaMaskVar(hmord,hmord,low,high,tpara):
    result:=SolveEqsVar(D1HermiteCond(tmp[1],dilation),tmp[2],para,tmp[1],0):
    result:=result:
end proc:


##Parametric a matrix-valued Hermite interpolatory filters
##with symmetry and supported on [low,high] with variable z
D1ParaHermiteMaskSym:=proc(hmord,dilation,low,high,syma,para)
    local tmp,EQ,tpara,result:
    tpara:=evaln(t||para):
    tmp:=D1ParaMaskVar(hmord,hmord,low,high,tpara):
    EQ:=D1HermiteCond(tmp[1],dilation):
    EQ:=EQ union D1SymMaskEqs(tmp[1],dilation,syma,syma):
    result:=SolveEqsVar(EQ,tmp[2],para,tmp[1],0): #show free parameters
    result:=result:
end proc:


##Output: a dual filter for a given filter with prescribed sum rules
##but no symmetry is imposed
D1DualMask:=proc(apoly,dilation,low,high,dsrord,para)
    local j,multi,dualvgu,tmppoly,tmp,result:
    multi:=RowDimension(apoly):
    if D1IsRealFilter(apoly) then
        tmp:=D1ParaMaskVar(multi,multi,low,high,para):
    else
        tmp:=D1ParaMaskVarC(multi,multi,low,high,para):
    end if:
    #add sum rule conditions
    dualvgu:=D1Dualvgu(apoly,dilation,dsrord):
    tmppoly:=D1SumRulePoly(tmp[1],dilation,dualvgu,dsrord):
    #add biorthogonal conditions
    result:=SolveEqs(D1BiorthEqs(apoly,tmp[1],dilation,0),tmp[1],1): #show free parameters
    result:=result:
end proc:


##Output: a dual filter for a given filter with prescribed sum rules
##but no symmetry is imposed
##for complex-valued dual filters
D1DualMaskC:=proc(apoly,dilation,low,high,dsrord,para)
    local j,multi,dualvgu,tmppoly,result:
    multi:=RowDimension(apoly):
    tmppoly:=D1ParaMaskC(multi,multi,low,high,para):
    #add sum rule conditions
    dualvgu:=D1Dualvgu(apoly,dilation,dsrord):
    tmppoly:=D1SumRulePoly(tmppoly,dilation,dualvgu,dsrord):
    #add biorthogonal conditions
    result:=SolveEqs(D1BiorthEqs(apoly,tmppoly,dilation,0),tmppoly):
    result:=result:
end proc:



####Old and not efficient routines####

##Output: set of equations for sum rules of a matrix-valued filter
##Matrix-valued poly must be a square matrix of Laurent polys.
##For scalar filters, row poly vector vgu of xi is not needed.
##srord is the order of sum rules.
##For definition of sum rules, see (5.6.28) in page 425 of the book.
##Remark: this implementation is not very efficient
D1SumRuleEqsOriginal:=proc(mpoly,dilation,vgu,srord)
    local j,multi,tmp,mvgu,mvgu2,result:
    multi:=RowDimension(mpoly):
    if (multi=1) then #special case for a scalar filter
        result:={subs({z=1},mpoly[1,1])-1}:
        for j from 1 to (abs(dilation)-1) do
            tmp:=D1TaylorPoly(mpoly,2*Pi*j/dilation,srord):
            result:=result union D1BigOEqs(tmp,xi,srord-1,-1)
        end do:
    else
        mvgu:=Matrix(1,multi,0):
        mvgu[1,1..multi]:=vgu:
        mvgu2:=Subs({xi=dilation*xi},mvgu):
        tmp:=mvgu-mvgu2.D1TaylorPoly(mpoly,0,srord):
        result:=D1BigOEqs(tmp,xi,srord-1,-1):
        for j from 1 to (abs(dilation)-1) do
            tmp:=mvgu2.D1TaylorPoly(mpoly,2*Pi*j/dilation,srord):
            result:=result union D1BigOEqs(tmp,xi,srord-1,-1):
        end do:
    end if:
    result:=result:
end proc:



#############################################################
##Routines ending with C for complex filters with
##complex-valued unknowns/parameters
#############################################################

##Parametrize a Matrix of Laurent polynomials supported on [low,high] with variable z and complex coefficients
##ouput: result[1]=matrix polynomials
##       result[2]=set of parameters/unknowns in result[1]
D1ParaMaskVarC:=proc(rowdim,coldim,low,high,para)
    local varr,vari,tmp,result;
    varr:=evaln(para||r):
    vari:=evaln(para||i):
    result:=Vector(2):
    tmp:=D1ParaMaskVar(rowdim,coldim,low,high,varr):
    result[1]:=tmp[1]:
    result[2]:=tmp[2]:
    tmp:=D1ParaMaskVar(rowdim,coldim,low,high,vari):
    result[1]:=result[1]+I*tmp[1]:
    result[2]:=result[2] union tmp[2]:
    result:=result:
end proc:


##Parametrize a Matrix of Laurent polynomials supported on [low,high] with variable z and complex coefficients
D1ParaMaskC:=proc(rowdim,coldim,low,high,para)
    local tmp,result;
    tmp:=D1ParaMaskVarC(rowdim,coldim,low,high,para):
    result:=tmp[1]:
    result:=result:
end proc:

##Parametrize a Matrix of polynomials in variable xi with  O(|\xi|^degord) and with complex coefficients
##ouput: result[1]=matrix polynomials
##       result[2]=set of parameters/unknowns in result[1]
D1ParaPolyVarC:=proc(rowdim,coldim,degord,para)
    local varr,vari,tmp,result;
    result:=Vector(2):
    varr:=evaln(para||r):
    vari:=evaln(para||i):
    tmp:=D1ParaPolyVar(rowdim,coldim,degord,varr):
    result[1]:=tmp[1]:
    result[2]:=tmp[2]:
    tmp:=D1ParaPolyVar(rowdim,coldim,degord,vari):
    result[1]:=result[1]+I*tmp[1]:
    result[2]:=result[2] union tmp[2]:
    result:=result:
end proc:


##Parametrize a Matrix of polynomials in variable xi with  O(|\xi|^degord) and with complex coefficients
D1ParaPolyC:=proc(rowdim,coldim,degord,para)
    local tmp,result;
    tmp:=D1ParaPolyVarC(rowdim,coldim,degord,para):
    result:=tmp[1]:
    result:=result:
end proc:


##Output:set of equations for prescribed symmetry of matrix filters
##Matrix-valued mpoly can be any size of matrix of Laurent polys.
##Symmetry:\hat{u}(\xi)=symb(dilation*xi)\overline{\hat{u}(\xi)}syma(\xi)^{-1}
#for simplicity, symphi/sympsi are Vectors (instead diagonal matrices)
##For symmetry of filters, see Exercise 6.42(items(c), (g)) on pages 575-576 of the book.
##choice<0 for complex symmetry; and choice>=0 for symmetry
D1SymMaskEqsC:=proc(mpoly,dilation,syma,symb,choice)
    local mpolyri,result;
    mpolyri:=D1RealImagPoly(mpoly):
    result:=D1SymMaskEqs(mpolyri[1],dilation,syma,symb):
    if (choice<0) then
        result:=result union D1SymMaskEqs(mpolyri[2],dilation,-syma,symb):
    else
        result:=result union D1SymMaskEqs(mpolyri[2],dilation,syma,symb):
    end if:
    result:=result:
end proc:


##Parametrize a Matrix of Laurent polynomials in z with (complex) symmetry
##choice<0 for complex symmetry; and choice>=0 for symmetry
D1ParaMaskSymC:=proc(rowdim,coldim,low,high,syma,symb,para,choice)
    local tmp,solt,tpara,result:
    tpara:=evaln(t||para):
    tmp:=D1ParaMaskVarC(rowdim,coldim,low,high,tpara):
    result:=SolveEqsVar(D1SymMaskEqsC(tmp[1],dilation,syma,symb,choice),tmp[2],para,tmp[1],0): #show number of free parameters
    result:=result:
end proc:

##Output: a dual filter for a given filter with prescribed sum rules
##and with prescribed symmetry for complex filters
##choice<0 for complex symmetry; and choice>=0 for symmetry
D1DualMaskSymC:=proc(apoly,dilation,low,high,dsrord,syma,para,choice)
    local j,multi,dualvgu,tmppoly,result:
    multi:=RowDimension(apoly):
    tmppoly:=D1ParaMaskSymC(multi,multi,low,high,syma,syma,para,choice):
    dualvgu:=D1Dualvgu(apoly,dilation,dsrord):
    tmppoly:=D1SumRulePoly(tmppoly,dilation,dualvgu,dsrord):
    result:=SolveEqs(D1BiorthEqs(apoly,tmppoly,dilation,0),tmppoly):
    result:=result:
end proc:



###################################################
#########End of routines ending with C#############
###################################################



