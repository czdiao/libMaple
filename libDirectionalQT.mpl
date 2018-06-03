#########################################################
#
#   Tool functions to find quasi-tight framelet filter
#   banks with directionality/vanishing moments.
#
#   Used in the high-dimensional quasi-tight paper.
#
#
#   Chenzhe
#   Nov, 2017
#


HermSqMatrix := proc(A, zlist)
    description "A is a Hermitian matrix of Laurent polynomials, split into SOS/DOS: A=sum_j epsilon_j u_j.u_j^*, where u_j only supported on two points";

    local n, s, j, k, c, cj, u, epsilon, Atmp, Ajk, outPairs, deg:

    n := RowDimension(A):

    Atmp := simplify(A):
    outPairs := []:
    # off diagonal
    for j from 1 to (n-1) do
        for k from (j+1) to n do
            Ajk := collect(A[j, k], zlist, distributed):
            Ajk := getAllTerms(Ajk, zlist):
            
            #print(nops(Ajk)):
            if not evalb(A[j,k]=0) then
            for s from 1 to nops(Ajk) do
                u := Vector(n):
                c := coeffs(Ajk[s]):
                Atmp[j,j] := Atmp[j,j] + c:
                Atmp[k,k] := Atmp[k,k] + c:

                epsilon := -signum(c):
                cj := sqrt(abs(c)):
                u[j] := epsilon*Ajk[s]/cj:
                u[k] := cj:

                outPairs := [op(outPairs), [u, epsilon]]:
                #print(s):
            end do;      
            end if;

        end do;
    end do;

    # non-const terms on the diagonals
    for j from 1 to n do
        Ajk := getAllTerms(Atmp[j, j], zlist):
        Ajk := remove(x->isConstTerm(x, zlist), Ajk):
        #print(nops(Ajk)):
        #print(Ajk):
        while not evalb(nops(Ajk) = 0) do
                u := Vector(n):
                c := coeffs(Ajk[1]):

                epsilon := -signum(c):
                cj := sqrt(abs(c)):
                u[j] := epsilon*Ajk[1]/cj + cj:

                outPairs := [op(outPairs), [u, epsilon]]:

                Atmp[j,j] := simplify(Atmp[j, j] + 2 * c - Ajk[1] - hcND(Ajk[1], zlist)) :
                Ajk := getAllTerms(Atmp[j, j], zlist):
                Ajk := remove(x->isConstTerm(x, zlist), Ajk):
        end do;
    end do;


    # constant terms on the diagonals
    for j from 1 to n do
        c := Atmp[j, j]:
        if evalb(c=0) then
            break:
        end if;
        epsilon := signum(c):
        u := Vector(n):
        u[j] := sqrt(abs(c)):
        outPairs := [op(outPairs), [u, epsilon]]:
    end do;

    # scalar case, change all u to be a Lpoly, not a vector
    if evalb(n=1) then
        for j from 1 to nops(outPairs) do
            outPairs[j][1] := outPairs[j][1][1]:
        end do;
    end if;

    return outPairs:
    
end proc:


DividedDiffSp := proc(u, m)
    description "Given 2D filter u with vanishing moments >= m, separate u into summation of divided difference.":
    # output: Vmu[mu, m-mu]

    local ld1, ld2, upoly, tmp, deg, td, v1, v2, j1, j2, Vmu:

    ld1 := ldegree(u, z[1]):
    ld2 := ldegree(u, z[2]):
    upoly := simplify(u / z[1]^ld1 / z[2]^ld2):
    deg := degree(upoly):

    Vmu := table():
    for v1 from 0 to m do
        v2 := m - v1:
        Vmu[v1, v2]:= 0:
    end do;

    for td from m to deg do
        for j1 from td by (-1) to 0 do
            j2 := td - j1:
            tmp := coeftayl(upoly, [z[1], z[2]] = [1, 1], [j1, j2] ):
            v1 := min(j1, m):
            v2 := m - v1:
            tmp := tmp * (z[1]-1)^(j1 - v1) * (z[2]-1)^(j2 - v2):
            Vmu[v1, v2]:= Vmu[v1, v2] + tmp:
        end do;
    end do;

    # shift back to Laurent polynomial
    for v1 from 0 to m do
        v2 := m - v1:
        Vmu[v1, v2]:= simplify(Vmu[v1, v2] * z[1]^ld1 * z[2]^ld2):
    end do;

    # verify the result
    tmp := 0:
    for v1 from 0 to m do
        v2 := m - v1:
        tmp:= tmp + (z[1]-1)^v1 * (z[2]-1)^v2 * Vmu[v1, v2]:
    end do;
    tmp := simplify(tmp - u):
    if not evalb(tmp=0) then
        error("Divided Difference separation error!"):
    end if;

    #Vlist := [seq(simplify(Vmu[v1, m-v1]), v1 = 0..m)]:
    # return a table: Vmu[v1, v2] is the coeff Lpoly of (z[1]-1)^v2 * (z[2]-1)^v2
    return eval(Vmu):
end proc:


HermSqVM := proc(u, m)
    description "u is a Hermitian scalar Lpoly with 2m order VM and real coeff, split it into SOS/DOS of filters with m order VM.":

    local Unu, theta, v1, v2, alpha1, alpha2, beta1, beta2, epsilon_l, u_j, outPairs, tmp, tmp2, eta, mu1, mu2, etaSplit, ul:

    Unu := DividedDiffSp(u/2, 2*m):

    outPairs := []:
    for v1 from 1 by 2 to (2*m-1) do
        v2 := 2*m - v1:
        if not evalb(Unu[v1, v2]=0) then
            alpha1 := min(v1, m):
            alpha2 := m - alpha1:
            beta1 := v1 - alpha1:
            beta2 := v2 - alpha2:

            u_j := (1/z[1]-1)^alpha1 * (1/z[2]-1)^alpha2:
            u_j := u_j + (z[1]-1)^beta1 * (z[2]-1)^beta2 * Unu[v1, v2]:
            outPairs := [op(outPairs), [ u_j, 1]]:

            tmp := (-1)^(alpha1 + alpha2) /z[1]^alpha1/z[2]^alpha2:
            Unu[2*alpha1, 2*alpha2] := Unu[2*alpha1, 2*alpha2] - tmp/2:

            tmp := Unu[v1, v2]:
            tmp := tmp * hc2D(tmp):
            tmp := tmp * (-1)^(beta1 + beta2) /z[1]^beta1/z[2]^beta2 :
            Unu[2*beta1, 2*beta2] := Unu[2*beta1, 2*beta2] - tmp/2:
        end if;
    end do;

    for v1 from 0 by 2 to 2*m do
        v2 := 2*m - v1:
        mu1 := v1/2:
        mu2 := v2/2:
        if not evalb(simplify(Unu[v1, v2]=0)) then
            tmp := Unu[v1, v2] * z[1]^mu1 * z[2]^mu2:
            tmp := tmp + hc2D(tmp):
            eta := (-1)^(mu1 + mu2) * tmp:
            eta := simplify(eta):

            etaSplit := HermSqMatrix(Matrix([eta]), [z[1], z[2]]):

            ul := [seq(etaSplit[j][1]*(z[1]-1)^mu1*(z[2]-1)^mu2, j=1..nops(etaSplit))]:
            epsilon_l := [seq(etaSplit[j][2], j = 1..nops(etaSplit))]:

            tmp:= [seq([ul[j], epsilon_l[j]], j = 1..nops(etaSplit))]:
            outPairs:= [op(outPairs), op(tmp)]:
        end if;

    end do;

    # verify solution
    tmp := 0:
    for v1 from 1 to nops(outPairs) do
        tmp2 := outPairs[v1][1]:
        tmp := tmp + outPairs[v1][2] * tmp2 * hc2D(tmp2):
    end do;
    tmp := simplify(tmp - u):
    if not evalb(tmp = 0) then
        error("Error in HermSqVM: solution not equal!"):
    end if;

    return outPairs:

end proc:

# See Step (S2) of Thm8. Break Auv into A1^*.A2.
FactorAuv := proc(Auv)
    local j, nrow, ncol, rlist, clist, A1, A2, tmp:

    nrow := RowDimension(Auv):
    ncol := ColumnDimension(Auv):
    rlist := []:
    clist := []:
    # find the nonzero row indices
    for j from 1 to nrow do
        tmp := Auv[j, ..]:
        if not is0Matrix(tmp) then
            rlist := [op(rlist), j]:
        end if;
    end do;
    # find the nonzero col indices
    for j from 1 to ncol do
        tmp := Auv[.., j]:
        if not is0Matrix(tmp) then
            clist := [op(clist), j]:
        end if;
    end do;

    if evalb(nops(rlist)<=nops(clist)) then
        A1 := Matrix(nrow, nops(rlist)):
        for j from 1 to nops(rlist) do
            A1[rlist[j], j] := 1:
        end do;

        A1 := Transpose(A1):
        A2 := A1.Auv:
    else
        A2 := Matrix(nops(clist), ncol):
        for j from 1 to nops(clist) do
            A2[j, clist[j]] := 1:
        end do;

        A1 := Auv.Transpose(A2):
        A1 := hcMatrix2D(A1):
    end if;


    # Verify
    tmp := simplify(hcMatrix2D(A1).A2 - Auv):
    if not is0Matrix(tmp) then
        error("Error in Split Auv into hcMatrix2D(A1).A2"):
    end if;

    return A1, A2:
end proc:


