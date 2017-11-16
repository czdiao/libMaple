#########################################################
#
#   Tool functions to find quasi-tight framelet filter
#   banks with directionality/vanishing moments.
#
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
            
            print(nops(Ajk)):
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
                print(s):
            end do;
        end do;
    end do;

    # non-const terms on the diagonals
    for j from 1 to n do
        Ajk := getAllTerms(Atmp[j, j], zlist):
        Ajk := remove(x->isConstTerm(x, zlist), Ajk):
                print(nops(Ajk)):
                print(Ajk):
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
        epsilon := signum(c):
        u := Vector(n):
        u[j] := sqrt(abs(c)):
        outPairs := [op(outPairs), [u, epsilon]]:
    end do;


    return outPairs:
    
end proc: