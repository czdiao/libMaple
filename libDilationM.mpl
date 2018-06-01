###############################
#   Library from Prof Han proc
#
#   Used for 1D arbitrary dilation
#
#   May, 2018
#
###############################

read "fw1d.mpl":

# Plot Refinable function phi with any dilation
# Wrapper calling D1Plotphi
D1MPlotphi:= proc(a, dilation, level)
    # To use:
    #   a is the Laurent polynomials for lowpass
    #   level is the number of levels we perform subdivision
    local phidata, xmin, xmax, ymin, ymax, yeps, apoly:

    apoly := Matrix([[a]]):
    phidata := D1Plotphi(apoly, apoly, dilation, level):

    xmin := min(phidata[0]):
    xmax := max(phidata[0]):

    ymin := min(phidata[1]):
    ymax := max(phidata[1]):
    yeps := (ymax - ymin)/20:
    ymin := ymin - yeps:
    ymax := ymax + yeps:

    plot(phidata[0], phidata[1], xmin..xmax, ymin..ymax, axes=boxed, tickmarks=[default, 4], color="Black");
end proc:


# Plot wavelet function psi with any dilation
# Wrapper calling D1Plotphi
D1MPlotpsi:= proc(a, b, dilation, level)
    # To use:
    #   a, b are the Laurent polynomials for lowpass and highpass
    #   level is the number of levels we perform subdivision
    local psidata, xmin, xmax, ymin, ymax, yeps, apoly, bpoly:

    apoly := Matrix([[a]]):
    bpoly := Matrix([[b]]):
    psidata := D1Plotphi(apoly, bpoly, dilation, level):

    xmin := min(psidata[0]):
    xmax := max(psidata[0]):

    ymin := min(psidata[1]):
    ymax := max(psidata[1]):
    yeps := (ymax - ymin)/20:
    ymin := ymin - yeps:
    ymax := ymax + yeps:

    plot(psidata[0], psidata[1], xmin..xmax, ymin..ymax, axes=boxed, tickmarks=[default, 4], color="Black");
end proc:

# wrapper for D1L2SmoothnessPoly(apoly,dilation,maxsr,errtol):
# Computing the L2 smoothness exponent
D1Msm2 := proc(a, dilation, maxsr)
    # maxsr is the maximum sum rule we guessed
    local apoly, errtol, sm:

    errtol := 0.0001:
    apoly := Matrix([a]):
    sm := D1L2SmoothnessPoly(apoly,dilation,maxsr,errtol):

    return evalf(sm[0], 6):

end proc:

getGammaCoset_M := proc(a, dilation, gammaj)
    local ashift, ldeg, deg, j, c, agamma:

    ashift := collect(a/z^gammaj, z): # shift the gammaj coset to 0-coset
    ldeg := ldegree(ashift, z):
    deg := degree(ashift, z):

    ldeg := floor(ldeg/dilation) * dilation:
    deg := ceil(deg/dilation) * dilation:

    agamma := 0:
    for j from ldeg to deg by dilation do
        c := coeff(ashift, z, j):
        agamma := agamma + c * z^(j/dilation):
    end do;

    return agamma:
end proc:

getCoset_M := proc(a, dilation)
    description "Given a representative number gammaj, compute this coset sequence of dilation M":
    # return a (col) vector of all coset sequences

    local ind, c, dd, ldeg, deg, a0, j:

    a0 := Vector(dilation):
    ldeg := ldegree(a, z):
    deg := degree(a, z):

    for j from ldeg by 1 to deg do
        ind := modp(j, dilation):   # coset index
        c := coeff(a, z, j):    # coefficent
        dd := (j-ind)/dilation:  # downsampled degree in each coset
        a0[ind+1] := a0[ind+1] + c*z^dd:
    end do;

    return a0:
end proc:

getFilter_M := proc(aCoset, dilation)
    # get original filter from coset sequences.
    local a, j:
    a := 0:

    for j from 0 to (dilation-1) do
        a := a + z^j * eval(aCoset[j+1], z = z^dilation)
    end do;

    return a:
end proc:

getRankOne_M := proc(a, dilation)
    # get the rank one matrix for each filter
    local omega, aVec, j, aOne:
    omega := exp(I*2*Pi/dilation):
    aVec := Matrix(dilation, 1):

    for j from 1 to dilation do
        aVec[j, 1] := eval(a, z = omega^(j-1) * z):
    end do;
    aOne := aVec.hermConj(aVec):
    return aOne:
end proc:

getM0_M := proc(a, Theta, dilation)
    local omega, aVec, M0, j:
    omega := exp(I*2*Pi/dilation):
    aVec := Matrix(dilation, 1):

    for j from 1 to dilation do
        aVec[j, 1] := eval(a, z = omega^(j-1) * z):
    end do;
    M0 := - eval(Theta, z = z^dilation) * aVec.hermConj(aVec):

    for j from 1 to dilation do
        M0[j, j] := M0[j, j] + eval(Theta, z =  omega^(j-1) * z):
    end do;

    return M0:
end proc:

getMnbFactored_M := proc(M0, dilation, nb)
    # M0 is calculated from getM0_M(a, Theta, dilation)
    local omega, Mnb, P, j, k:
    omega := exp(I*2*Pi/dilation):
    Mnb := M0:
    for j from 1 to dilation do
        for k from 1 to dilation do
            P := (1 - omega^(k-1) * z)^nb * (1 - omega^(1-j) / z)^nb:
            if not Ldivide(Mnb[j,k], P) then
                print(j):
                print(k):
                print(P):
                error("Cannot divide the vanishing moment!"):
            end if;
            Mnb[j,k] := Lquo(Mnb[j,k], P):
        end do;
    end do;

    return Mnb:
end proc:

getNnbFactored_M := proc(Mnb, dilation)
    local omega, Finv, j, k, Nnb:
    omega :=  exp(I*2*Pi/dilation):
    Finv := Matrix(dilation, dilation, 0):
    for j from 1 to dilation do
        for k from 1 to dilation do
            Finv[j,k] := 1/dilation * omega^(-(j-1)*(k-1)) * z^(1-j):
        end do;
    end do;

    Nnb := Finv . Mnb . hermConj(Finv):
    return expand~(rationalize~(simplify~(Nnb))):
end proc:

CheckPR_M := proc(FB, Sig, dilation)
    # FB is a Vector with all filters, lowpass and highpass, Sig is a list/Vector of size (s+1)
    # Check PR on polyphase, to avoid the complex calculations
    # Should return dilation*dilation 0 matrix
    local j, n, output, fVec:
    n := nops(Sig):
    output := 1/dilation * IdentityMatrix(dilation):
    for j from 1 to n do
        fVec := Matrix(getCoset_M(FB[j], dilation)):
        output := output - Sig[j] * fVec.hermConj(fVec):
    end do;

    return simplify~(output):
end proc:


(*
getEmu_M := proc( vm, dilation)
    # Used to factor out vanishing moments to for polyphase matrices.
    # See Theorem 8 in
    #   "QUASI-TIGHT FRAMELETS WITH DIRECTIONALITY OR HIGH VANISHING MOMENTS DERIVED FROM ARBITRARY REFINABLE FUNCTIONS"
    # vm: vanishing moment to be factored out
    local deltaMu, deltaCoset, j, k, Emu:
    deltaMu := (1-z)^vm:

    Emu := Matrix(dilation, dilation):

    for j from 1 to dilation do
        for k from 1 to dilation do
            Emu[j, k] := getGammaCoset_M(deltaMu, dilation, k-j):
        end do;
    end do;

    return hermConj(Emu):   # Notation is different with the paper.
end proc:

getFactoredN:= proc(a, dilation, vm)
    local Na, aCoset, EmuInv, EmuInvStar:

    aCoset := getCoset_M(a, dilation):
    Na := 1/dilation * IdentityMatrix(dilation) - Matrix(aCoset).hermConj(Matrix(aCoset)):
    EmuInv := MatrixInverse(getEmu_M(vm, dilation)):
    EmuInvStar := MatrixInverse(hermConj(getEmu_M(vm, dilation))):
    Na := EmuInv.Na.EmuInvStar:

    return simplify~(Na):
end proc: *)




