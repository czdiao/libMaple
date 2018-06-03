#########################################################
#
#   Tool functions to plot z domain Laurent polynomials
#
#   Define z=exp(-I\xi)
#
#   This library relies on libLpoly.mpl
#
#   Chenzhe
#   Jun, 2018
#

with(LinearAlgebra):
with(plots):
with(PolynomialTools):



PlotPolyZreal := proc(poly)
description "Plot real part value of z-domain Polynomial poly in -Pi..Pi";
## poly should be in z-domain

	local t, polyt, p1, p2, yreal, ymin, ymax, yeps, xiVal:
	polyt := eval(simplify(poly), z = exp(-I*xi)):

    xiVal := [seq((j*2*Pi/100)-Pi, j = 0..50 )]:
    yreal := [seq( Re(eval(polyt, xi = xiVal[j])), j = 1..nops(xiVal) )]:
    ymin := min(yreal):
    ymax := max(yreal):
    yeps := (ymax - ymin)/20:

	#p1:= plot([Re(polyt)], xi = -Pi..Pi, axes = boxed, labels = [``, ``]):
    #p2:= plot(0, x=-Pi..Pi, linestyle=dash, color="Blue", labels = [``, ``]):
    p1:= plot([Re(polyt)], xi=-(Pi*1.01)..Pi, tickmarks = [5, 4], axes = boxed, color="Black", labels = [``, ``]):
    p2:= plot(0, x=-(Pi*1.01)..Pi, y = -yeps..yeps, tickmarks = [5, default], linestyle=dash, color="Black", labels = [``, ``], axes = boxed):

    display({p1, p2});
end proc:



PlotPolyZ := proc(poly)
description "Plot abs value of z-domain Polynomial poly in -Pi..Pi";
## poly should be in z-domain

	local xi, polyt:
	polyt := eval(simplify(poly), z = exp(-I*xi)):
	plot([abs(polyt)], xi = -Pi..Pi, color="Black");
end proc:


PlotPolyZn := proc(polylist)
    description "Plot multiple z-domain Laurent Polynomials simultaneously":
    # input should be a list: [p1, p2, ..., pn]
    # linestype in order of: solid, dot, dash, dashdot, longdash, spacedash, spacedot. 

	local polyt, n, iter:

	n := nops(polylist):
    polyt := eval(simplify(polylist), z=exp(-I*xi)):
    #plot(evalf(abs~(polyt)), xi=-Pi..Pi, tickmarks = [spacing((1/2)*Pi), default], legend = legendlist, linestyle=[seq(i, i=1..n)], axes=boxed );
    plot(evalf(abs~(polyt)), xi=-(Pi*1.01)..Pi, tickmarks = [5, default], linestyle=[seq(i, i=1..n)], axes=boxed, labels = [``, ``] );


end proc:

PlotOrthCheck := proc(a)
    description "Check orthogonality of filter, blue dashed line is y=1";
    # Plot (1) y=1, (2) |a(\xi)|^2 + |a(\xi + Pi)|^2

    local astar, ainv, ainvstar, p, p0, p1, p2:

    astar:= hc(a):
    ainv:= eval(a, z=-z):
    ainvstar:= hc(ainv):

    p:=a*astar + ainv*ainvstar:

    #p0:=pointplot([[0, 0]]):
    p1:=plot(1, x=-Pi..Pi, linestyle=dash, color="Blue"):
    p2:=PlotPolyZ(p):

    display({p1, p2}):
    #display(p2):
end proc:


PlotRootsZ := proc(poly)
    description "Plot all roots of Laurent polynomial poly in complex plane, also label their multiplicity";
    # Input Laurent polynomial should have fractional coeffs, rather than floating point coeffs.
    # Otherwise Splits() will fail. We could use factors(poly, complex) instead
    # We also add the plot of unit cirlce

	local f, r, nr, t, proots, pcircle, ptext, ii, mult, textpt, ldeg:
	
	f:=Splits(collect(poly,z), z):
	#ldeg:=ldegree(poly, z): f:=factors(collect(poly*z^(-ldeg),z), complex):
	
	f:= f[2]:   # all factors
	nr:= nops(f):   # number of distinct factors
	r:=[seq(-eval(f[ii][1],z=0), ii=1..nr)]:    # list all roots
	mult:=[seq(f[ii][2], ii=1..nr)]:    # multiplicity of the roots
	textpt:=[seq( [Re(r[ii]), Im(r[ii]), typeset(mult[ii])] , ii=1..nr)]:
	
	proots:=complexplot(r, style = point, symbol =  soliddiamond, symbolsize = 20):
	pcircle:=complexplot(cos(t)+I*sin(t), t=0 .. 2*Pi):
	ptext:= textplot(textpt, 'align'={'below', 'left'}, font = [TIMES, BOLD, 20]):
	
	display({proots, pcircle, ptext}):
	
end proc:

PlotRootsZ_num := proc(poly)
    description "Plot all roots of Laurent polynomial poly in complex plane, also label their multiplicity";
    # Input Laurent polynomial could have floating point coeffs.
    # use factors(poly, complex) instead of Splits
    # We also add the plot of unit cirlce

	local f, r, nr, t, proots, pcircle, ptext, ii, mult, textpt, ldeg:
	
	#f:=Splits(collect(poly,z), z):
	ldeg:=ldegree(poly, z): f:=factors(collect(poly*z^(-ldeg),z), complex):
	
	f:= f[2]:   # all factors
	nr:= nops(f):   # number of distinct factors
	r:=[seq(-eval(f[ii][1],z=0), ii=1..nr)]:    # list all roots
	mult:=[seq(f[ii][2], ii=1..nr)]:    # multiplicity of the roots
	textpt:=[seq( [Re(r[ii]), Im(r[ii]), typeset(mult[ii])] , ii=1..nr)]:
	
	proots:=complexplot(r, style = point, symbol =  soliddiamond, symbolsize = 20):
	pcircle:=complexplot(cos(t)+I*sin(t), t=0 .. 2*Pi):
	ptext:= textplot(textpt, 'align'={'below', 'left'}, font = [TIMES, BOLD, 20]):
	
	display({proots, pcircle, ptext}):
	
end proc:


PlotPhi := proc(a, n)
    description "Plot scalar refinable function":

    local xval, yval, iter, xmin, xmax, k, aRound, ymin, ymax, yeps:
    yval:= 1:   # choose initial function to be delta in cascade algorithm
    aRound := RoundingLpoly(a):

    for iter from 1 to n do
        yval:= 2 * aRound * eval(yval, z=z^2):
        yval:= collect(yval, z):
    od:

    xmin:= ldegree(yval, z):
    xmax:= degree(yval, z):
    xval:= [seq( k/2^n, k=xmin..xmax )]:

    yval:=CoefficientVector(collect(yval/z^xmin, z), z):

    xmin := floor(min(xval)):
    xmax := ceil(max(xval)):
    ymax := max(yval):
    ymin := min(yval):
    yeps := (ymax - ymin)/20:

    plot(xval, yval, xmin..xmax, (ymin-yeps)..(ymax + yeps), axes=boxed, tickmarks=[default, 4], color="Black");

end proc:


PlotPsi := proc(a, b, n)
    description "Plot scalar wavelet function":

    local xval, yval, iter, xmin, xmax, k, aRound, bRound, ymin, ymax, yeps:
    aRound := RoundingLpoly(a):
    bRound := RoundingLpoly(b):
    yval:= bRound:   

    for iter from 1 to (n-1) do
        yval:= 2 * aRound * eval(yval, z=z^2):
        yval:= collect(yval, z):
    od:

    xmin:= ldegree(yval, z):
    xmax:= degree(yval, z):
    xval:= [seq( k/2^n, k=xmin..xmax )]:

    yval:=CoefficientVector(collect(yval/z^xmin, z), z):

    xmin := floor(2*min(xval))/2:
    xmax := ceil(2*max(xval))/2:
    ymax := max(yval):
    ymin := min(yval):
    yeps := (ymax - ymin)/20:

    plot(xval, yval, xmin..xmax, (ymin-yeps)..(ymax + yeps), axes=boxed, tickmarks=[default, 4], color="Black");
end proc:

PlotPhiShift := proc(a, u, n)
    description "Plot linear combination of scalar refinable function":
    # phi_0 is defined by lowpass a, plot phi(x) = sum_k u(k) phi_0(x-k)
    # This is used when the integer shifts of phi is not stable.
    # Special Case:
    #       PlotPhi(a, n) = PlotPhiShift(a, 1, n)

    local xval, yval, iter, xmin, xmax, k, aRound, ymax, ymin, yeps:
    yval:= RoundingLpoly(u):   # choose initial function to be delta in cascade algorithm
    aRound := RoundingLpoly(a):

    for iter from 1 to n do
        yval:= 2 * aRound * eval(yval, z=z^2):
        yval:= collect(yval, z):
    od:

    xmin:= ldegree(yval, z):
    xmax:= degree(yval, z):
    xval:= [seq( k/2^n, k=xmin..xmax )]:

    yval:=CoefficientVector(collect(yval/z^xmin, z), z):

    xmin := floor(min(xval)):
    xmax := ceil(max(xval)):
    ymax := max(yval):
    ymin := min(yval):
    yeps := (ymax - ymin)/20:

    plot(xval, yval, xmin..xmax,  (ymin-yeps)..(ymax + yeps),  axes=boxed, tickmarks=[default, 4], color="Black");
end proc:


PlotPsiShift := proc(a, u, b, n)
    description "Plot scalar wavelet function":
    # phi is defined by u*phi_0, see PlotPhiShift() above
    # This is used when the integer shifts of phi is not stable.
    # Special Case:
    #       PlotPsi(a, b, n) = PlotPsiShift(a, 1, b, n)

    local xval, yval, iter, xmin, xmax, k, aRound, bRound, ymax, ymin, yeps:
    aRound := RoundingLpoly(a):
    bRound := RoundingLpoly(b):
    yval:= RoundingLpoly(bRound * u):

    for iter from 1 to (n-1) do
        yval:= 2 * aRound * eval(yval, z=z^2):
        yval:= collect(yval, z):
    od:

    xmin:= ldegree(yval, z):
    xmax:= degree(yval, z):
    xval:= [seq( k/2^n, k=xmin..xmax )]:

    yval:=CoefficientVector(collect(yval/z^xmin, z), z):

    xmin := floor(2*min(xval))/2:
    xmax := ceil(2*max(xval))/2:
    ymax := max(yval):
    ymin := min(yval):
    yeps := (ymax - ymin)/20:

    plot(xval, yval, xmin..xmax, (ymin-yeps)..(ymax + yeps), axes=boxed, tickmarks=[default, 4], color="Black");
end proc:


StemPlot := proc(a)
    description "stem plot of a filter in Laurent Polynomial":

    local ldeg, deg, j, xx, yy, ymin, ymax, meps:

    ldeg := ldegree(a, z):
    deg := degree(a, z):

    xx := [seq(j, j=ldeg..deg)]:
    yy := [seq(coeff(a, z, j), j=ldeg..deg)]:

    ymin := min(min(yy), 0):
    ymax := max(max(yy), 0):
    meps := (ymax - ymin)/20:

    #display( DynamicSystems:-DiscretePlot(xx, yy,style=stem), 
    #    plot(0, x=(ldeg-1)..(deg+1), y=min(yy)-meps..max(yy)+meps, axes=boxed) );
    display([plot(xx, yy, style=point, symbol=circle, axes=boxed),
                plot(0, x=(ldeg-1)..(deg+1), y=ymin-meps..ymax+meps, axes=boxed, 
                    tickmarks=[[seq(i, i=ldeg..deg)], default]),
                plot([seq([[xx[i],yy[i]],[xx[i],0]],i=1..nops(xx))],colour=black)])

end proc:


PlotHermEigenVal := proc(A)
    description "point plot of the eigenvalues of a Hermite matrix of Laurent polynomials":
    # A is assumed to be Hermite, so the eigenvalues are real. Only plot the real part.
    local dim, dim2, eiv, k, points, N, xi, p1, p2, ymin, ymax, meps:

    N := 200:
    xi := [seq((j*2*Pi/N)-Pi, j = 0..N )]:
    dim, dim2 := Dimension(A):

    eiv := [seq(Eigenvalues(evalf(eval(A, z = cos(xi[j])+I*sin(xi[j])))), j = 1..nops(xi))]:

    points := {}:
    for k from 1 to dim do
        points := points union {seq([xi[j], Re(eiv[j][k])], j = 1..nops(xi))}:
    end do;

    ymax := max(max([seq(points[j][2], j = 1..nops(points))]), 0):
    ymin := min(min([seq(points[j][2], j = 1..nops(points))]), 0):
    meps := (ymax - ymin)/20:

    p1:= pointplot(points, symbol = point, axes=boxed):
    p2:= plot(0, x=-Pi..Pi, y=ymin-meps..ymax+meps, linestyle=dash, tickmarks = [5, default], color="Blue", labels = [``, ``]):

    display({p1, p2}):

end proc:

############# 2D wavelet/mask plot ##############
# 2-dimensional frequency plot
PlotPolyZD2 := proc(a)
description "Plot a 2d mask in frequency domain, a is a polynomial of z[1] and z[2]";
	local xi1, xi2, polyt:
	polyt := eval(simplify(a), {z[1] = exp(-I*xi1), z[2] = exp(-I*xi2)}):
	plot3d([abs(polyt)], xi1 = -Pi..Pi, xi2= -Pi..Pi);
end proc:

PlotPolyZD2real := proc(a)
description "Plot a 2d mask in frequency domain, a is a polynomial of z[1] and z[2]";
	local xi1, xi2, polyt:
	polyt := eval(simplify(a), {z[1] = exp(-I*xi1), z[2] = exp(-I*xi2)}):
	plot3d([Re(polyt)], xi1 = -Pi..Pi, xi2= -Pi..Pi);
end proc:

# General dilation matrix
PlotPhi2D := proc(a, MM, n)
    description "Plot 2D refinable function with general dilation matrix MM":

    local xval, yval, iter, k1, k2, aRound, Minv, dM, X1, X2, Y, MinvN, ld1, d1, ld2, d2, mm, nn, di1, di2, c, di1m, di2m, x1, x2, xvec:
    Minv := MatrixInverse(MM):
    dM := Determinant(MM):
    MinvN := Minv^n:
    yval:= 1:   # choose initial function to be delta in cascade algorithm
    aRound := a: #RoundingLpoly(a):

    for iter from 1 to n do
        yval:= dM * aRound * upsampleM(yval, MM):
        yval:= collect(yval, [z[1], z[2]], distributed):
    od:

    # find the point coords, put into 3 matrices (X1, X2, Y)
    # notice that the matrix Y start from bottom left corner
    # horizontal is X1 direction (z[1])
    # vertical is X2 direction (z[2])
    ld1 := ldegree(yval, z[1]):
    d1 := degree(yval, z[1]):
    ld2 := ldegree(yval, z[2]):
    d2 := degree(yval, z[2]):

    mm := d1 - ld1 +1:
    nn := d2 - ld2 +1:
    X1 := Matrix(nn, mm):
    X2 := Matrix(nn, mm):
    Y  := Matrix(nn, mm):

    for di1 from ld1 to d1 do
        for di2 from ld2 to d2 do
            c := coeff(coeff(yval, z[1], di1), z[2], di2):
	        di1m := di1 - ld1 +1:
	        di2m := di2 - ld2 +1:
	        Y[nn-di2m+1, di1m] := c:
            xvec := MinvN.<di1, di2>:
            X1[nn-di2m+1, di1m] := xvec[1]:
            X2[nn-di2m+1, di1m] := xvec[2]:
        end do:
    end do:

    return X1, X2, Y:
    #plot(xval, yval, xmin..xmax, axes=boxed);

end proc:

PlotPsi2D := proc(a, b, MM, n)
    description "Plot 2D refinable function with general dilation matrix MM":

    local xval, yval, iter, k1, k2, aRound, Minv, dM, X1, X2, Y, MinvN, ld1, d1, ld2, d2, mm, nn, di1, di2, c, di1m, di2m, x1, x2, xvec:
    Minv := MatrixInverse(MM):
    dM := Determinant(MM):

    MinvN := Minv^n:
    yval:= b:   # choose initial function to be delta in cascade algorithm
    aRound := a: #RoundingLpoly(a):

    for iter from 1 to n do
        yval:= dM * aRound * upsampleM(yval, MM):
        yval:= collect(yval, [z[1], z[2]], distributed):
    od:

    # find the point coords, put into 3 matrices (X1, X2, Y)
    # notice that the matrix Y start from bottom left corner
    # horizontal is X1 direction (z[1])
    # vertical is X2 direction (z[2])
    ld1 := ldegree(yval, z[1]):
    d1 := degree(yval, z[1]):
    ld2 := ldegree(yval, z[2]):
    d2 := degree(yval, z[2]):

    mm := d1 - ld1 +1:
    nn := d2 - ld2 +1:
    X1 := Matrix(nn, mm):
    X2 := Matrix(nn, mm):
    Y  := Matrix(nn, mm):

    for di1 from ld1 to d1 do
        for di2 from ld2 to d2 do
            c := coeff(coeff(yval, z[1], di1), z[2], di2):
	        di1m := di1 - ld1 +1:
	        di2m := di2 - ld2 +1:
	        Y[nn-di2m+1, di1m] := c:
            xvec := MinvN.<di1, di2>:
            X1[nn-di2m+1, di1m] := xvec[1]:
            X2[nn-di2m+1, di1m] := xvec[2]:
        end do:
    end do:

    return X1, X2, Y:
    #plot(xval, yval, xmin..xmax, axes=boxed);

end proc:

PlotPhi2D_QCX := proc(a, n)
    description "Plot refinable function with quincunx dilation":
    local MM, X1, X2, Y:
    MM := Matrix([[1, 1], [1, -1]]):

    X1, X2, Y := PlotPhi2D(a, MM, n):
    return X1, X2, Y:
end proc:

PlotPsi2D_QCX := proc(a, b, n)
    description "Plot wavelet function with quincunx dilation":
    local MM, X1, X2, Y:
    MM := Matrix([[1, 1], [1, -1]]):

    X1, X2, Y := PlotPsi2D(a, b, MM, n):
    return X1, X2, Y:
end proc:

PlotPhi2D_I2 := proc(a, n)
    description "Plot refinable function with 2I_2 dilation":
    local MM, X1, X2, Y:
    MM := Matrix([[2, 0], [0, 2]]):

    X1, X2, Y := PlotPhi2D(a, MM, n):
    return X1, X2, Y:
end proc:


PlotPsi2D_I2 := proc(a, b, n)
    description "Plot wavelet function with 2I_2 dilation":
    local MM, X1, X2, Y:
    MM := Matrix([[2, 0], [0, 2]]):

    X1, X2, Y := PlotPsi2D(a, b, MM, n):
    return X1, X2, Y:
end proc:



############# Export function ##################
ExportFunctionPlot := proc(p::evaln, pname)
    local name, place, opts:
    name := cat(pname, ".eps"):
    opts := `landscape,width=768,height=768,noborder,axes=boxed,color="Black"`:
    plotsetup('eps', 'plotoutput'=name, 'plotoptions'=opts):
    print( plots:-display( eval(p), 'axesfont' = [ TIMES, 30 ],
                        'labelfont' = [ TIMES, ROMAN, 30] ) ):
    plotsetup(default):
end proc:
