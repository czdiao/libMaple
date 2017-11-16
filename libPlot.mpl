#########################################################
#
#   Tool functions to plot z domain Laurent polynomials
#
#   Define z=exp(-I\xi)
#
#   This library relies on libLpoly.mpl
#
#   Chenzhe
#   Oct, 2016
#

with(LinearAlgebra):
with(plots):
with(PolynomialTools):



PlotPolyZreal := proc(poly)
description "Plot real part value of z-domain Polynomial poly in -Pi..Pi";
## poly should be in z-domain

	local t, polyt:
	polyt := eval(simplify(poly), z = exp(-I*xi)):
	plot([Re(polyt)], xi = -Pi..Pi);
end proc:


PlotPolyZ := proc(poly)
description "Plot abs value of z-domain Polynomial poly in -Pi..Pi";
## poly should be in z-domain

	local xi, polyt:
	polyt := eval(simplify(poly), z = exp(-I*xi)):
	plot([abs(polyt)], xi = -Pi..Pi, color="Black");
end proc:

PlotPolyZn := proc(polylist, legendlist)
    description "Plot multiple z-domain Laurent Polynomials simultaneously":
    # input should be a list: [p1, p2, ..., pn]
	local polyt, n, iter:

	n := nops(polylist):
    polyt := eval(simplify(polylist), z=exp(-I*xi)):
    plot(evalf(abs~(polyt)), xi=-Pi..Pi, legend = legendlist, linestyle=[seq(i, i=1..n)] );


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

    local xval, yval, iter, xmin, xmax, k:
    yval:= 1:   # choose initial function to be delta in cascade algorithm

    for iter from 1 to n do
        yval:= 2 * a * eval(yval, z=z^2):
        yval:= collect(yval, z):
    od:

    xmin:= ldegree(yval, z):
    xmax:= degree(yval, z):
    xval:= [seq( k/2^n, k=xmin..xmax )]:

    yval:=CoefficientVector(collect(yval/z^xmin, z), z):
    plot(xval, yval);

end proc:

PlotPsi := proc(a, b, n)
    description "Plot scalar wavelet function":

    local xval, yval, iter, xmin, xmax, k:
    yval:= b:   

    for iter from 1 to (n-1) do
        yval:= 2 * a * eval(yval, z=z^2):
        yval:= collect(yval, z):
    od:

    xmin:= ldegree(yval, z):
    xmax:= degree(yval, z):
    xval:= [seq( k/2^n, k=xmin..xmax )]:

    yval:=CoefficientVector(collect(yval/z^xmin, z), z):
    plot(xval, yval);
end proc:


#############
PlotPolyZD2 := proc(a)
description "Plot a 2d mask in frequency domain, a is a polynomial of z[1] and z[2]";
	local xi1, xi2, polyt:
	polyt := eval(simplify(a), {z[1] = exp(-I*xi1), z[2] = exp(-I*xi2)}):
	plot3d([abs(polyt)], xi1 = -Pi..Pi, xi2= -Pi..Pi);
end proc:

