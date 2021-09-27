/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

print "Running self check \"AlgCheRes\"";
zeit := Cputime();

//============================================================================
//with this test I found this strange Magma bug with scalar multiplication in group algebras.
W:=ComplexReflectionGroup(2,1,2); c:=CherednikParameter(W,[1,5]); H:=RestrictedRationalCherednikAlgebra(W,c);
Initialize(~H);
M:=MatrixAlgebra(H);
R:=RModule(M);
M0:=MatrixAlgebraOfDegreeZeroPart(H);
R0 := RModule(M0);

P<z>:=PolynomialRing(Rationals());
K:=NumberField(z^2+1);
W:=ChangeRing(ComplexReflectionGroup(2,1,2),K); c:=CherednikParameter(W,[1,5]); H:=RestrictedRationalCherednikAlgebra(W,c);
Initialize(~H);
MK:=MatrixAlgebra(H);
RK:=RModule(MK);
MK0:=MatrixAlgebraOfDegreeZeroPart(H);
RK0 := RModule(MK0);
assert ActionGenerators(R0) eq ActionGenerators(RK0);

print "Time: "*Sprint(Cputime(zeit));


quit;
