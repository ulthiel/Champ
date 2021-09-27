/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Testing if Verma modules work for restricted generic rational Cherednik algebra of G5. This takes some time.
*/
print "Running self check \"G5_Verma\"";
zeit := Cputime();

G1:=ComplexReflectionGroup(5);
G := ChangeRing(G1, CyclotomicField(12));
G`DBDir := G1`DBDir;
Representations(~G,0);
LiftRepresentationsToCommonBaseField(~G,0);
H:=RestrictedRationalCherednikAlgebra(G);
StandardModules(~H);
V:=[*H`StandardModules[i] : i in [19,20,21]*];

for i:=1 to 3 do
    assert IsModule(H,V[i]);
end for;

res,L,D:=HeadsOfLocalModules(V:pRange:={10^5..10^6});

for i:=1 to 3 do
    assert IsModule(H,L[i]);
end for;

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();

quit;
