/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Testing if Verma modules work for restricted generic rational Cherednik algebra of G4.
*/
print "Running self check \"G4_Verma\"";
zeit := Cputime();

G:=ExceptionalComplexReflectionGroup(4);
c:=CherednikParameter(G);
Representations(~G,0);
H:=RestrictedRationalCherednikAlgebra(G);
StandardModules(~H);
V:=[*H`StandardModules[i] : i in [1..#H`StandardModules]*];

for i:=1 to 7 do
    assert IsModule(H,V[i]);
end for;

res,L,D:=HeadsOfLocalModules(V:pRange:={10^5..10^6});

for i:=1 to 7 do
    assert IsModule(H,L[i]);
end for;

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();

quit;
