/*
    CHAMP (CHerednik Algebra Magma Package)
    Copyright (C) 2013, 2014 Ulrich Thiel
    Licensed under GNU GPLv3, see COPYING.
    thiel@mathematik.uni-stuttgart.de
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
V:=[* GradedModuleOld(VermaModule(H,G`Representations[0][i])) : i in [1..7] *];

for i:=1 to 7 do
    assert IsModule(H,V[i]);
end for;

res,L,D:=HeadsOfLocalModules(V:pRange:={10^5..10^6});

for i:=1 to 7 do
    assert IsModule(H,L[i]);
end for;

//reversed ones
V:=[* VermaModule(G,c,G`Representations[0][i] : Reversed:=true) : i in [1..7] *];

for i:=1 to 7 do
    assert IsModuleForRRCA(G,c,V[i]);
end for;

res,L,D:=HeadsOfLocalModules(V:pRange:={10^5..10^6});

for i:=1 to 7 do
    assert IsModuleForRRCA(G,c,L[i]);
end for;

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();

quit;
