/*
    CHAMP (CHerednik Algebra Magma Package)
    Copyright (C) 2013, 2014 Ulrich Thiel
    Licensed under GNU GPLv3, see COPYING.
    thiel@mathematik.uni-stuttgart.de
*/

/*
    Testing if Verma modules work for restricted generic rational Cherednik algebra of G5.
*/
print "Running self check \"G5_Verma\"";
zeit := Cputime();

G:=ExceptionalComplexReflectionGroup(5); c:=CherednikParameter(G); Representations(~G,0);
DiagonalizeRepresentations(~G,0,1); LiftRepresentationsToCommonBaseField(~G,0);
V:=[* VermaModule(G,c,G`Representations[0][i]) : i in [19,20,21] *];

for i:=1 to 3 do
    assert IsModuleForRRCA(G,c,V[i]);
end for;

res,L,D:=HeadsOfLocalModules(V:pRange:={10^5..10^6});

for i:=1 to 3 do
    assert IsModuleForRRCA(G,c,L[i]);
end for;

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();

quit;
