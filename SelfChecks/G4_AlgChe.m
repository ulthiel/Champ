/*
    CHAMP (CHerednik Algebra Magma Package)
    Copyright (C) 2013, 2014 Ulrich Thiel
    Licensed under GNU GPLv3, see COPYING.
    thiel@mathematik.uni-stuttgart.de
*/
/*
    Testing some relations in the generic rational Cherednik algebra of G4.
*/
print "Running self check \"G4_AlgChe\"";
zeit := Cputime();

G:=ExceptionalComplexReflectionGroup(4);
H:=RationalCherednikAlgebra(G);
eu:=EulerElement(H);
x1:=H.5;
x2:=H.6;
y1:=H.3;
y2:=H.4;
t:=H`tParameter;
assert x1*eu - eu*x1 eq t*x1;
assert x2*eu - eu*x2 eq t*x2;
assert y1*eu - eu*y1 eq -t*y1;
assert y2*eu - eu*y2 eq -t*y2;

H:=RationalCherednikAlgebra(G,0);
eu:=EulerElement(H);
assert IsCentral(eu);
assert IsCentral(eu^2);
assert IsCentral(eu^3);

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();

quit;
