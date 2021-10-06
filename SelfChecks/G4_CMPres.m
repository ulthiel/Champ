// Checking the presentation of the CM space for G4
// (this is in the paper with Bonnafe)

print "Running self check \"G4_CMPres\"";
zeit := Cputime();

W := ComplexReflectionGroup(4);
H:=RationalCherednikAlgebra(W,0);
CenterGenerators(~H);
Zpresentation:=CenterPresentation(H);

R := H`BaseRing;
k11 := R.1;
k12 := R.2;

K1 := -1/3*k11 + 2/3*k12;
K2 := 2/3*k11 - 1/3*k12;

X1 := H`CenterGenerators[2];
X2 := H`CenterGenerators[6];
Y1 := H`CenterGenerators[5];
Y2 := H`CenterGenerators[8];
eu := H`CenterGenerators[1];
A  := H`CenterGenerators[3];
B  := H`CenterGenerators[4];
C :=  H`CenterGenerators[7];

assert eu eq EulerElement(H);

Zdegree:=func<f|Bidegree(f)[2]-Bidegree(f)[1]>;

assert Zdegree(X1) eq 4;
assert Zdegree(X2) eq 6;
assert Zdegree(Y1) eq -4;
assert Zdegree(Y2) eq -6;
assert Zdegree(eu) eq 0;
assert Zdegree(A) eq 2;
assert Zdegree(B) eq -2;
assert Zdegree(C) eq 0;


// Check relations in the Cherednik algebra
print "Checking relation 1";
assert 2*Y1*X1 - 15*eu^4 + 702*(K1^2 + K1*K2 + K2^2)*eu^2 + 12*C*eu + 2592*K1*K2*(K1+K2)*eu + B*A eq H`Zero;

print "Checking relation 2";
assert 2*Y2*X1 + 3*A*eu*Y1 - 9*B*eu^3 + 378*(K1^2 + K1*K2 + K2^2)*B*eu + 4*C*B eq H`Zero;

print "Checking relation 3";
assert 9*eu^3*X1 - 324*(K1^2 + K1*K2 + K2^2)*eu*X1 - 8*C*X1 + 2*B*X2 - 3*A^2*eu eq H`Zero;

print "Checking relation 4";
assert 3*B*eu*X1 + 2*Y1*X2 - 9*A*eu^3 + 378*(K1^2 + K1*K2 + K2^2)*A*eu + 4*C*A eq H`Zero;

print "Checking relation 5";
assert 9*eu^3*Y1 - 324*(K1^2 + K1*K2 + K2^2)*eu*Y1 - 8*C*Y1 + 2*A*Y2 - 3*B^2*eu eq H`Zero;

print "Checking relation 6";
assert 2*B*X1^2 - 3*A*eu^2*X1 + 144*(K1^2 + K1*K2 + K2^2)*A*X1 + 10*eu^3*X2 - 468*(K1^2 + K1*K2 + K2^2)*eu*X2 \
  - 8*C*X2 - 1728*K1*K2*(K1+K2)*X2 + 2*A*Y1^2 - 3*B*eu^2*Y1 + 144*(K1^2 + K1*K2 + K2^2)*B*Y1 \
  + 10*eu^3*Y2 - 468*(K1^2+K1*K2+K2^2)*eu*Y2 - 8*C*Y2 - 1728*K1*K2*(K1+K2)*Y2 - A^3 - B^3 eq H`Zero;

print "Checking relation 7";
assert 9*eu^2*Y1*X1 + 2*Y2*X2 - 27*eu^6 + 11664*K1*K2*(K1+K2)*eu^3 + 61236*(K1^2+K1*K2+K2^2)^2*eu^2 \
  + 2160*(K1^2+K1*K2+K2^2)*C*eu + 16*C^2 eq H`Zero;

print "Checking relation 8";
assert 2*A*Y1^2 - 3*B*eu^2*Y1 + 144*(K1^2 + K1*K2 + K2^2)*B*Y1 + 10*eu^3*Y2 - 468*(K1^2 + K1*K2 + K2^2)*eu*Y2 - 8*C*Y2 \
  - 1728*K1*K2*(K1+K2)*Y2 - B^3 eq H`Zero;

print "Checking relation 9";
assert 60*eu^2*Y1*X1 + 1944*(K1^2+K1*K2+K2^2)*Y1*X1 + 5*B^2*X1 + 10*Y2*X2 + 5*A^2*Y1 - 360*eu^6 + 280*C*eu^3 \
  + 97200*K1*K2*(K1+K2)*eu^3 + 798984*(K1^2+K1*K2+K2^2)^2*eu^2 + 14544*(K1^2+K1*K2+K2^2)*C*eu \
  + 1819584*K1*K2*(K1+K2)*(K1^2 + K1*K2 + K2^2)*eu + 1332*(K1^2+K1*K2+K2^2)*B*A \
  - 17280*K1*K2*(K1+K2)*C eq H`Zero;


// Check whether relations match
X1 := H`CenterSpace.2;
X2 := H`CenterSpace.6;
Y1 := H`CenterSpace.5;
Y2 := H`CenterSpace.8;
eu := H`CenterSpace.1;
A  := H`CenterSpace.3;
B  := H`CenterSpace.4;
C :=  H`CenterSpace.7;

rels := H`CenterPresentation;

// Relation 1
assert -rels[1] eq 2*Y1*X1 - 15*eu^4 + 702*(K1^2 + K1*K2 + K2^2)*eu^2 + 12*C*eu + 2592*K1*K2*(K1+K2)*eu + B*A;

// Relation 2
assert 1/2*rels[4] eq 2*Y2*X1 + 3*A*eu*Y1 - 9*B*eu^3 + 378*(K1^2 + K1*K2 + K2^2)*B*eu + 4*C*B;

// Relation 3
assert -2*rels[5] eq 9*eu^3*X1 - 324*(K1^2 + K1*K2 + K2^2)*eu*X1 - 8*C*X1 + 2*B*X2 - 3*A^2*eu;

// Relation 4
assert 1/2*rels[2] eq 3*B*eu*X1 + 2*Y1*X2 - 9*A*eu^3 + 378*(K1^2 + K1*K2 + K2^2)*A*eu + 4*C*A;

// Relation 5
assert 2*rels[3] eq 9*eu^3*Y1 - 324*(K1^2 + K1*K2 + K2^2)*eu*Y1 - 8*C*Y1 + 2*A*Y2 - 3*B^2*eu;

// Relation 6
assert rels[6] eq 2*B*X1^2 - 3*A*eu^2*X1 + 144*(K1^2 + K1*K2 + K2^2)*A*X1 + 10*eu^3*X2 - 468*(K1^2 + K1*K2 + K2^2)*eu*X2 \
  - 8*C*X2 - 1728*K1*K2*(K1+K2)*X2 + 2*A*Y1^2 - 3*B*eu^2*Y1 + 144*(K1^2 + K1*K2 + K2^2)*B*Y1 \
  + 10*eu^3*Y2 - 468*(K1^2+K1*K2+K2^2)*eu*Y2 - 8*C*Y2 - 1728*K1*K2*(K1+K2)*Y2 - A^3 - B^3;

// Relation 7
assert -1/4*(rels[7] - 8*rels[9]) eq 9*eu^2*Y1*X1 + 2*Y2*X2 - 27*eu^6 + 11664*K1*K2*(K1+K2)*eu^3 + 61236*(K1^2+K1*K2+K2^2)^2*eu^2 \
  + 2160*(K1^2+K1*K2+K2^2)*C*eu + 16*C^2;

// Relation 8
assert rels[8] + rels[6] eq 2*A*Y1^2 - 3*B*eu^2*Y1 + 144*(K1^2 + K1*K2 + K2^2)*B*Y1 + 10*eu^3*Y2 - 468*(K1^2 + K1*K2 + K2^2)*eu*Y2 - 8*C*Y2 \
  - 1728*K1*K2*(K1+K2)*Y2 - B^3;

// Relation 9
assert 10*rels[9] eq 60*eu^2*Y1*X1 + 1944*(K1^2+K1*K2+K2^2)*Y1*X1 + 5*B^2*X1 + 10*Y2*X2 + 5*A^2*Y1 - 360*eu^6 + 280*C*eu^3 \
  + 97200*K1*K2*(K1+K2)*eu^3 + 798984*(K1^2+K1*K2+K2^2)^2*eu^2 + 14544*(K1^2+K1*K2+K2^2)*C*eu \
  + 1819584*K1*K2*(K1+K2)*(K1^2 + K1*K2 + K2^2)*eu + 1332*(K1^2+K1*K2+K2^2)*B*A \
  - 17280*K1*K2*(K1+K2)*C;

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();

quit;
