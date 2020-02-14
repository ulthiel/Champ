/*
    CHAMP (CHerednik Algebra Magma Package)
    Copyright (C) 2013, 2014 Ulrich Thiel
    Licensed under GNU GPLv3, see COPYING.
    thiel@mathematik.uni-stuttgart.de
*/
/*
    Computation in [BR13], Chapter 6.
    
    This is a very serious test of AlgChe.
*/

/*
	Timing: 
		* Version 1.5-228: 80.820 on Mac
*/ 

print "Running self check \"G4_BR\"";
zeit := Cputime();

W:=CHAMP_GetFromDB("ReflectionGroups/G4_BR", "GrpMat");

C:=CherednikParameter(W : Type:="BR", Rational:=false);

H:=RationalCherednikAlgebra(W,<0,C> : PoissonModtSquare:=true);

zeta := RootOfUnity(3);

C1:=C(1)/(zeta-1); //this is correct! we get C_s in this way!
C2:=C(2)/(-zeta-2);

K0 := 1/3*zeta*C1 + 1/3*(-zeta - 1)*C2;
K1 := 1/3*C1 + 1/3*C2;
K2 := 1/3*(-zeta - 1)*C1 + 1/3*zeta*C2;

//6.2.A (note, x and y are switched here; seems to be a mistake)
SymplecticDoublingFundamentalInvariants(~W);
A := W`SymplecticDoublingCoordinateAlgebra;
y1 := A.1;
y2 := A.2;
x1 := A.3;
x2 := A.4;

sigma := x1^4 + 4/3*(1+2*zeta)*x2*x1^3 - 2*x2^2*x1^2 - 4/3*(1+2*zeta)*x2^3*x1 + x2^4;
pi := x1^6 + 2*(1+2*zeta)*x2*x1^5 - 5*x2^2*x1^4 - 5*x2^4*x1^2 - 2*(1+2*zeta)*x2^5*x1 + x2^6;
Sigma := y2*y1^3 + (1+2*zeta)*y2^2*y1^2 - y2^3*y1;
Pi := y1^6 + 2*(1+2*zeta)*y2*y1^5 - 5*y2^2*y1^4 - 5*y2^4*y1^2 - 2*(1+2*zeta)*y2^5*y1 + y2^6;

eu0 := x1*y1 + x2*y2;
a0 := x1*y1^3 + (1+2*zeta)*x2*y1^3 + (1+2*zeta)*x1*y2*y1^2 - x2*y2*3*y1^2 - 3*x1*y2^2*y1 - (1+2*zeta)*x2*y2^2*y1 - (1+2*zeta)*x1*y2^3 + x2*y2^3;
b0 := x1^3*y1 + 2*(1+2*zeta)*x2*x1^2*y1 - 3*x2^2*x1*y1 - 3*x2*x1^2*y2 - 2*(1+2*zeta)*x2^2*x1*y2 + x2^3*y2;
delta0 := x2*x1^2*y1^3 + (1+2*zeta)*x2^2*x1*y1^3 - x2^3*y1^3 - x1^3*y2*y1^2 - 3*x2^2*x1*y2*y1^2 - (1+2*zeta)*x2^3*y2*y1^2 - (1+2*zeta)*x1^3*y2^2*y1 + 3*x2*x1^2*y2^2*y1 + x2^3*y2^2*y1 + x1^3*y2^3 + (1+2*zeta)*x2*x1^2*y2^3 - x2^2*x1*y2^3;

assert SymplecticDoublingAction(sigma,W.1) eq sigma;
assert SymplecticDoublingAction(sigma,W.2) eq sigma;

assert SymplecticDoublingAction(pi,W.1) eq pi;
assert SymplecticDoublingAction(pi,W.2) eq pi;

assert SymplecticDoublingAction(Sigma,W.1) eq Sigma;
assert SymplecticDoublingAction(Sigma,W.2) eq Sigma;

assert SymplecticDoublingAction(Pi,W.1) eq Pi;
assert SymplecticDoublingAction(Pi,W.2) eq Pi;

assert SymplecticDoublingAction(eu0,W.1) eq eu0;
assert SymplecticDoublingAction(eu0,W.2) eq eu0;

assert SymplecticDoublingAction(a0,W.1) eq a0;
assert SymplecticDoublingAction(a0,W.2) eq a0;

assert SymplecticDoublingAction(b0,W.1) eq b0;
assert SymplecticDoublingAction(b0,W.2) eq b0;

assert SymplecticDoublingAction(delta0,W.1) eq delta0;
assert SymplecticDoublingAction(delta0,W.2) eq delta0;

assert SequenceToSet(W`SymplecticDoublingFundamentalInvariants) eq {sigma,pi,Sigma,Pi,eu0,a0,b0,delta0};

//6.2.2
assert eu0^4 + 4*(1+2*zeta)*sigma*Sigma - b0*a0 + 3*(1+2*zeta)*delta0*eu0 eq 0;
assert (1+2*zeta)*a0*eu0^3 - 4*b0*eu0*Sigma - delta0*a0 - (1+2*zeta)*Pi*sigma eq 0;
assert a0*pi - b0^2*eu0 + (1+2*zeta)*delta0*sigma eq 0;
assert (1+2*zeta)*b0*eu0^3 - delta0*b0 - 4*pi*Sigma - (1+2*zeta)*a0*eu0*sigma eq 0;
assert a0^2*eu0 - 4*delta0*Sigma - Pi*b0 eq 0;
assert eu0^3*pi + 3*a0*sigma^2 - 3*b0*eu0^2*sigma - b0^3 + 3*(1+2*zeta)*delta0*pi eq 0;
assert 4*(1+2*zeta)*a0*eu0^2*Sigma - a0^3 - 16*b0*Sigma^2 + eu0^3*Pi + 3*(1+2*zeta)*delta0*Pi eq 0;
assert eu0^6 + 4*(1+2*zeta)*eu0^2*sigma*Sigma + 2*(1+2*zeta)*delta0*eu0^3 + delta0^2 - Pi*pi eq 0;
assert (1+2*zeta)*a0^2*sigma + 4*b0^2*Sigma + 4*delta0*eu0^3 - (1+2*zeta)*Pi*pi + 4*(1+2*zeta)*delta0^2 eq 0;

//6.2.B
eu := TruncationInverse(H,eu0);
a := TruncationInverse(H,a0);
b := TruncationInverse(H,b0);
delta := TruncationInverse(H,delta0);

assert TruncationInverse(H,sigma) eq H!sigma;
assert TruncationInverse(H,pi) eq H!pi;
assert TruncationInverse(H,Sigma) eq H!Sigma;
assert TruncationInverse(H,Pi) eq H!Pi;

sigma := H!sigma;
pi := H!pi;
Sigma := H!Sigma;
Pi := H!Pi;

assert IsCentral(eu);
assert IsCentral(a);
assert IsCentral(b);
assert IsCentral(delta);
assert IsCentral(sigma);
assert IsCentral(pi);
assert IsCentral(Sigma);
assert IsCentral(Pi);

//(Z1)
assert eu^4 + 4*(1+2*zeta)*Sigma*sigma - b*a + 3*(1+2*zeta)*delta*eu + 18*eu^2*(K1^2 + K2*K1 + K2^2) + 756*eu*(K2*K1^2 + K2^2*K1) eq Zero(H);

//(Z2)
assert (1+2*zeta)*a*eu^3 - 4*b*eu*Sigma - a*delta - (1+2*zeta)*Pi*sigma - 30*(1+2*zeta)*(K1^2+K1*K2+K2^2)*a*eu - 108*(1+2*zeta)*(K1^2*K2 + K1*K2^2)*a eq Zero(H);

//(Z3)
assert a*pi - b^2*eu + (1+2*zeta)*delta*sigma + 54*(K1^2+K1*K2+K2^2)*eu*sigma - 324*(K1^2*K2+K1*K2^2)*sigma eq Zero(H);

//(Z4)
assert (1+2*zeta)*b*eu^3 - 4*pi*Sigma - (1+2*zeta)*a*eu*sigma - delta*b - 30*(1+2*zeta)*(K1^2 + K1*K2 + K2^2)*b*eu - 108*(1+2*zeta)*(K1^2*K2 + K1*K2^2)*b eq Zero(H);

//(Z5)
assert a^2*eu - 4*delta*Sigma - b*Pi + 72*(1+2*zeta)*(K1^2+K1*K2+K2^2)*eu*Sigma - 432*(1+2*zeta)*(K1^2*K2 + K1*K2^2)*Sigma eq Zero(H);

//(Z6)
assert eu^3*pi + 3*a*sigma^2 - 3*b*eu^2*sigma - b^3 + 3*(1+2*zeta)*delta*pi + 144*(K1^2+K1*K2+K2^2)*b*sigma + 18*(K1^2+K1*K2+K2^2)*eu*pi + 756*(K1^2*K2+K1*K2^2)*pi eq Zero(H);

//(Z7)
assert 4*(1+2*zeta)*a*eu^2*Sigma - a^3 - 16*b*Sigma^2 + eu^3*Pi + 3*(1+2*zeta)*delta*Pi + 18*(K1^2 + K1*K2 + K2^2)*eu*Pi - 192*(1+2*zeta)*(K1^2 + K1*K2 + K2^2)*a*Sigma + 756*(K1^2*K2 + K1*K2^2)*Pi eq Zero(H);

//(Z8)
assert eu^6 + 4*(1+2*zeta)*eu^2*Sigma*sigma + 2*(1+2*zeta)*eu^3*delta - Pi*pi + delta^2 - 36*(K1^2 + K1*K2 + K2^2)*eu^4 + 12*(1+2*zeta)*(K1^2+K1*K2+K2^2)*delta*eu + 216*(1+2*zeta)*(K1^2*K2 + K1*K2^2)*delta + 1080*(K1^2*K2 + K1*K2^2)*eu^3 + 1620*((K1^2+K1*K2+K2^2))^2*eu^2 - 3888*(K1^2*K2+K1*K2^2)*(K1^2+K1*K2+K2^2)*eu - 34992*(K1^2*K2 + K1*K2^2)^2 eq Zero(H);

//(Z9)
assert (1+2*zeta)*a^2*sigma + 4*b^2*Sigma + 4*delta*eu^3 - (1+2*zeta)*Pi*pi + 4*(1+2*zeta)*delta^2 - 24*(1+2*zeta)*(K1^2+K1*K2+K2^2)*eu^4 - 288*(K1^2+K1*K2+K2^2)*delta*eu - 576*(K1^2+K1*K2+K2^2)*Sigma*sigma + 2160*(1+2*zeta)*(K1^2+K1*K2+K2^2)^2*eu^2 + 432*(1+2*zeta)*(K1^2*K2 + K1*K2^2)*eu^3 - 864*(K1^2*K2 + K1*K2^2)*delta + 20736*(1+2*zeta)*(K1^2*K2 + K1*K2^2)*(K1^2 + K1*K2 + K2^2)*eu + 46656*(1+2*zeta)*(K1^2*K2 + K1*K2^2)^2 eq Zero(H);

//6.5.A
assert PoissonBracket(Sigma,sigma) eq -4/3*(1+2*zeta)*eu^3 + 4*delta + 24*(1+2*zeta)*(K1^2 + K1*K2 + K2^2)*eu - 144*(1+2*zeta)*(K1^2*K2 + K1*K2^2);

assert PoissonBracket(Sigma,pi) eq -4*(1+2*zeta)*eu^2*b + 2*(1+2*zeta)*a*sigma + 48*(1+2*zeta)*(K1^2 + K1*K2 + K2^2)*b;

assert PoissonBracket(Sigma,eu) eq -4*Sigma;

assert PoissonBracket(Sigma,a) eq -(1+2*zeta)*Pi;

assert PoissonBracket(Sigma,b) eq -2*(1+2*zeta)*a*eu;

assert PoissonBracket(Sigma,delta) eq -a^2 - 72*(1+2*zeta)*(K1^2 + K1*K2 + K2^2)*Sigma;

assert PoissonBracket(Pi,sigma) eq -16*a*eu^2 - (32/3)*(1+2*zeta)*b*Sigma + 192*(K1^2 + K1*K2 + K2^2)*a;

//corrected after talking to cedric
//(1440 eu^3 + 96(1+2\zeta)\delta) was replaced by (1440 eu^3 - 96(1+2\zeta)\delta)
assert PoissonBracket(Pi,pi) eq -24*b*a*eu + 32*(1+2*zeta)*delta*eu^2 + +16*(1+2*zeta)*eu*sigma*Sigma + (K1^2+K1*K2+K2^2)*(1440*eu^3 - 96*(1+2*zeta)*delta)-19008*(K1^2+K1*K2+K2^2)^2*eu + 31104*(K1^2 + K1*K2+K2^2)*(K1^2*K2 + K1*K2^2);

assert PoissonBracket(Pi,eu) eq - 6*Pi;

assert PoissonBracket(Pi,a) eq 32*Sigma^2;

assert PoissonBracket(Pi,b) eq 8*(1+2*zeta)*eu^2*Sigma - 6*a^2 - 384*(1+2*zeta)*(K1^2 + K1*K2 + K2^2)*Sigma;

assert PoissonBracket(Pi,delta) eq -2*(1+2*zeta)*eu^2*Pi + 16*a*eu*Sigma - 12*(1+2*zeta)*(K1^2+K1*K2+K2^2)*Pi;

assert PoissonBracket(sigma,eu) eq 4*sigma;

assert PoissonBracket(sigma,a) eq 8*b*eu;

assert PoissonBracket(sigma,b) eq 4*pi;

assert PoissonBracket(sigma,delta) eq -(4/3)*(1+2*zeta)*b^2 + 72*(1+2*zeta)*(K1^2+K1*K2+K2^2)*sigma;

assert PoissonBracket(pi,eu) eq 6*pi;

assert PoissonBracket(pi,a) eq 6*eu^2*sigma + 6*b^2 - 288*(K1^2+K1*K2+K2^2)*sigma;

assert PoissonBracket(pi,b) eq 6*sigma^2;

assert PoissonBracket(pi,delta) eq 2*(1+2*zeta)*eu^2*pi - 4*(1+2*zeta)*b*eu*sigma + 12*(1+2*zeta)*(K1^2+K1*K2+K2^2)*pi;

assert PoissonBracket(a,eu) eq -2*a;

assert PoissonBracket(a,b) eq -2*eu^3 - 10*(1+2*zeta)*delta - 252*(K1^2+K1*K2+K2^2)*eu - 216*(K1^2*K2 + K1*K2^2);

assert PoissonBracket(a,delta) eq -2*(1+2*zeta)*a*eu^2 + 16*b*Sigma + 60*(1+2*zeta)*(K1^2+K1*K2+K2^2)*a;

assert PoissonBracket(b,eu) eq 2*b;

assert PoissonBracket(b,delta) eq 2*(1+2*zeta)*b*eu^2 - 4*(1+2*zeta)*a*sigma - 60*(1+2*zeta)*(K1^2+K1*K2+K2^2)*b;

assert PoissonBracket(delta,eu) eq Zero(H);

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
print "Memory usage: "*Sprint(GetMemoryUsageInMB());
IndentPop();


quit;