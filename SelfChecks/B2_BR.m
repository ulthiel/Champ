/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Computation in [BR13], Section 19.

    This is a very serious test of AlgChe.
*/
//StartProfile();

print "Running self check \"B2_BR\"";

G:=CHAMP_GetFromDB("ReflectionGroups/B2_BR", "GrpMat");

C:=CherednikParameter(G : Type:="BR", Rational:=false);

H:=RationalCherednikAlgebra(G,0,C);


zeit := Cputime();
/*
x:=GetVariable(H, "y1"); //notation is not great actually
y:=GetVariable(H, "y2");
X:=GetVariable(H, "x1");
Y:=GetVariable(H, "x2");
s:=GetVariable(H, "g1");
t:=GetVariable(H, "g2");
*/
x := H.3;
y := H.4;
X := H.5;
Y := H.6;
s := H.1;
t := H.2;

//19.1

//note that everything is opposite!
w:=t*s;
wprime:=s*t;
sprime:=t*s*t;
tprime:=s*t*s;
w0:=sprime*s;

A:=C(1)*(-1/2); //this is correct! we get C_s in this way!
B:=C(2)*(-1/2);

//A:=C(1);
//B:=C(2);

//19.1.2
assert X*x-x*X eq -A*(s+sprime) - 2*B*t;
assert Y*x-x*Y eq A*(s-sprime);
assert X*y-y*X eq A*(s-sprime);
assert Y*y-y*Y eq -A*(s+sprime)-2*B*tprime;

//19.1.3
assert X^2*x-x*X^2 eq -X*(s+sprime)*A-Y*(s-sprime)*A;
assert Y*X*x-x*Y*X eq -Y*t*B*2;
assert Y^2*x-x*Y^2 eq X*(s+sprime)*A + Y*(s-sprime)*A;
assert X^2*y-y*X^2 eq Y*(s+sprime)*A + X*(s-sprime)*A;
assert Y*X*y-y*Y*X eq -X*tprime*B*2;
assert Y^2*y-y*Y^2 eq -Y*(s+sprime)*A - X*(s-sprime)*A;

//19.3
sigma:=x^2+y^2;
pi:=y^2*x^2;
Sigma:=X^2+Y^2;
Pi:=Y^2*X^2;

//19.4
eu:=X*x+Y*y+A*(s+sprime)+B*(t+tprime);
euprime:=Y*X*(Y*x+X*y) - Y*X*(s-sprime)*A + Y^2*t*B + X^2*tprime*B;
euprimeprime:=(Y*x+X*y)*y*x - (s-sprime)*y*x*A + t*y^2*B + tprime*x^2*B;
delta:=Y*X*y*x + X*tprime*x*B + Y*t*y*B + (1+w0)*B^2 + (w+wprime)*B*A;

//19.4.1
assert IsCentral(eu);
assert IsCentral(euprime);
assert IsCentral(euprimeprime);
assert IsCentral(delta);

//19.4.2
assert euprime*eu eq Pi*sigma + delta*Sigma;   //Z1
assert euprimeprime*eu eq pi*Sigma + delta*sigma;  //Z2
assert euprime*delta eq euprimeprime*Pi + eu*Sigma*B^2;    //Z3
assert euprimeprime*delta eq euprime*pi + eu*sigma*B^2;    //Z4
assert delta^2 eq Pi*pi + eu^2*B^2; //Z5
assert euprime^2 eq (4*delta - eu^2 + Sigma*sigma + 4*A^2 - 4*B^2)*Pi + Sigma^2*B^2;   //Z6
assert euprimeprime^2 eq (4*delta - eu^2 + Sigma*sigma + 4*A^2 - 4*B^2)*pi + sigma^2*B^2;   //Z7
assert euprimeprime*euprime eq (4*delta - eu^2 + Sigma*sigma + 4*A^2 - 4*B^2)*delta - Sigma*sigma*B^2; //Z8
assert (4*delta - eu^2 + Sigma*sigma + 4*A^2 - 4*B^2)*eu eq euprime*sigma + euprimeprime*Sigma;    //Z9

//19.4.5
//the following is the minimal polynomial of eu.
assert eu^8 - 2*eu^6*(Sigma*sigma + 4*A^2 + 4*B^2) + eu^4*( Sigma^2*sigma^2 + 2*(Pi*sigma^2 + pi*Sigma^2 - 8*Pi*pi) + 8*Sigma*sigma*(A^2+B^2) + 16*(A^2-B^2)^2) - 2*eu^2*( (Pi*sigma^2 + pi*Sigma^2)*(Sigma*sigma + 4*A^2 - 4*B^2) - 8*Pi*pi*Sigma*sigma + Sigma^2*sigma^2*B^2*2) + (Pi*sigma^2 - pi*Sigma^2)^2 eq Zero(H);

//Poisson brackets (Part II, 2.1)
deltaplus := 12*delta - 2*eu^2 + 2*Sigma*sigma + 8*(A^2-B^2);

assert PoissonBracket(sigma,sigma) eq Zero(H);
assert PoissonBracket(sigma,pi) eq Zero(H);
assert PoissonBracket(Sigma,sigma) eq 4*eu;
assert PoissonBracket(Pi,sigma) eq 4*euprime;
assert PoissonBracket(eu,sigma) eq 2*sigma;
assert PoissonBracket(delta,sigma) eq 2*euprimeprime;
assert PoissonBracket(euprime,sigma) eq deltaplus;
assert PoissonBracket(euprimeprime,sigma) eq 4*pi;

assert PoissonBracket(sigma,pi) eq Zero(H);
assert PoissonBracket(pi,pi) eq Zero(H);
assert PoissonBracket(Sigma,pi) eq 4*euprimeprime;
assert PoissonBracket(Pi,pi) eq 4*eu*delta - 8*B^2*eu;
assert PoissonBracket(eu,pi) eq 4*pi;
assert PoissonBracket(delta,pi) eq 2*pi*eu;
assert PoissonBracket(euprime,pi) eq 4*delta*sigma - 4*B^2*sigma + 2*Sigma*pi;
assert PoissonBracket(euprimeprime,pi) eq 2*pi*sigma;

assert PoissonBracket(sigma,Sigma) eq -4*eu;
assert PoissonBracket(pi,Sigma) eq -4*euprimeprime;
assert PoissonBracket(Sigma,Sigma) eq Zero(H);
assert PoissonBracket(Pi,Sigma) eq Zero(H);
assert PoissonBracket(eu,Sigma) eq -2*Sigma;
assert PoissonBracket(delta,Sigma) eq -2*euprime;
assert PoissonBracket(euprime,Sigma) eq -4*Pi;
assert PoissonBracket(euprimeprime,Sigma) eq -deltaplus;

assert PoissonBracket(sigma,Pi) eq -4*euprime;
assert PoissonBracket(pi,Pi) eq -4*eu*delta + 8*B^2*eu;
assert PoissonBracket(Sigma,Pi) eq Zero(H);
assert PoissonBracket(Pi,Pi) eq Zero(H);
assert PoissonBracket(eu,Pi) eq -4*Pi;
assert PoissonBracket(delta,Pi) eq -2*eu*Pi;
assert PoissonBracket(euprime,Pi) eq -2*Pi*Sigma;
assert PoissonBracket(euprimeprime,Pi) eq -4*delta*Sigma + 4*B^2*Sigma - 2*Pi*sigma;

assert PoissonBracket(sigma,eu) eq -2*sigma;
assert PoissonBracket(pi,eu) eq -4*pi;
assert PoissonBracket(Sigma,eu) eq 2*Sigma;
assert PoissonBracket(Pi,eu) eq 4*Pi;
assert PoissonBracket(eu,eu) eq Zero(H);
assert PoissonBracket(delta,eu) eq Zero(H);
assert PoissonBracket(euprime,eu) eq 2*euprime;
assert PoissonBracket(euprimeprime,eu) eq -2*euprimeprime;

assert PoissonBracket(sigma,delta) eq -2*euprimeprime;
assert PoissonBracket(pi,delta) eq -2*eu*pi;
assert PoissonBracket(Sigma,delta) eq 2*euprime;
assert PoissonBracket(Pi,delta) eq 2*eu*Pi;
assert PoissonBracket(eu,delta) eq Zero(H);
assert PoissonBracket(delta,delta) eq Zero(H);
assert PoissonBracket(euprime,delta) eq 2*(Pi*sigma + B^2*Sigma);
assert PoissonBracket(euprimeprime,delta) eq -2*(Sigma*pi + B^2*sigma);

assert PoissonBracket(sigma,euprime) eq -deltaplus;
assert PoissonBracket(pi,euprime) eq -4*delta*sigma + 4*B^2*sigma - 2*Sigma*pi;
assert PoissonBracket(Sigma,euprime) eq 4*Pi;
assert PoissonBracket(Pi,euprime) eq 2*Pi*Sigma;
assert PoissonBracket(eu,euprime) eq -2*euprime;
assert PoissonBracket(delta,euprime) eq -2*(Pi*sigma+B^2*Sigma);
assert PoissonBracket(euprime,euprime) eq Zero(H);
assert PoissonBracket(euprimeprime,euprime) eq -4*B^2*eu - 2*euprime*sigma - 2*euprimeprime*Sigma;

assert PoissonBracket(sigma,euprimeprime) eq -4*pi;
assert PoissonBracket(pi,euprimeprime) eq -2*sigma*pi;
assert PoissonBracket(Sigma,euprimeprime) eq deltaplus;
assert PoissonBracket(Pi,euprimeprime) eq 4*delta*Sigma - 4*B^2*Sigma + 2*Pi*sigma;
assert PoissonBracket(eu,euprimeprime) eq 2*euprimeprime;
assert PoissonBracket(delta,euprimeprime) eq 2*(Sigma*pi + B^2*sigma);
assert PoissonBracket(euprime,euprimeprime) eq 4*B^2*eu + 2*euprime*sigma + 2*euprimeprime*Sigma;
assert PoissonBracket(euprimeprime,euprimeprime) eq Zero(H);

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();


//P := StopProfile();
//ProfileHTMLOutput(P,"test");
//ProfilePrintByTotalTime(P);

quit;
