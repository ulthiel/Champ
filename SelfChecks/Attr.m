/*
    CHAMP (CHerednik Algebra Magma Package)
    Copyright (C) 2013, 2014 Ulrich Thiel
    Licensed under GNU GPLv3, see COPYING.
    thiel@mathematik.uni-stuttgart.de
*/

//==========================================================================
/*
    Check attribute assignment without reference
*/
print "Running self check \"Attr\"";

Attach("Attr.i.m");

AddAttribute(AlgChtr, "MyAttribute");
G := ShephardTodd(4);
C := CharacterRing(G);
CheckAttributeAssignment(C);
if not assigned C`MyAttribute then
    error "Attribute assignment without reference not working for AlgChtr.";
end if;

AddAttribute(AlgFr, "MyAttribute");
A := FreeAlgebra(Rationals(),3);
CheckAttributeAssignment(A);
if not assigned A`MyAttribute then
    error "Attribute assignment without reference not working for AlgFr.";
end if;

AddAttribute(AlgFP, "MyAttribute");
A := FreeAlgebra(Rationals(),3);
I := ideal<A|A.1>;
B := A/I;
CheckAttributeAssignment(B);
if not assigned B`MyAttribute then
    error "Attribute assignment without reference not working for AlgFP.";
end if;

AddAttribute(AlgGrp, "MyAttribute");
R := GroupAlgebra(Rationals(), ShephardTodd(4));
CheckAttributeAssignment(R);
if not assigned R`MyAttribute then
    error "Attribute assignment without reference not working for AlgGrp.";
end if;

AddAttribute(AlgMat, "MyAttribute");
A := MatrixAlgebra(Rationals(),3);
CheckAttributeAssignment(A);
if not assigned A`MyAttribute then
    error "Attribute assignment without reference not working for AlgMat.";
end if;

AddAttribute(FldRat, "MyAttribute");
K := Rationals();
CheckAttributeAssignment(K);
if not assigned K`MyAttribute then
    error "Attribute assignment without reference not working for FldRat.";
end if;

AddAttribute(FldCyc, "MyAttribute");
K := CyclotomicField(7);
CheckAttributeAssignment(K);
if not assigned K`MyAttribute then
    error "Attribute assignment without reference not working for FldCyc.";
end if;

AddAttribute(FldFin, "MyAttribute");
K := GF(8);
CheckAttributeAssignment(K);
if not assigned K`MyAttribute then
    error "Attribute assignment without reference not working for FldFin.";
end if;

AddAttribute(FldFunRat, "MyAttribute");
K := RationalFunctionField(Rationals(),3);
CheckAttributeAssignment(K);
if not assigned K`MyAttribute then
    error "Attribute assignment without reference not working for FldFunRat.";
end if;

AddAttribute(FldNum, "MyAttribute");
P<X> := PolynomialRing(Rationals());
K := NumberField(X^2+2);
CheckAttributeAssignment(K);
if not assigned K`MyAttribute then
    error "Attribute assignment without reference not working for FldNum.";
end if;

AddAttribute(GrpAb, "MyAttribute");
G := AbelianGroup([2,3]);
CheckAttributeAssignment(G);
if not assigned G`MyAttribute then
    error "Attribute assignment without reference not working for GrpAb.";
end if;

//Does not work from 2.20-5 any more!?
AddAttribute(GrpAuto, "MyAttribute");
G := AutomorphismGroup(SymmetricGroup(3));
CheckAttributeAssignment(G);
if not assigned G`MyAttribute then
    error "Attribute assignment without reference not working for GrpAuto.";
end if;


AddAttribute(GrpFP, "MyAttribute");
G := Group<X|X^2>;
CheckAttributeAssignment(G);
if not assigned G`MyAttribute then
    error "Attribute assignment without reference not working for GrpFP.";
end if;


AddAttribute(GrpMat, "MyAttribute");
G := ShephardTodd(4);
CheckAttributeAssignment(G);
if not assigned G`MyAttribute then
    error "Attribute assignment without reference not working for GrpMat.";
end if;

AddAttribute(GrpPC, "MyAttribute");
G := PCGroup(ShephardTodd(4));
CheckAttributeAssignment(G);
if not assigned G`MyAttribute then
    error "Attribute assignment without reference not working for GrpPC.";
end if;

AddAttribute(GrpPerm, "MyAttribute");
G := SymmetricGroup(5);
CheckAttributeAssignment(G);
if not assigned G`MyAttribute then
    error "Attribute assignment without reference not working for GrpPerm.";
end if;

AddAttribute(GrpRWS, "MyAttribute");
G := RWSGroup(Group<X|X^2>);
CheckAttributeAssignment(G);
if not assigned G`MyAttribute then
    error "Attribute assignment without reference not working for GrpRWS.";
end if;

AddAttribute(HomGrp, "MyAttribute");
G := SymmetricGroup(3);
f := hom<G->G | [G.1, G.2]>;
CheckAttributeAssignment(f);
if not assigned f`MyAttribute then
    error "Attribute assignment without reference not working for HomGrp.";
end if;

AddAttribute(Map, "MyAttribute");
f := map<{1,2}->{1,2}|[<1,1>,<2,2>]>;
CheckAttributeAssignment(f);
if not assigned f`MyAttribute then
    error "Attribute assignment without reference not working for Map.";
end if;

AddAttribute(ModFld, "MyAttribute");
V := KModule(Rationals(),3);
CheckAttributeAssignment(V);
if not assigned V`MyAttribute then
    error "Attribute assignment without reference not working for ModFld.";
end if;

AddAttribute(ModGrp, "MyAttribute");
//G := SymmetricGroup(3);
//Irr := AbsolutelyIrreducibleModules(G, Rationals());
//M := Irr[3];
M := PermutationModule(SymmetricGroup(4), Rationals());
CheckAttributeAssignment(M);
if not assigned M`MyAttribute then
    error "Attribute assignment without reference not working for ModGrp.";
end if;

AddAttribute(ModTupFld, "MyAttribute");
V := VectorSpace(Rationals(),3);
CheckAttributeAssignment(V);
if not assigned V`MyAttribute then
    error "Attribute assignment without reference not working for ModTupFld.";
end if;

AddAttribute(ModRng, "MyAttribute");
V := RModule([IdentityMatrix(Rationals(),2)]);
CheckAttributeAssignment(V);
if not assigned V`MyAttribute then
    error "Attribute assignment without reference not working for ModRng.";
end if;

AddAttribute(RngCyc, "MyAttribute");
K := CyclotomicField(3);
O := RingOfIntegers(K);
CheckAttributeAssignment(O);
if not assigned O`MyAttribute then
    error "Attribute assignment without reference not working for RngCyc.";
end if;

AddAttribute(RngOrd, "MyAttribute");
P<X> := PolynomialRing(Rationals());
K := NumberField(X^2+2);
O := RingOfIntegers(K);
CheckAttributeAssignment(O);
if not assigned O`MyAttribute then
    error "Attribute assignment without reference not working for RngOrd.";
end if;

AddAttribute(RngMPol, "MyAttribute");
P := PolynomialRing(Rationals(),3);
CheckAttributeAssignment(P);
if not assigned P`MyAttribute then
    error "Attribute assignment without reference not working for RngMPol.";
end if;

AddAttribute(RngSerLaur, "MyAttribute");
P := LaurentSeriesRing(Rationals(),3);
CheckAttributeAssignment(P);
if not assigned P`MyAttribute then
    error "Attribute assignment without reference not working for RngSerLaur.";
end if;

//==========================================================================
/*
    Check independent attribute assignment
*/

AddAttribute(AlgFr, "MyAttribute");
A1 := FreeAlgebra(Rationals(),2);
A2 := FreeAlgebra(Rationals(),2);
A1`MyAttribute := 1;
A2`MyAttribute := 0;
if A1`MyAttribute eq A2`MyAttribute then
    error "Independent attribute assignment not working for AlgFr.";
end if;


AddAttribute(GrpMat, "MyAttribute");
G1 := MatrixGroup<2, Rationals() | Matrix(Rationals(),2,[0,1,1,0])>;
G2 := MatrixGroup<2, Rationals() | Matrix(Rationals(),2,[0,1,1,0])>;
G1`MyAttribute := 1;
G2`MyAttribute := 0;
if G1`MyAttribute eq G2`MyAttribute then
    error "Independent attribute assignment not working for GrpMat.";
end if;


AddAttribute(ModFld, "MyAttribute");
V1 := KModule(Rationals(),3);
V2 := KModule(Rationals(),3);
V1`MyAttribute := 1;
V2`MyAttribute := 0;
if V1`MyAttribute eq V2`MyAttribute then
    error "Independent attribute assignment not working for ModTupFld (VectorSpace).";
end if;

AddAttribute(ModRng, "MyAttribute");
V1 := RModule([IdentityMatrix(Rationals(),2)]);
V2 := RModule([IdentityMatrix(Rationals(),2)]);
V1`MyAttribute := 1;
V2`MyAttribute := 0;
if V1`MyAttribute eq V2`MyAttribute then
    error "Independent attribute assignment not working for ModRng.";
end if;


AddAttribute(ModTupFld, "MyAttribute");
V1 := VectorSpace(Rationals(),3);
V2 := VectorSpace(Rationals(),3);
V1`MyAttribute := 1;
V2`MyAttribute := 0;
if V1`MyAttribute ne V2`MyAttribute then
    error "Independent attribute assignment not working for ModTupFld (VectorSpace).";
end if;

AddAttribute(ModTupFld, "MyAttribute");
V1 := KSpace(Rationals(),3);
V2 := KSpace(Rationals(),3);
V1`MyAttribute := 1;
V2`MyAttribute := 0;
if V1`MyAttribute ne V2`MyAttribute then
    error "Independent attribute assignment not working for ModTupFld (KSpace).";
end if;

AddAttribute(RngMPol, "MyAttribute");
P1 := PolynomialRing(Rationals(),2);
P2 := PolynomialRing(Rationals(),2);
P1`MyAttribute := 1;
P2`MyAttribute := 0;
if P1`MyAttribute eq P2`MyAttribute then
    error "Independent attribute assignment not working for RngMPol.";
end if;


quit;
