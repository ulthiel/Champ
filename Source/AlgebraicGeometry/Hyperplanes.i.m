/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Some basic intrinsics dealing with hyperplanes.
*/

//============================================================================
intrinsic Subspace(f::RngMPolElt) -> ModTupFld
{
	If f is a linear polynomial in several variables, returns the zero set of f as a sub vector space.
}

    if Degree(f) gt 1 then
        error "Polynomial has to be linear to define a subspace";
    end if;
    P := Parent(f);
    gens := [ P.i : i in [1..Ngens(P)]];
    equ := [ Zero(BaseRing(P)) : i in [1..Ngens(P)] ];
    for t in Terms(f) do
        mon := Monomials(t)[1];
        coeff := Coefficients(t)[1];
        equ[Position(gens, mon)] := coeff;
    end for;
    A := Matrix(BaseRing(P), Ngens(P), 1, equ);
    ker := Kernel(A);
    return ker;

end intrinsic;

//============================================================================
intrinsic NormalizeRationalHyperplaneEquation(f::RngMPolElt) -> RngMPolElt
{
	Produces a normalized (unique) equation for a hyperplane with rational coefficients.
}
    if Degree(f) gt 1 then
        error "Polynomial has to be linear.";
    end if;
    R := PolynomialRing(Rationals(), Rank(Parent(f)));
    AssignNames(~R, Names(Parent(f)));
    f := R!f;
    Coeffs := Coefficients(f);
    Mons := Monomials(f);
    Dens := [ Integers()!Denominator(Coeffs[i]) : i in [1..#Coeffs] ];
    LCM := LeastCommonMultiple(Dens)*Sign(Integers()!LeadingCoefficient(f));
    for i:=1 to #Coeffs do
        Coeffs[i] := Coeffs[i]*LCM;
    end for;

    return Polynomial(Coeffs,Mons);

end intrinsic;

//============================================================================
intrinsic RationalHyperplaneTuple(f::RngMPolElt) -> Tup
{The tuple of coefficients returned by NormalizeRationalHyperplaneEquation.}

    fnor := NormalizeRationalHyperplaneEquation(f);
    P := Parent(f);
    tup := <Zero(P) : i in [1..Ngens(P)]>;

    coeffs := Coefficients(fnor);
    mons := Monomials(fnor);
    gens := [ P.i : i in [1..Ngens(P)]];

    for i:=1 to #mons do
        j := Position(gens, mons[i]);
        if j ne 0 then
            tup[j] := coeffs[i];
        end if;
    end for;

    return tup;

end intrinsic;
