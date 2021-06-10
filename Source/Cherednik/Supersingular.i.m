/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Supersingular characters (see my counter-example paper)
*/

declare attributes GrpMat:
	FakeDegreeResidues,
	SupersingularRepresentations;

//============================================================================
intrinsic FakeDegreeResidues(~G::GrpMat)
{}
	if assigned G`FakeDegreeResidues then
        return;
    end if;

    FakeDegrees(~G);
    CoinvariantAlgebraPoincareSeries(~G);

    G`FakeDegreeResidues := [];
    P<t> := PolynomialRing(Integers():Global:=false);

    p := Characteristic(BaseRing(G));
    Modules(~G,p);

    for i:=1 to #G`Modules[p] do
        num := (Integers()!Dimension(G`Modules[p][i]))*t^Valuation(G`FakeDegrees[i])*G`CoinvariantAlgebraPoincareSeries;
        res := num mod G`FakeDegrees[i];
        Append(~G`FakeDegreeResidues, res);
    end for;

end intrinsic;

//============================================================================
intrinsic SupersingularRepresentations(~G::GrpMat)
{The list of supersingular representations of G.}

    if assigned G`SupersingularRepresentations then
        return;
    end if;

    FakeDegreeResidues(~G);
    p := Characteristic(BaseRing(G));

    G`SupersingularRepresentations := [i : i in [1..#G`Representations[p]] | G`FakeDegreeResidues[i] ne 0 ];

end intrinsic;

//============================================================================
intrinsic SupersingularRepresentations(G::GrpMat) -> SeqEnum
{The list of supersingular representations of G.}

	SupersingularRepresentations(~G);
	return G`SupersingularRepresentations;

end intrinsic;

//============================================================================
intrinsic IsSupersingular(f::RngUPolElt[RngInt], g::RngUPolElt[RngInt]) -> SetEnum
{True iff (f,g) is supersingular.}

    fact := Factorization(f);
    mults := [ fact[i][2] : i in [1..#fact] ];
    dom := CartesianProduct([ [0..mults[i]] : i in [1..#fact] ]);

    for d in dom do
        q := &*[ fact[i][1]^d[i] : i in [1..#fact] ];
        print q;
    end for;

    print fact;
    print dom;
    return {};

end intrinsic;
