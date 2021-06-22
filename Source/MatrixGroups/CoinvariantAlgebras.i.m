/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Intrinsics for coinvariant algebras of matrix groups.
*/


//===========================================================================
declare attributes GrpMat:
	CoinvariantAlgebraPoincareSeries,
	CoinvariantAlgebra,
    ToCoinvariantAlgebra,
    FromCoinvariantAlgebra,
    CoinvariantAlgebraGradedGModule, //The graded G-module defined by the action of +G+ on its coinvariant algebra.
    CoinvariantAlgebraGradedModule; //The graded module defined by the action of the standard generators +x_i+ on the coinvariant algebra.

//now, same for the dual
declare attributes GrpMat:
	DualCoinvariantAlgebra,
    ToDualCoinvariantAlgebra,
    FromDualCoinvariantAlgebra;

//============================================================================
intrinsic CoinvariantAlgebra(~G::GrpMat : Verbose:=false)
{
	Attaches the coinvariant algebra of G.
}

    if assigned G`CoinvariantAlgebra then
        return;
    end if;

    if Verbose then
    	print "Computing coinvariant algebra.";
    end if;

    HilbertIdeal(~G);
    Groebner(G`HilbertIdeal);

    G`CoinvariantAlgebra, G`ToCoinvariantAlgebra := quo<G`CoordinateAlgebra|Basis(G`HilbertIdeal)>;

    G`FromCoinvariantAlgebra := hom<G`CoinvariantAlgebra->G`CoordinateAlgebra|[G`CoordinateAlgebra.i : i in [1..Dimension(G)]]>;

    G`CoinvariantAlgebra`Zero := Zero(G`CoinvariantAlgebra);
    G`CoinvariantAlgebra`One := One(G`CoinvariantAlgebra);
    G`CoinvariantAlgebra`Generators := [G`CoinvariantAlgebra.i : i in [1..Ngens(G`CoinvariantAlgebra)]];

end intrinsic;



//============================================================================
intrinsic CoinvariantAlgebra(G::GrpMat : Verbose:=false) -> RngMPolRes, Map
{
	Returns the coinvariant algebra of G.
}

    CoinvariantAlgebra(~G : Verbose:=Verbose);
    return G`CoinvariantAlgebra, G`ToCoinvariantAlgebra;

end intrinsic;

//============================================================================
intrinsic DualCoinvariantAlgebra(~G::GrpMat : Verbose:=false)
{
	Attaches the dual coinvariant algebra of G.
}

    if assigned G`DualCoinvariantAlgebra then
        return;
    end if;

    if Verbose then
    	print "Computing dual coinvariant algebra.";
    end if;

	DualHilbertIdeal(~G);
	Groebner(G`DualHilbertIdeal);

	G`DualCoinvariantAlgebra, G`ToDualCoinvariantAlgebra := quo<G`DualCoordinateAlgebra|Basis(G`DualHilbertIdeal)>;

    G`FromDualCoinvariantAlgebra := hom<G`DualCoinvariantAlgebra->G`DualCoordinateAlgebra|[G`DualCoordinateAlgebra.i : i in [1..Dimension(G)]]>;

    G`DualCoinvariantAlgebra`Zero := Zero(G`DualCoinvariantAlgebra);
    G`DualCoinvariantAlgebra`One := One(G`DualCoinvariantAlgebra);
    G`DualCoinvariantAlgebra`Generators := [G`DualCoinvariantAlgebra.i : i in [1..Ngens(G`DualCoinvariantAlgebra)]];

end intrinsic;

//============================================================================
intrinsic DualCoinvariantAlgebra(G::GrpMat : Verbose:=false) -> RngMPolRes
{
	The dual coinvariant algebra of G.
}

	DualCoinvariantAlgebra(~G : Verbose:=Verbose);
	return G`DualCoinvariantAlgebra, G`ToDualCoinvariantAlgebra;

end intrinsic;

//============================================================================
intrinsic CoinvariantAlgebraGradedGModule(~G::GrpMat : Rep:="Sparse", Verbose:=false)
{
	The graded G-module defined by the action of G on its coinvariant algebra.
}

    if assigned G`CoinvariantAlgebraGradedGModule then
        return;
    end if;

    CoinvariantAlgebra(~G : Verbose:=Verbose);

    if Verbose then
    	print "Computing graded G-module of coinvariant algebra.";
    end if;

    G`CoinvariantAlgebraGradedGModule := GradedGModule(G`CoinvariantAlgebra, G :  Rep:=Rep, Verbose:=Verbose, Dual:=true);

end intrinsic;

//============================================================================
intrinsic CoinvariantAlgebraGradedGModule(G::GrpMat : Rep:="Sparse", Verbose:=false) -> RngMPolRes
{}

	CoinvariantAlgebraGradedGModule(~G : Rep:=Rep, Verbose:=Verbose);
	return G`CoinvariantAlgebraGradedGModule;

end intrinsic;

//============================================================================
intrinsic CoinvariantAlgebraGradedModule(~G::GrpMat : Rep:="Sparse", Verbose:=false)
{}

	 if assigned G`CoinvariantAlgebraGradedModule then
        return;
    end if;

    CoinvariantAlgebra(~G);

    if Verbose then
    	print "Computing CoinvariantAlgebraGradedModule.";
    end if;

    G`CoinvariantAlgebraGradedModule := GradedModule(G`CoinvariantAlgebra : Rep:=Rep, Verbose:=Verbose);

end intrinsic;

//============================================================================
intrinsic CoinvariantAlgebraGradedModule(G::GrpMat : Rep:="Sparse", Verbose:=false) -> ModGr
{}

	CoinvariantAlgebraGradedModule(~G : Rep:=Rep, Verbose:=Verbose);
	return G`CoinvariantAlgebraGradedModule;

end intrinsic;

//============================================================================
intrinsic CoinvariantAlgebraPoincareSeries(~G::GrpMat)
{
	The Poincare series of the coinvariant algebra of G. If the invariant ring of G is polynomial, we use a formula, otherwise we compute it directly.
}

    if not assigned G`CoinvariantAlgebraPoincareSeries then
        P<t> := PolynomialRing(Integers());
        IsPolynomial(~G);
        if G`IsPolynomial then
            Degrees(~G);
            G`CoinvariantAlgebraPoincareSeries := &*[ &+[ t^j : j in [0..G`Degrees[i]-1]] : i in [1..#G`Degrees] ] ;
        else
            CoinvariantAlgebra(~G);
            Basis(~G`CoinvariantAlgebra);
            p := 0;
            for i:=1 to #G`CoinvariantAlgebra`Basis do
                p +:= t^Degree(G`CoinvariantAlgebra`Basis[i]);
            end for;
        end if;
    end if;

end intrinsic;
