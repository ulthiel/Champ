/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions for invariant rings.

    *CAUTION:* Although in the Magma documentation the notation K[V] is used, this is (as explained in the documentation) equal to S(V). As we want understand the invariant theory of a matrix group G to be the invariant theory of K[V] = S(V^*) we have to dualize everything. So, we provide here some wrappers so that the corresponding invariant theoretic attributes for matrix group are relative to K[V] = S(V^*).

    History:
        * Switching x and y notation again, to do the same as in BR.
        * Sunday, August 25, 2013 16:12:23: Initial
*/


//===========================================================================
declare attributes GrpMat:
    CoordinateAlgebra, //The coordinate algebra K[V] = S(V^*)
    PrimaryInvariants,
    SecondaryInvariants,
    FundamentalInvariants, //Fundamental invariants (non-unique)
    HilbertIdeal,
    Degrees, //Degrees of the fundamental invariants (unique)
    InvariantRing, //The invariant ring (of the dual group!)
    IsPolynomial,
    IsReflectionGroup;

//now, everything again for the dual representation
declare attributes GrpMat:
    DualCoordinateAlgebra,
    DualPrimaryInvariants,
    DualSecondaryInvariants,
    DualFundamentalInvariants,
    DualHilbertIdeal,
    DualDegrees,
    DualInvariantRing,
    IsDualPolynomial;

//===========================================================================
intrinsic InvariantRing(~G::GrpMat)
{The invariant theory of G on K[V] = S(V^*). Note the dual!}

    if assigned G`InvariantRing then
        return;
    end if;

    DualGroup(~G); //this attaches the correct invariant ring already

end intrinsic;

//===========================================================================
intrinsic DualGroup(~G::GrpMat)
{The dual group of a matrix group. This also attaches the correct invariant ring to G`DualGroup.}

    if assigned G`InvariantRing then //yes, only in this way we can see that G is not the dual group of something in which case we would mess up x and y
        return;
    end if;

    G`DualGroup := MatrixGroup<Dimension(G), BaseRing(G) | [ Transpose(G.i^-1) : i in [1..Ngens(G)]]>;

    //we have to set correct stuff already here.
    G`InvariantRing := InvariantRing(G`DualGroup);
    G`DualGroup`InvariantRing := InvariantRing(G);
    G`CoordinateAlgebra := PolynomialRing(G`InvariantRing);
    Generators(~G`CoordinateAlgebra);
    G`CoordinateAlgebra`Zero := Zero(G`CoordinateAlgebra);
    G`CoordinateAlgebra`One := One(G`CoordinateAlgebra);
    G`DualGroup`CoordinateAlgebra := PolynomialRing(G`DualGroup`InvariantRing);
    Generators(~G`DualGroup`CoordinateAlgebra);
    AssignNames(~G`CoordinateAlgebra, ["x"*Sprint(i) : i in [1..Dimension(G)]]);
    AssignNames(~G`DualGroup`CoordinateAlgebra, ["y"*Sprint(i) : i in [1..Dimension(G)]]);

    //now, everything for the dual
    G`DualInvariantRing := InvariantRing(G);
    G`DualCoordinateAlgebra := G`DualGroup`CoordinateAlgebra;
    Generators(~G`DualCoordinateAlgebra);
    G`DualCoordinateAlgebra`Zero := Zero(G`DualCoordinateAlgebra);
    G`DualCoordinateAlgebra`One := One(G`DualCoordinateAlgebra);

    if assigned G`IsReflectionGroup then
    	G`DualGroup`IsReflectionGroup := G`IsReflectionGroup;
    end if;

end intrinsic;

//===========================================================================
intrinsic CoordinateAlgebra(~G::GrpMat)
{
	Attaches the coordinate algebra K[V] = S(V^*).
}

	DualGroup(~G);

	G`Dimension := Dimension(G);
	G`CoordinateAlgebra`Zero := Zero(G`CoordinateAlgebra);
	G`CoordinateAlgebra`One := One(G`CoordinateAlgebra);
	G`CoordinateAlgebra`Generators := [G`CoordinateAlgebra.i : i in [1..G`Dimension] ];

end intrinsic;

//===========================================================================
intrinsic CoordinateAlgebra(G::GrpMat) -> RngMPolElt
{
	The coordinate algebra K[V] = S(V^*).
}

	CoordinateAlgebra(~G);
	return G`CoordinateAlgebra;

end intrinsic;

//===========================================================================
intrinsic DualCoordinateAlgebra(G::GrpMat) -> RngMPolElt
{}

	DualGroup(~G);
	return G`DualCoordinateAlgebra;

end intrinsic;


//===========================================================================
intrinsic FundamentalInvariants(~G::GrpMat)
{Attaches the fundamental invariants of a matrix group of G.}

    if assigned G`FundamentalInvariants then
        return;
    end if;

    DualGroup(~G);
    G`FundamentalInvariants := FundamentalInvariants(G`InvariantRing);

end intrinsic;

//===========================================================================
intrinsic FundamentalInvariants(G::GrpMat) -> SeqEnum
{}

    FundamentalInvariants(~G);
    return G`FundamentalInvariants;

end intrinsic;

//===========================================================================
intrinsic DualFundamentalInvariants(~G::GrpMat)
{Attaches the dual fundamental invariants of a matrix group of G.}

    if assigned G`DualFundamentalInvariants then
        return;
    end if;

	DualGroup(~G);
    G`DualFundamentalInvariants := FundamentalInvariants(G`DualInvariantRing);

end intrinsic;

//===========================================================================
intrinsic DualFundamentalInvariants(G::GrpMat) -> SeqEnum
{}

    DualFundamentalInvariants(~G);
    return G`DualFundamentalInvariants;

end intrinsic;

//===========================================================================
intrinsic FundamentalInvariants(g::GrpMatElt) -> SeqEnum
{}

    H := MatrixGroup<Ncols(g),BaseRing(g)|[g]>;
    return FundamentalInvariants(H);

end intrinsic;

//===========================================================================
intrinsic PrimaryInvariants(~G::GrpMat)
{}

    if not assigned G`PrimaryInvariants then
        DualGroup(~G);
        G`PrimaryInvariants := PrimaryInvariants(G`InvariantRing);
    end if;

end intrinsic;

//===========================================================================
intrinsic DualPrimaryInvariants(~G::GrpMat)
{}

	if not assigned G`DualPrimaryInvariants then
        DualGroup(~G);
        G`DualPrimaryInvariants := PrimaryInvariants(G`DualInvariantRing);
    end if;

end intrinsic;

//===========================================================================
intrinsic HilbertIdeal(~G::GrpMat)
{Attaches the Hilbert ideal of a matrix group of G.}

    if assigned G`HilbertIdeal then
        return;
    end if;

    DualGroup(~G);

    FundamentalInvariants(~G);

    G`HilbertIdeal := ideal<G`CoordinateAlgebra|G`FundamentalInvariants>;

end intrinsic;

//===========================================================================
intrinsic HilbertIdeal(G::GrpMat) -> SeqEnum
{The Hilbert ideal of G.}

	HilbertIdeal(~G);
	return G`HilbertIdeal;

end intrinsic;

//===========================================================================
intrinsic DualHilbertIdeal(~G::GrpMat)
{Attaches the dual Hilbert ideal of a matrix group of G.}

	if assigned G`DualHilbertIdeal then
        return;
    end if;

    DualGroup(~G);

    DualFundamentalInvariants(~G);

    G`DualHilbertIdeal := ideal<G`DualCoordinateAlgebra|G`DualFundamentalInvariants>;

end intrinsic;

//===========================================================================
intrinsic DualHilbertIdeal(G::GrpMat) -> SeqEnum
{The dual Hilbert ideal of a matrix group of G.}

	DualHilbertIdeal(~G);
	return G`DualHilbertIdeal;

end intrinsic;

//===========================================================================
intrinsic Degrees(~G::GrpMat : UseDB:=true)
{
	Attaches the degrees of a matrix group of G. If UseDB is set to true, then these will be loaded from the data base if available.
}

    if assigned G`Degrees then
        return;
    end if;

    DualGroup(~G);
    foundindb := false;
    if UseDB and assigned G`DBDir and CHAMP_ExistsInDB(G`DBDir, "Invariants/Degrees") then
        G`Degrees := CHAMP_GetFromDB(G`DBDir, "Invariants/Degrees");
        return;
    end if;

    FundamentalInvariants(~G);

    G`Degrees := [ Degree(f) : f in G`FundamentalInvariants ];

	//for non-modular reflection groups the dual degrees are the same as the degrees
	if assigned G`IsReflectionGroup and G`IsReflectionGroup and IsNonModular(G) then
		G`DualDegrees := G`Degrees;
	end if;

end intrinsic;

//===========================================================================
intrinsic Degrees(G::GrpMat) -> SeqEnum
{}

    Degrees(~G);
    return G`Degrees;

end intrinsic;

//===========================================================================
intrinsic IsPolynomial(~G::GrpMat)
{True iff the invariant ring of G is polynomial.}

    if assigned G`IsPolynomial then
        return;
    end if;

    IsNonModular(~G);
    IsReflectionGroup(~G);

    if G`IsNonModular and G`IsReflectionGroup then
        G`IsPolynomial := true; //nonmodular reflection groups have polynomial invariants

    elif assigned G`FundamentalInvariants then

        //see [Thi14] for this condition
        if #G`FundamentalInvariants eq Dimension(G) then
            G`IsPolynomial := true;
        else
            G`IsPolynomial := false;
        end if;

    else
        FundamentalInvariants(~G);
        IsPolynomial(~G);
    end if;

end intrinsic;

//===========================================================================
intrinsic IsPolynomial(G::GrpMat) -> BoolElt
{True iff G has polynomial invariant ring.}

    IsPolynomial(~G);
    return G`IsPolynomial;

end intrinsic;


//============================================================================
intrinsic FundamentalInvariantsSingularCode(G::GrpMat) -> MonStgElt
{
	Returns the Singular code to compute the fundamental invariants of G.
}
	if Type(BaseRing(G)) ne FldRat then
		error "Only for rational field as base ring.";
	end if;
	CoordinateAlgebra(~G);
	str := "LIB \"finvar.lib\";\n";
	str *:= "ring R=0,(";
	for i:=1 to #Names(G`CoordinateAlgebra) do
		str *:= Names(G`CoordinateAlgebra)[i];
		if i lt #Names(G`CoordinateAlgebra) then
			str *:= ",";
		end if;
	end for;
	str *:= "),dp;\n";
	n := Dimension(G);
	for i:=1 to Ngens(G) do
		str *:= "matrix G"*Sprint(i)*"["*Sprint(n)*"]["*Sprint(n)*"]="*Replace(Replace(Sprint(Eltseq(Transpose(G.i))),"\\[",""),"\\]", "")*";\n";
	end for;
	str *:= "int t=timer;\n";
	str *:= "list reyn=group_reynolds(";
	for i:=1 to Ngens(G) do
		str *:= "G"*Sprint(i);
		if i lt Ngens(G) then
			str *:=",";
		end if;
	end for;
	str *:= ");\n";
	str *:= "matrix inv = invariant_algebra_reynolds(reyn[1],1);\n";
	str *:= "timer-t;\n";

	return str;

end intrinsic;
