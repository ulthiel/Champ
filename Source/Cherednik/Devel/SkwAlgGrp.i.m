/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

declare type SkwAlgGrp[SkwAlgGrpElt];

declare attributes SkwAlgGrp:
    BaseAlgebra,
    BaseGroup,
    BaseRing,
    CoveringFreeAlgebra,
    FromCoveringFreeAlgebraToBaseAlgebra,
    FromBaseAlgebraToSkewGroupAlgebra,
    Generators,
    AlwaysNormalize;

declare attributes SkwAlgGrpElt:
    Parent,
    ElementInCoveringFreeAlgebra;

//==============================================================================
intrinsic SkewGroupAlgebra(G::Grp, A::Alg : AlwaysNormalize:=true, UseCommutative:=false) -> SkwAlgGrp
/*
    History:
        Sunday, July 07, 2013 11:30:48: Initial.
*/
{Returns the Skw group algebra G#A. For this the action ^ of G on A has to be defined.}

    R := BaseRing(A);
    StrongGenerators(~G);

    if IsCommutative(A) or UseCommutative then
        F := PolynomialRing(R, Ngens(A) + G`NumberOfStrongGenerators);
    else
        F := FreeAlgebra(R, Ngens(A) + G`NumberOfStrongGenerators);
    end if;

    S := New(SkwAlgGrp);

    S`AlwaysNormalize := AlwaysNormalize;
    S`BaseAlgebra := A;
    S`BaseGroup := G;
    S`BaseRing := R;
    S`CoveringFreeAlgebra := F;

    S`ToBaseAlgebra := hom<S`FromCoveringFreeAlgebraToBaseAlgebra->S`BaseAlgebra | [ Zero(S`BaseAlgebra) : i in [1..G`NumberOfStrongGenerators]] cat [ S`BaseAlgebra.i : i in [1..Ngens(A)] ]>;



    S`Generators := [* *];

    for i:=1 to Ngens(F) do
        x := New(SkwAlgGrpElt);
        x`Parent := S;
        x`ElementInCoveringFreeAlgebra := F.i;
        Append(~S`Generators, x);
    end for;

    return S;

end intrinsic;

//==============================================================================
intrinsic AssignNames(~S::SkwAlgGrp, Names::SeqEnum[MonStgElt])
/*
    History:
        Sunday, July 07, 2013 12:18:35: Initial.
*/
{}

    AssignNames(~S`CoveringFreeAlgebra, Names);

end intrinsic;

//==============================================================================
intrinsic Ngens(S::SkwAlgGrp) -> RngIntElt
/*
    History:
        Monday, July 08, 2013 11:25:51: Initial.
*/
{}

    return Ngens(S`CoveringFreeAlgebra);

end intrinsic;


//==============================================================================
intrinsic Print(S::SkwAlgGrp)
/*
    History:
        Sunday, July 07, 2013 11:40:00: Initial.
*/
{}

    printf "Skew group algebra over %o\n", S`BaseRing;
    printf "Underlying group:\n";
    IndentPush();
    printf "%o\n", S`BaseGroup;
    IndentPop();
    printf "Underlying algebra:\n";
    IndentPush();
    printf "%o\n", S`BaseAlgebra;
    IndentPop();
    printf "Variables: %o", Names(S`CoveringFreeAlgebra);
    IndentPush();
    //printf "Strong generators for group: %o",

end intrinsic;

//==============================================================================
intrinsic Print(x::SkwAlgGrpElt)
/*
    History:
        Sunday, July 07, 2013 11:56:22: Initial.
*/
{}

    printf "%o", x`ElementInCoveringFreeAlgebra;

end intrinsic;

//==============================================================================
intrinsic Parent(x::SkwAlgGrpElt) -> SkwAlgGrp
/*
    History:
        Sunday, July 07, 2013 12:00:05: Initial.
*/
{}

    return x`Parent;

end intrinsic;

//==============================================================================
intrinsic '.'(S::SkwAlgGrp, i::RngIntElt) -> SkwAlgGrpElt
/*
    History:
        Sunday, July 07, 2013 20:42:39: Initial.
*/
{}

    return S`Generators[i];

end intrinsic;

//==============================================================================
intrinsic '*'(x::SkwAlgGrpElt, y::SkwAlgGrpElt) -> SkwAlgGrpElt
/*
    History:
        Sunday, July 07, 2013 12:11:12: Initial.
*/
{}
    z := New(SkwAlgGrpElt);
    z`ElementInCoveringFreeAlgebra := x`ElementInCoveringFreeAlgebra*y`ElementInCoveringFreeAlgebra;
    z`Parent := x`Parent;
    return z;

end intrinsic;

//==============================================================================
intrinsic '+'(x::SkwAlgGrpElt, y::SkwAlgGrpElt) -> SkwAlgGrpElt
/*
    History:
        Sunday, July 07, 2013 12:11:12: Initial.
*/
{}
    z := New(SkwAlgGrpElt);
    z`ElementInCoveringFreeAlgebra := x`ElementInCoveringFreeAlgebra+y`ElementInCoveringFreeAlgebra;
    z`Parent := x`Parent;
    return z;

end intrinsic;

//==============================================================================
intrinsic '^'(x::SkwAlgGrpElt, n::RngIntElt) -> SkwAlgGrpElt
/*
    History:
        Sunday, July 07, 2013 12:11:12: Initial.
*/
{}
    z := New(SkwAlgGrpElt);
    z`ElementInCoveringFreeAlgebra := x`ElementInCoveringFreeAlgebra^n;
    z`Parent := x`Parent;
    return z;

end intrinsic;


//==============================================================================
intrinsic Monomials(x::SkwAlgGrpElt) -> SeqEnum
/*
    History:
        Monday, July 08, 2013 11:30:25: Initial.
*/
{}

    return Monomials

//==============================================================================
//intrinsic NormalForm(x::SkwAlgGrpElt) -> SkwAlgGrpElt
/*
    History:
        Monday, July 08, 2013 11:29:39: Initial.
*/
//{}



//end intrinsic;
