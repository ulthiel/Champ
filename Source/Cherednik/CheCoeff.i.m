/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Intrinsics for computing the commutator [y_i,x^mu]_{s}
*/

declare attributes GrpMat:
	CherednikCoefficients; //table of the Coefficients (y_i,x_j)_{s_k}, indexed by [i][j][k]

//=========================================================================
intrinsic CanonicalPairing(x::ModTupFldElt, y::ModTupFldElt) -> FldElt
{If x and y are elements of a vector space V, we can interpret y as an element of V^* with respect to the dual basis of V and then compute <x,y> = y(x)}

    if Parent(x) ne Parent(y) then
        error "Both vectors have to lie in the same space";
    end if;
    n := Dimension(Parent(x));
    return &+[x[i]*y[i] : i in [1..n]];

end intrinsic;


//=========================================================================
intrinsic CherednikCoefficient(x::ModTupFldElt, y::ModTupFldElt, s::Rec) -> FldElt
/*
    (y,x)_s
*/
{Compute the Cherednik coefficient (x,y)_s for vectors x in V, y in V^* and a reflection record s.}

    if not s`IsReflection then
        error "Element has to be a reflection.";
    end if;

    if s`IsDiagonalizable eq false then
        error "Reflection has to be diagonalizable.";
    end if;

    return CanonicalPairing(x,s`Coroot)*CanonicalPairing(s`Root,y)/CanonicalPairing(s`Root,s`Coroot);

end intrinsic;

//=========================================================================
intrinsic CherednikCoefficient(i::RngIntElt, j::RngIntElt, s::Rec) -> FldElt
/*
    (y_i,x_j)_s
*/
{Compute the Cherednik coefficient (x,y)_s for vectors x in V, y in V^* and a reflection record s.}

    W := Parent(s`Element);
    try
    	k := s`ReflectionNumber; //may not be assigned
    	return W`CherednikCoefficients[i][j][k];
    catch e
    	V := VectorSpace(W);
    	return CherednikCoefficient(V.i,V.j,s);
	end try;

end intrinsic;

//=========================================================================
intrinsic CherednikCoefficient(x::ModTupFldElt, y::ModTupFldElt, s::GrpMatElt) -> FldElt
/*
    (y,x)_s
*/
{Compute the Cherednik coefficient (x,y)_s for vectors x in V, y in V^* and a reflection s.}

    return CherednikCoefficient(x,y,ReflectionData(s));

end intrinsic;


//=========================================================================
intrinsic CherednikCoefficient(i::RngIntElt, j::RngIntElt, s::GrpMatElt) -> FldElt
/*
    (y_i,x_j)_s
*/
{Compute the Cherednik coefficient (x,y)_s for vectors x in V, y in V^* and a reflection s.}

    return CherednikCoefficient(i,j,ReflectionData(s));

end intrinsic;

//============================================================================
intrinsic CherednikCoefficients(~W::GrpMat)
{
	Attaches the library of Cherednik coefficients.
}

	if assigned W`CherednikCoefficients	then
		return;
	end if;

	ReflectionLibrary(~W);
	N := #W`ReflectionLibraryFlat;
	d := W`Dimension;

	// just a safe check
	for k:=1 to N do
		assert W`ReflectionLibraryFlat[k]`ReflectionNumber eq k;
	end for;

	W`CherednikCoefficients := < < < CherednikCoefficient(i,j,W`ReflectionLibraryFlat[k]) : k in [1..N]> : j in [1..d]> : i in [1..d]>;

end intrinsic;
