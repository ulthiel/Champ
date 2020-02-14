/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
    Simple extensions for structure constant algebras
*/

//==============================================================================
intrinsic StructureConstants(~A::AlgAss)
/*
    Intrinsic: StructureConstants

    Attaches the set of structure constants

    Declaration:
        :intrinsic StructureConstants(~A::AlgAss)
        :intrinsic StructureConstants(A::AlgAss) -> SetEnum

    Parameters:
       A - An algebra of type +AlgAss+

    Description:
        Assigns to the corresponding attribute the set of structure constants. This can be used for example to search for denominators in the structure constants.

*/
{}

    if assigned A`StructureConstants then
        return;
    end if;

    A`StructureConstants := { X[4] : X in BasisProducts(A : Rep:="Sparse") };

end intrinsic;

//==============================================================================
intrinsic StructureConstants(A::AlgAss) -> SetEnum
{}

    StructureConstants(~A);
    return A`StructureConstants;

end intrinsic;

//==============================================================================
/*
    Intrinsic: VectorSpace

    The underlying vector space

    Declaration:
        :intrinsic VectorSpace(~A::AlgAss)

    Parameters:
       A - An algebra of type +AlgAss+

    Description:
        Assigns to the corresponding attribute the underlying vector space of +A+. Also sets <AlgAss.VectorSpaceMap>.
*/
intrinsic VectorSpace(~A::AlgAss)
{}

    if assigned A`StructureConstants then
        return;
    end if;

    A`VectorSpace, A`VectorSpaceMap := VectorSpace(A);

end intrinsic;







//==============================================================================
//=============================================================================
/*
    Namespace: AlgAss

    Additions to the category +AlgAss+.
*/
declare attributes AlgAss:
    StructureConstants,
    /*
    	Attribute: StructureConstants

   	 	The structure constants
	*/
    VectorSpace,
    /*
    	Attribute: VectorSpace

    	The underlying vector space
	*/
    VectorSpaceMap;
    /*
    	Attribute: VectorSpaceMap

    	Map into the underlying vector space
	*/
