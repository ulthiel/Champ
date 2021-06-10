/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions for matrix algebras.
*/


//==============================================================================
intrinsic RModule(~A::AlgMat)
{}

    if assigned A`RModule then
        return;
    end if;

    A`RModule := RModule(A); 	//REPAIR!

end intrinsic;

//==============================================================================
intrinsic LeftNaturalModule(~A::AlgMat)
{}

    StructureConstantAlgebra(~A);


end intrinsic;

//==============================================================================
intrinsic SimpleModules(~A::AlgMat)
{}

    if assigned A`SimpleModules then
        return;
    end if;

    RModule(~A);
    A`SimpleModules := Constituents(A`RModule);	//this is WRONG! Need the natural module for this!

end intrinsic;

//==============================================================================
intrinsic SimpleModules(A::AlgMat) -> SeqEnum
{}

    SimpleModules(~A);
    return A`SimpleModules;

end intrinsic;

//==============================================================================
intrinsic JacobsonRadical(~A::AlgMat)
{}

    if assigned A`JacobsonRadical then
        return;
    end if;

    A`JacobsonRadical := JacobsonRadical(A);

end intrinsic;

//==============================================================================
intrinsic IsSplit(A::AlgMat) -> BoolElt
{}

    JacobsonRadical(~A);
    SimpleModules(~A);
    return Dimension(A) eq Dimension(A`JacobsonRadical) + &+[Dimension(S)^2 : S in A`SimpleModules ];

end intrinsic;


/*
    Namespace: AlgMat

    Additions to the category +AlgMat+.
*/
declare attributes AlgMat:
    SimpleModules,
    /*
    	Attribute: SimpleModules

    	The simple modules
    */
    JacobsonRadical,
    /*
    	Attribute: JacobsonRadical

    	The Jacobson radical
    */
    RModule,
    /*
    	Attribute: RModule

    	The natural R-module
    */
    StructureConstantAlgebra;
    /*
    	Attribute: StructureConstantAlgebra

    	The structure constant algebra.
    */
