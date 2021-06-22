/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Simple extensions for free algebras.
*/



//==============================================================================
intrinsic Generators(~A::AlgFr)
/*
    Intrinsic: Generators

    Attaches the sequence of generators.

    Declaration:
        :intrinsic Generators(~A::AlgFr)
*/
{}

    if not assigned A`Generators then
        A`Generators := [ A.i : i in [1..Ngens(A)] ];
    end if;

end intrinsic;

//==============================================================================
intrinsic One(~A::AlgFr)
/*
    Intrinsic: One

    Attaches the unit element.

    Declaration:
        :intrinsic One(~A::AlgFr)
*/
{}

    if not assigned A`One then
        A`One := One(A);
    end if;

end intrinsic;

//==============================================================================
intrinsic Zero(~A::AlgFr)
/*
    Intrinsic: Zero

    Attaches the zero element.

    Declaration:
        :intrinsic Zero(~A::AlgFr)
*/
{}

    if not assigned A`Zero then
        A`Zero := Zero(A);
    end if;

end intrinsic;



//==============================================================================
/*
    Namespace: AlgFr

    Additions to the category +AlgFr+.
*/
declare attributes AlgFr:
    Generators,
/*
	Attribute: Generators

	Contains the sequence of generators.
*/
	One,
/*
	Attribute: One

	Contains the unit element.
*/
	Zero;
/*
	Attribute: Zero

	Contains the zero element.
*/
