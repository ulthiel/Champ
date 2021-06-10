/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Simple extensions for modules over algebras (type ModRng)
*/

//============================================================================
intrinsic Specialize(M::ModRng, P::RngOrdIdl) -> ModRng
{
	Let M be a finite-dimensional module for a finite-dimensional algebra A
    over a number field K. If P is a prime ideal of the ring of integers of K,
    return the reduction M/PM, which is an (A/PA)-module.
}

    k,q := ResidueClassFieldGF(P);
    specmats := [ Specialize(ActionGenerator(M,i),P) : i in [1..#ActionGenerators(M)]];
    return RModule(specmats);

end intrinsic;

//============================================================================
intrinsic Specialize(M::ModRng, p::RngIntElt) -> ModRng
{
	Let M be a finite-dimensional module for a finite-dimensional algebra A
    over a number field K. If P is a prime ideal of the ring of integers of K,
    return the reduction M/PM, which is an (A/PA)-module.
}


    K := BaseRing(M);
    if Type(K) eq FldRat then
        L, quotmor := ResidueClassFieldGF(p);
    else
        P := RandomPrimeIdeal(K,p);
        L, quotmor := ResidueClassFieldGF(P);
    end if;
    return ChangeRing(M, L);

end intrinsic;
