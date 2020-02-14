/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
    Simple extensions for vector spaces.
*/

//============================================================================
intrinsic Quotient(V::ModTupFld, U::ModTupFld : Check:=false) -> ModTupFld, Map, Map
{Returns the quotient V/U, the quotient morphism and a section of the quotient morphism}

    C:=Complement(V,U);
    Q,q:=quo<V|U>;
    s := hom<Q->C|[C.i : i in [1..Dimension(Q)]]>;
    if Check then
        if exists{i : i in [1..Dimension(Q)] | q(s(Q.i)) ne Q.i } then
            error "Something wrong with the section morphism.";
        end if;
    end if;
    return Q,q,s;

end intrinsic;
