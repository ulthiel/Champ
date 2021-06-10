/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Intrinsics for rewrite algebras

    History:
        Sunday, October 27, 2013 00:02:07: Initial.
*/

/*
	Minimal Magma Version: [2,19,0]
*/

declare type AlgRew[AlgRewElt];

declare attributes AlgRew:
    BaseRing,
    CoveringFreeAlgebra,
    Generators,
    AlwaysNormalize,
    RewriteRules;

declare attributes AlgRewElt:
    Parent,
    CoveringFreeAlgebraElement,
    IsInNormalForm;


//==============================================================================
intrinsic RewriteAlgebra(R::Rng, N::RngIntElt : AlwaysNormalize:=true) -> AlgRew
/*
    History:
        Sunday, October 27, 2013 00:03:54: Initial.
*/
{}
    H := New(AlgRew);

    H`CoveringFreeAlgebra := FreeAlgebra(R, N);

    H`Generators := [**];

    for i:=1 to N do
        x := New(AlgRewElt);
        x`Parent := H;
        x`CoveringFreeAlgebraElement := H`CoveringFreeAlgebra.i;
        x`IsInNormalForm := true;
        Append(~H`Generators,x);
    end for;

    H`BaseRing := R;

    H`RewriteRules := AssociativeArray({[1,2,3]});

    H`AlwaysNormalize := AlwaysNormalize;

    return H;

end intrinsic;




//==============================================================================
intrinsic Print(H::AlgRew)
/*
    History:
        Friday, September 20, 2013 16:27:50: Initial.
*/
{}

    printf "Rewrite algebra of rank %o over %o\n", Ngens(H), H`BaseRing;
    printf "Number of rewrite relations: %o\n", #H`RewriteRules;
    printf "Always normalize: %o", H`AlwaysNormalize;

end intrinsic;

//==============================================================================
intrinsic Print(x::AlgRewElt)
/*
    History:
        Friday, September 20, 2013 16:27:50: Initial.
*/
{}

    printf "%o", x`CoveringFreeAlgebraElement;

end intrinsic;

//==============================================================================
intrinsic AddRewriteRule(~H::AlgRew, code::SeqEnum, x::AlgRewElt)
/*
    History:
        Sunday, October 27, 2013 16:43:48: Initial.
*/
{}

    H`RewriteRules[code] := x;

end intrinsic;

//==============================================================================
intrinsic Ngens(H::AlgRew) -> RngIntElt
/*
    History:
        Friday, September 20, 2013 16:27:31: Initial.
*/
{}

    return Ngens(H`CoveringFreeAlgebra);

end intrinsic;

//==============================================================================
intrinsic Zero(H::AlgRew) -> AlgRewElt
/*
    History:
        Sunday, October 27, 2013 16:14:05: Initial.
*/
{}

    x := New(AlgRewElt);
    x`Parent := H;
    x`CoveringFreeAlgebraElement := Zero(H`CoveringFreeAlgebra);
    x`IsInNormalForm := true;

    return x;

end intrinsic;

//==============================================================================
intrinsic One(H::AlgRew) -> AlgRewElt
/*
    History:
        Sunday, October 27, 2013 16:14:05: Initial.
*/
{}

    x := New(AlgRewElt);
    x`Parent := H;
    x`CoveringFreeAlgebraElement := One(H`CoveringFreeAlgebra);
    x`IsInNormalForm := true;

    return x;

end intrinsic;

//==============================================================================
intrinsic Parent(x::AlgRewElt) -> AlgCheElt
/*
    History:
        Friday, September 20, 2013 16:29:07: Initial.
*/
{}

    return x`Parent;

end intrinsic;

//==============================================================================
intrinsic '.'(H::AlgRew, i::RngIntElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:30:08: Initial.
*/
{}

    return H`Generators[i];

end intrinsic;

//==============================================================================
intrinsic '*'(x::AlgRewElt, y::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement*y`CoveringFreeAlgebraElement;
    if x`Parent`AlwaysNormalize then
        return Normalize(z);
    else
        return z;
    end if;

end intrinsic;

//==============================================================================
intrinsic '*'(x::RngElt, y::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := y`Parent;
    z`CoveringFreeAlgebraElement := x*y`CoveringFreeAlgebraElement;
    if assigned y`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic '*'(x::AlgRewElt, y::RngElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement*y;
    if assigned x`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic '+'(x::AlgRewElt, y::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement+y`CoveringFreeAlgebraElement;
    if assigned x`IsInNormalForm and x`IsInNormalForm and assigned y`IsInNormalForm and y`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;
    return z;

end intrinsic;


//==============================================================================
intrinsic '+'(x::AlgRewElt, y::RngElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement+y;
    if assigned x`IsInNormalForm and x`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic '+'(y::RngElt, x::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement+y;
    if assigned x`IsInNormalForm and x`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(x::AlgRewElt, y::RngElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement-y;
    if assigned x`IsInNormalForm and x`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(y::RngElt, x::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement-y;
    if assigned x`IsInNormalForm and x`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(x::AlgRewElt, y::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement-y`CoveringFreeAlgebraElement;
    if assigned x`IsInNormalForm and x`IsInNormalForm and assigned y`IsInNormalForm and y`IsInNormalForm then
        z`IsInNormalForm := true;
    end if;

    return z;

end intrinsic;

//==============================================================================
intrinsic '-'(x::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:51:22: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := -x`CoveringFreeAlgebraElement;
    if assigned x`IsInNormalForm then
        z`IsInNormalForm := x`IsInNormalForm;
    end if;

    return z;

end intrinsic;

//==============================================================================
intrinsic '^'(x::AlgRewElt, i::RngIntElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:49:29: Initial.
*/
{}

    z := New(AlgRewElt);
    z`Parent := x`Parent;
    z`CoveringFreeAlgebraElement := x`CoveringFreeAlgebraElement^i;
    return z;

end intrinsic;

//==============================================================================
intrinsic Monomials(x::AlgRewElt) -> List
/*
    History:
        Sunday, October 27, 2013 00:38:22: Initial.
*/
{}

    mons := Monomials(x`CoveringFreeAlgebraElement);
    monscorrect := [* *];
    for i:=1 to #mons do
        y := New(AlgRewElt);
        y`Parent := x`Parent;
        y`CoveringFreeAlgebraElement := mons[i];
        if assigned x`IsInNormalForm then
            y`IsInNormalForm := x`IsInNormalForm;
        end if;
        Append(~monscorrect, y);
    end for;

    return monscorrect;

end intrinsic;

//==============================================================================
intrinsic Terms(x::AlgRewElt) -> List
/*
    History:
        Sunday, October 27, 2013 00:38:22: Initial.
*/
{}

    terms := Terms(x`CoveringFreeAlgebraElement);
    termscorrect := [* *];
    for i:=1 to #terms do
        y := New(AlgRewElt);
        y`Parent := x`Parent;
        y`CoveringFreeAlgebraElement := terms[i];
        if assigned x`IsInNormalForm then
            y`IsInNormalForm := x`IsInNormalForm;
        end if;
        Append(~termscorrect, y);
    end for;

    return termscorrect;

end intrinsic;

//==============================================================================
intrinsic Coefficients(x::AlgRewElt) -> List
/*
    History:
        Sunday, October 27, 2013 00:47:12: Initial
*/
{}

    return Coefficients(x`CoveringFreeAlgebraElement);

end intrinsic;

//==============================================================================
intrinsic ApplyRuleToTerm(x::AlgRewElt) -> AlgRewElt
/*
    History:
        Friday, September 20, 2013 16:53:13: Initial.
*/
{}
    if assigned x`IsInNormalForm and x`IsInNormalForm then
        return x;
    end if;

    mons := Monomials(x);

    if #mons gt 1 then
        error "Can apply rules only to terms.";
    end if;

    H := Parent(x);
    mon := mons[1];
    coeff := Coefficients(x)[1];
    code := Eltseq(mon`CoveringFreeAlgebraElement);
    normal := true;
    for key in Keys(H`RewriteRules) do
        if IsSubsequence(key, code) then
            normal := false;
            //find start of key in code
            p:=1;
            while p le #code-#key do
                if key eq [code[i] : i in [p..p+#key-1] ] then
                    break;
                end if;
                p +:= 1;
            end while;

            leftpart := ArrayProduct([ H`CoveringFreeAlgebra.code[i] : i in [1..p-1]] : OneElement:=One(H`CoveringFreeAlgebra));
            rightpart := ArrayProduct([ H`CoveringFreeAlgebra.code[i] : i in [p+#key..#code]] : OneElement:=One(H`CoveringFreeAlgebra));
            newelt := coeff*leftpart*H`RewriteRules[key]`CoveringFreeAlgebraElement*rightpart;

            y := New(AlgRewElt);
            y`Parent := H;
            y`CoveringFreeAlgebraElement := newelt;
            y`IsInNormalForm := false;

            break;

        end if;
    end for;

    if normal then
        x`IsInNormalForm := true;
        return x;
    else
        return y;
    end if;

end intrinsic;


//==============================================================================
intrinsic Normalize(x::AlgRewElt) -> AlgRewElt
/*
    History:
        Sunday, October 27, 2013 01:43:36: Initial.
*/
{}

    badterms := Terms(x);
    y := Zero(x`Parent);

    while not IsEmpty(badterms) do
        t := badterms[1];
        Remove(~badterms, 1);
        trewrite := ApplyRuleToTerm(t);
        if assigned t`IsInNormalForm and t`IsInNormalForm then
            y +:= trewrite;
        else
            badterms cat:=Terms(trewrite);
        end if;
    end while;

    return y;

end intrinsic;
