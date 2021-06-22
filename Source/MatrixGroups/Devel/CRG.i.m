/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    Intrinsics around complex reflection groups.
*/

//==============================================================================
intrinsic ImprimitiveCRGOrderCandidates(o::RngIntElt : Verbose:=true) -> SetEnum
/*
    History:
        Monday, August 12, 2013 15:06:25: Initial.
*/
{All triples (m,p,n) such that (m^n*n!)/p = o.}

    triples := {@ <o, 1, 1> @};

    nbound := 2;
    while Factorial(nbound) le o do
        nbound +:= 1;
    end while;

    for n:=2 to nbound do
        mupperbound := Ceiling( (o / Factorial(n))^(1/(n-1)));
        mlowerbound := Floor( (o / Factorial(n))^(1/(n)));
        print n,mlowerbound, mupperbound;
        for m:=mlowerbound to mupperbound do
            if m^(n-1)*Factorial(n) gt o then
                continue;
            end if;
            if m^n*Factorial(n) lt o then
                continue;
            end if;
            if Verbose and m mod 10000 eq 0 then
                PrintPercentage(m, mupperbound);
            end if;
            for p in Divisors(m) do
                if (m^n*Factorial(n))/p eq o then
                    triples join:={@<m,p,n>@};
                end if;
            end for;
        end for;
    end for;
    return triples;

end intrinsic;

//==============================================================================
intrinsic CheckIsomorphismBetweenExceptionalAndInfiniteSeries(i::RngIntElt)
/*
    History:
        Monday, August 12, 2013 21:32:42: Initial.
*/
{Checks which G(m,p,n) are isomorphic to the exceptional group G_i (hopefully none!).}

    G := ShephardTodd(i);
    cands := ImprimitiveCRGOrderCandidates(#G);
    for j:=1 to #cands do
        T := cands[j];
        m := T[1];
        p := T[2];
        n := T[3];
        //if IsIsomorphic(G,ShephardTodd(T[1],T[2],T[3])) then
        //H := ShephardTodd(T[1],T[2],T[3]);
        //if PrimaryInvariants(Center(G)) eq PrimaryInvariants(Center(H)) and PrimaryInvariants(G/CommutatorSubgroup(G)) eq PrimaryInvariants(H/CommutatorSubgroup(H)) then
        Zdeg := ImprimitiveReflectionGroupCenterOrder(m,p,n);
        if Zdeg eq #Center(G) then
            print "Could be isomorphic to "*Sprint(T);
        end if;
        PrintPercentage(j, #cands);
    end for;

    print "Everything OK.";

end intrinsic;
