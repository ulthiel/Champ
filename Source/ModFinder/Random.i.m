/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

//==============================================================================
intrinsic RandomPrime(S::SetEnum) -> RngIntElt
/*
    History:
        Monday, July 01, 2013 11:15:51: Initial.

        Sunday, September 22, 2013 15:29:45: Ported from GrpMat.i.m
*/
{A random prime number in S.}

    //this is very slow and should be improved
    P := SequenceToSet(PrimesInInterval(Minimum(S), Maximum(S))) meet S;
    if IsEmpty(P) then
        error "There is no prime number in this range.";
    end if;
    return Random(P);

end intrinsic;

//==============================================================================
intrinsic RandomPrimeIdeal(K::Fld : PreferTotallySplit:=true, pRange:={2..10^6}) -> RngOrdIdl, RngIntElt
/*
    History:
        Monday, July 01, 2013 11:13:14: Initial.

        Sunday, September 22, 2013 15:29:51: Ported from GrpMat.i.m
*/
{A random prime ideal in K lying over a (totally split) prime number p between in pRange.}

    if Type(K) eq FldRat then
        p := RandomPrime(pRange);
        return p,p;
    elif Type(K) eq FldNum or Type(K) eq FldCyc then
        O := RingOfIntegers(K);
        primesinrange := SequenceToSet(PrimesInInterval(Minimum(pRange), Maximum(pRange))) meet pRange;
        while true do
            p := Random(primesinrange);
            primesinrange diff:={p};
            fact := Factorization(ideal<O|p>);
            if PreferTotallySplit and not IsTotallySplit(fact[1][1]) and not IsEmpty(primesinrange) then
                continue;
            end if;
            primes := [ fact[i][1] : i in [1..#fact] ];
            P := Random(primes);
            return P,p;
        end while;
    else
        error "No method implemented.";
    end if;

end intrinsic;

//==============================================================================
intrinsic RandomPrimeIdeal(K::FldNum, p::RngIntElt) -> RngOrdIdl
/*
    History:
        Sunday, September 22, 2013 15:39:33: Initial.
*/
{A random prime ideal in K lying over p.}

    return RandomPrimeIdeal(K : pRange:={p});

end intrinsic;
