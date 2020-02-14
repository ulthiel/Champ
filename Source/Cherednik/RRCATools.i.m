/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
Intrinsics around restricted rational Cherednik algebras
*/


//============================================================================
intrinsic BadPrimesForRRCA(G::GrpMat) -> SetEnum
{The bad primes for the restricted rational Cherednik algebra of G as defined in [Thi14]. If bases of the coinvariant algebras of G are already set, these are used, otherwise they are computed automatically.}

    print "Computing coinvariant algebra.";
    A := CoinvariantAlgebra(G);
    print "Transforming into structure constant algebra.";
    AS := StructureConstantAlgebra(A);
    if not {@ A.i : i in [1..Ngens(A)] @} subset A`Basis then
        error "There is a generator not in the basis.";
    end if;
    print "Computing structure constants.";
    D1 := PrimeFactors({ Denominator(x) : x in StructureConstants(AS)});
    IndentPush();
    print D1;
    IndentPop();

    print "Computing dual coinvariant algebra.";
    DualGroup(~G);
    A := CoinvariantAlgebra(G`DualGroup);
    print "Transforming into structure constant algebra.";
    AS := StructureConstantAlgebra(A);
    if not {@ A.i : i in [1..Ngens(A)] @} subset A`Basis then
        error "There is a generator not in the basis.";
    end if;
    print "Computing structure constants.";
    D2 := PrimeFactors({ Denominator(x) : x in StructureConstants(AS)});
    IndentPush();
    print D2;
    IndentPop();

    Gentries := SequenceToSet(FlatFixed( [ Eltseq(g) : g in Generators(G) ] ));
    Gdualentries := SequenceToSet(FlatFixed( [ Eltseq(g^-1) : g in Generators(G) ] ));
    D3 := PrimeFactors({ Denominator(x) : x in Gentries});
    D4 := PrimeFactors({ Denominator(x) : x in Gentries});

    print "Action of G.";
    IndentPush();
    print D3;
    IndentPop();
    print "Dual action of G.";
    IndentPush();
    print D4;
    IndentPop();

    print "Cherednik coefficients.";
    ReflectionLibrary(~G);
    V := VectorSpace(G);
    D5 := PrimeFactors({ Denominator(CherednikCoefficient(V.i,V.j,s)) : i,j in [1..Dimension(V)], s in G`ReflectionLibraryFlat});
    IndentPush();
    print D5;
    IndentPop();

    print "Total: ";
    IndentPush();
    D:= D1 join D2 join D3 join D4 join D5;
    print D;
    IndentPop();

    return D;

end intrinsic;
