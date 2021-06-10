/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

procedure SieveBydbPair(G, ~labels)

    CharacterTable(~G);
    FakeDegrees(~G);

    dbpairs := [ <G`CharacterTable[i][1],Valuation(G`FakeDegrees[i])> : i in [1..#G`CharacterTable]];
    baddbpairs := {};

    numberofgoodcharacters := 0;

    for i in [1..#G`CharacterTable] do
        m := NumberOfOccurrences(dbpairs, dbpairs[i]);
        if m eq 1 then
            labels[i] := "\\phi_{"*Sprint(dbpairs[1])*","*Sprint(dbpairs[2])*"}";
        else
            baddbpairs join:={dbpairs[i])};
            numberofgoodcharacters+:=1;
        end if;
    end for;

    baddbpairs := Sort(SetToSequence(baddbpairs));

    printf "\\noindent\\textbf{Sieve by db-pairs:} %o characters (leaving %o bad) are already uniquely determined by their db-pair. The bad db-pairs with multiplicity are:", numberofgoodcharacters, #G`CharacterTable - numberofgoodcharacters;

    printf "\\begin{center}\n";
    for i:=1 to #baddbpairs do
        x := baddbpairs[i];
        "$(%o,%o)^{%o}$ ", x[1], x[2], NumberOfOccurrences(dbpairs, x);
        if i lt #baddbpairs then
            printf ", \\; ";
        end if;
    end for;
    printf "\\;.\n\\end{center}\n";

end procedure;

procedure CRGLabelingForGroup(n);

    G:=ExceptionalComplexReflectionGroup(n);
    CharacterTable(~G);
    FakeDegrees(~G);
    labels := ["\\phi" : i in [1..#G`CharacterTable]];
    SieveBydbPair(G, ~labels);

end procedure;
