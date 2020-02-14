/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


/*
	Branching of irreducible characters of a subgroup
*/


//============================================================================
intrinsic BranchingIndex(G::Grp, H::Grp) -> RngIntElt
{}
	assert H subset G;

	CharacterTable(~G);
	CharacterTable(~H);
	b := 0;
	for chi in H`CharacterTable do
		chiG := Induction(chi, G);
		D := Decomposition(chiG);
		bchi := Integers()!(&+D);
		if b eq 0 or bchi lt b then
			b := bchi;
		end if;
	end for;
	return b;

end intrinsic;

//============================================================================
intrinsic ParabolicBranchingIndex(G::Grp) -> RngIntElt
{}

	ParabolicSubgroups(~G);
	b := 0;
	for i in G`ParabolicSubgroups do
		H := G`Subgroups[i]`subgroup;
		if #H eq #G then
			continue;
		end if;
		bH := BranchingIndex(G,H);
		if b eq 0 or bH lt b then
			b := bH;
		end if;
	end for;

	return b;

end intrinsic;
