/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Some tools for cuspidal Calogero-Moser families
*/

//============================================================================
intrinsic CuspidalCalogeroMoserFamiliesOld(~H::AlgChe : SaveToDB:=false, UseDB:=true)
{}

	//try to determine the cuspidal families in the fastes way, not computing all the stuff if possible
	CalogeroMoserFamilies(~H);
	W := H`Group;
	CharacterTable(~W);
	Representations(~W); //needed for rigidity check
	SymplecticDoublingFundamentalInvariants(~W);
	inv := W`SymplecticDoublingFundamentalInvariants;
	gen, type := IsGeneric(H);
	if not gen then
		error "Need generic algebra as input";
	end if;

	//pairs where poisson bracket is of degree zero
	deg0pairs := { <i,j> : i,j in [1..#inv] | i lt j and -Bidegree(inv[i])[1] + Bidegree(inv[i])[2] - Bidegree(inv[j])[1] + Bidegree(inv[j])[2] eq 0 };
	//find already computed poisson brackets
	computedpoissonbrackets := {};
	if not assigned H`CenterGeneratorsPoissonBrackets then
		H`CenterGeneratorsPoissonBrackets := AssociativeArray({<k,l> : k,l in [1..#inv]});
	end if;
	for pair in deg0pairs do
		i := pair[1];
		j := pair[2];
		if gen and assigned type and CHAMP_ExistsInDB(W`DBDir, "CherednikGeneric"*type*"CenterGeneratorsPoissonBracket"*Sprint(i)*"_"*Sprint(j)) then
				computedpoissonbrackets join:={<i,j>};
		end if;
	end for;

	print "There are "*Sprint(#inv)*" center generators";
	print "There are "*Sprint(#deg0pairs)*" pairs of distinct center generators whose Poisson bracket is of degree zero";
	print "There are "*Sprint(#computedpoissonbrackets)*" of these Poisson brackets in the data base already";

	//GENERIC CASE
	print "Dealing with generic case";
	//singleton families are not cuspidal, so can ignore them
	fams := {fam : fam in H`CalogeroMoserFamilies[1] | #fam gt 1};
	print "\tThere are "*Sprint(#fams)*" non-singleton families";
	cuspfams := {};
	c := H`cParameter;
	for fam in fams do
		//if we find a rigid character in fam, it's cuspidal
		print "\tDealing with family "*Sprint(fam);
		IndentPush();
		print "\t\tSearching for rigid representation";
		IndentPush();
		for i in fam do
			if IsRigid(W`Representations[0][i], c) then
				print "\t\t\tFound a rigid character, so family is cuspidal";
				IndentPop();
				IndentPop();
				cuspfams join:={fam};
				continue fam;
			end if;
		end for;
		print "\t\t\tThere are no rigid representations.";
		IndentPop();
		print "\t\tContinuing with already computed Poisson brackets";
		for pair in computedpoissonbrackets do
			CenterGeneratorsPoissonBracket(~H, pair[1],pair[2]);
			z := H`CenterGeneratorsPoissonBrackets[pair];
		 	for i in fam do
		 		alpha := ActionOnVermaModule(W`CharacterTable[i], z);
		 		if alpha ne 0 then
		 			print "\t\t\tFound a representation and Poisson bracket acting by non-zero, so family is not cuspidal";
		 			continue fam;
		 		end if;
		 	end for;
		 end for;
		print "\t\tAll computed Poisson brackets act by non-zero on any representation";
		if #computedpoissonbrackets eq #deg0pairs then
			print "\t\tSince we have computed all Poisson brackets, the family is cuspidal";
			cuspfams join:={fam};
			continue fam;
		end if;
		print "\t\tStill have to decide whether family is cuspidal. Computing more Poisson brackets.";
		while #computedpoissonbrackets ne #deg0pairs do
			print "\t\t\tLoading or computing remaining center generators";
			remainingpairs := deg0pairs diff computedpoissonbrackets;
			gens := { pair[1] : pair in remainingpairs} join { pair[2] : pair in remainingpairs };
			for i in gens do
				CenterGenerator(~H, i);
			end for;
			//find "small" pair
			remainingpairs := SetToSequence(remainingpairs);
			complexities := [ #Sprint(H`CenterGenerators[pair[1]]) + #Sprint(H`CenterGenerators[pair[2]]) : pair in remainingpairs];
			goodpair := remainingpairs[Position(complexities, Minimum(complexities))];
			print "\t\t\t\Computing Poisson bracket for "*Sprint(goodpair);
			CenterGeneratorsPoissonBracket(~H, goodpair[1], goodpair[2]);
			computedpoissonbrackets join:={goodpair};
			print "\t\t\tChecking whether action is non-zero";
			for i in fam do
		 		alpha := ActionOnVermaModule(W`CharacterTable[i], H`CenterGeneratorsPoissonBrackets[goodpair]);
		 		if alpha ne 0 then
		 			print "\t\t\tAction is non-zero, so family is not cuspidal";
		 			continue fam;
		 		end if;
		 	end for;
		end while;
		print "\t\tAll actions are zero, so family is cuspidal";
		cuspfams join:={fam};
	end for;

end intrinsic;

intrinsic CuspidalCalogeroMoserFamilies(~H::AlgChe : SaveToDB:=false, UseDB:=true)
{}

	if assigned H`CuspidalCalogeroMoserFamilies then
		return;
	end if;

	CalogeroMoserFamilies(~H);
	W := H`Group;
	CharacterTable(~W);
	SymplecticDoublingFundamentalInvariants(~W);
	inv := W`SymplecticDoublingFundamentalInvariants;
	gen, type := IsGeneric(H);
	if not gen then
		error "Need generic algebra as input";
	end if;

	//pairs where poisson bracket is of degree zero
	deg0pairs := { <i,j> : i,j in [1..#inv] | i lt j and -Bidegree(inv[i])[1] + Bidegree(inv[i])[2] - Bidegree(inv[j])[1] + Bidegree(inv[j])[2] eq 0 };

	CenterGenerators(~H);

	print "Computing Poisson brackets (there are "*Sprint(#deg0pairs)*" degree-0 pairs)";
	for pair in deg0pairs do
		print pair;
		CenterGeneratorsPoissonBracket(~H, pair[1],pair[2]);
		z := H`CenterGeneratorsPoissonBrackets[pair];
	end for;

	print "Computing actions of Poisson brackets on Verma modules of generic families";
	//We compute \Omega_F({z_i,z_j}) for the generic families F.
	//We index the generic families by the minimum character index they
	//contain.
	genfamsmins := { Minimum(fam) : fam in H`CalogeroMoserFamilies[1] };
	genactions := AssociativeArray(genfamsmins);
	for i in genfamsmins do
		genactions[i] := {};
		for pair in deg0pairs do
			z := H`CenterGeneratorsPoissonBrackets[pair];
			alpha := ActionOnVermaModule(W`CharacterTable[i],z);
			genactions[i] join:= {alpha};
		end for;
	end for;

	R := Universe(H`CalogeroMoserFamilies);
	H`CuspidalCalogeroMoserFamilies := AssociativeArray(R);

	for hyp in Keys(H`CalogeroMoserFamilies) do
		H`CuspidalCalogeroMoserFamilies[hyp] := {};
		RH := R/ideal<R|hyp>;
		for fam in H`CalogeroMoserFamilies[hyp] do
			i := Minimum(fam);
			//\Omega_F^c is given by the reduction of the generic \Omega_F
			for alpha in genactions[i] do
				if RH!(R!alpha) ne 0 then
					continue fam;
				end if;
			end for;
			H`CuspidalCalogeroMoserFamilies[hyp] join:={fam};
		end for;
	end for;
end intrinsic;

intrinsic CuspidalCalogeroMoserFamilies(H::AlgChe : SaveToDB:=false, UseDB:=true) -> Assoc
{}

	CuspidalCalogeroMoserFamilies(~H);
	return H`CuspidalCalogeroMoserFamilies;

end intrinsic;
