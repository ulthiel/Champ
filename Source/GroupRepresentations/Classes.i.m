//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Handling conjugacy classes
//
//==============================================================================

//==============================================================================
// Magma likes to reorder classes when manually settimg them (because it wants
// a particular ordering).
// The following functions help to determine the permutation.
// This is important when importing classes from elsewhere.
//==============================================================================
intrinsic ClassPermutation(G::Grp, Y::SeqEnum[MonStgElt]) -> SeqEnum
{If G has ClassNames assigned and Y is a sequence of strings (with class names), return the permutation sigma such that Y[sigma[i]] = ClassNames[i].}

	if assigned G`ClassNames then
		return FindPermutation(G`ClassNames, Y);
	else
		error "Group does not have ClassNames assigned.";
	end if;

end intrinsic;

intrinsic ClassPermutation(G::Grp, Y::SeqEnum[GrpElt]) -> SeqEnum
{If G has Classes assigned and Y is a sequence of group elements (giving representatives of the classes), return the permutation sigma such that Y[sigma[i]] = ClassRepresentative(G,i).}

	if assigned G`Classes then
		class_reps := [ x[3] : x in G`Classes ];
		sigma := [];
		remaining_indeces := {1..#Y};

		for i:=1 to #G`Classes do
			x := ClassRepresentative(G,i);

			//we now want to find the element in Y which is conjugate to x
			//we should first test equality before doing anything about conjugacy,
			//maybe we're lucky.

			//maybe we're really lucky
			if x eq Y[i] then
				Append(~sigma,i);
				remaining_indeces diff:={i};
				continue;
			end if;

			//next, do a full search with equality
			for j in remaining_indeces do
				if x eq Y[j] then
					Append(~sigma,j);
					remaining_indeces diff:={j};
					continue i;
				end if;
			end for;

			//finally, we need to work
			for j in remaining_indeces do
				if IsConjugate(G, x, Y[j]) then
					Append(~sigma,j);
					remaining_indeces diff:={j};
					continue i;
				end if;
			end for;

			error "There is no permutation";

		end for;

		// test if we really found a permutation
		if not #sigma eq #G`Classes and SequenceToSet(sigma) eq {1..#G`Classes} then
			error "There is no permutation";
		else
			return sigma;
		end if;
	else
		error "Group does not have Classes assigned.";
	end if;

end intrinsic;

intrinsic ClassPermutation(G::Grp, Y::SeqEnum[SeqEnum[RngIntElt]]) -> SeqEnum
{If G has ClassWords or Classes assigned and Y is a sequence of integers (giving words for class representatives), return the permutation sigma such that Y[sigma[i]] = ClassNames[i].}

	if assigned G`ClassWords then
		return FindPermutation(G`ClassWords, Y);
	elif assigned G`Classes then
		elts := [ WordToElement(G,x) : x in Y ];
		return ClassPermutation(G, elts);
	else
		error "Group does neither have ClassWord nor Classes assigned.";
	end if;

end intrinsic;

//==============================================================================
intrinsic ClassRepresentatives(G::Grp) -> SeqEnum
{}

	C := Classes(G);
	return [ ClassRepresentative(G,i) : i in [1..#G`Classes] ];

end intrinsic;

//==============================================================================
intrinsic SaveClasses(G::Grp)
{Saves all available information about conjugacy classes in the database.}

	if not assigned G`DBName then
		error "The group needs to have a the attribute DBName assigned.";
	end if;

	code := "ClassesRecordFormat := recformat<NumberOfClasses : RngIntElt, ClassOrders : SeqEnum, ClassLengths : SeqEnum, ClassRepresentatives : SeqEnum, ClassWords : SeqEnum, ClassNames : SeqEnum>;\n";

	code *:= "ClassesRecord := rec<ClassesRecordFormat|>;\n";

	code *:= Sprintf("ClassesRecord`NumberOfClasses := %o;\n", NumberOfClasses(G));

	code *:= Sprintf("ClassesRecord`ClassOrders := %o;\n", Sprint([x[1] : x in G`Classes], "Magma"));

	code *:= Sprintf("ClassesRecord`ClassLengths := %o;\n", Sprint([x[2] : x in G`Classes], "Magma"));

	if assigned G`ClassWords then
		code *:= Sprintf("ClassesRecord`ClassWords := %o;\n", Sprint(G`ClassWords, "Magma"));
	else
		code *:= Sprintf("ClassesRecord`ClassRepresentatives := %o;\n", [Sprintf(x,"Magma") : x in ClassRepresentatives(G) ]);
	end if;

	if assigned G`ClassNames then
		code *:= Sprintf("ClassesRecord`ClassNames := %o;\n", Sprint(G`ClassNames, "Magma"));
	end if;

	code *:= "return ClassesRecord";

	CHAMP_SaveToDB(code, "Classes", G`DBName);

end intrinsic;

intrinsic LoadClasses(~W::Grp)
{Loads all available information about conjugacy classes from the database.}

	;

end intrinsic;
