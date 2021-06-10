/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/


//============================================================================
intrinsic GlueSets(X::SetEnum[SetEnum]) -> SetEnum
{The gluing (meet) of set partitions.}

	if IsEmpty(X) then
		return {{}};
	end if;

	Y := {e : e in x, x in X};

	while exists(t){<e,f> : e,f in Y | e ne f and not IsEmpty(e meet f)} do
		e := t[1];
		f := t[2];
		Y diff:={e};
		Y diff:={f};
		Y join:={e join f};
	end while;

	return Y;
end intrinsic;

//============================================================================
intrinsic GlueSets(X::SetEnum[SetIndx]) -> SetIndx
{The gluing (meet) of set partitions.}

	if IsEmpty(X) then
		return {{@@}};
	end if;

	Y := {@e : e in x, x in X@};

	while exists(t){<e,f> : e,f in Y | e ne f and not IsEmpty(e meet f)} do
		e := t[1];
		f := t[2];
		Y diff:={@e@};
		Y diff:={@f@};
		Y join:={@e join f@};
	end while;

	return Y;
end intrinsic;


//============================================================================
intrinsic BlockGraph(Omega::SeqEnum : RestrictionIdeal:=false) -> Assoc
{}
	//number of central character values per character
	N := #Omega[1];

	for i:=2 to #Omega do
		if #Omega[i] ne N then
			error "Wrong size.";
		end if;
	end for;

	//base ring
	P := Parent(Omega[1][1]);

	//root
	root := {@@};
	todo := SequenceToIndexedSet([1..#Omega]);
	while not IsEmpty(todo) do
		i := todo[1];
		iblock := {@i@};
		Diff(~todo, {@i@});
		for j in {i+1..#Omega} meet todo do
			if Omega[i] eq Omega[j] then
				iblock join:={@j@};
				Diff(~todo, {@j@});
			end if;
		end for;
		root join:={@ iblock @};
	end while;

	print "Root is";
	IndentPush();
	print root;
	IndentPop();

	//just take central characters for each block of root, that's enough
	Omega := [ Omega[root[i][1]] : i in [1..#root] ];

	//gluing ideals (ideals where at least two families glue)
	gluingidealsarr := AssociativeArray({@<i,j> : i,j in [1..#Omega] | i lt j@});
	for pair in Universe(gluingidealsarr) do
		i := pair[1];
		j := pair[2];
		gluingidealsarr[<i,j>] := ideal<P|[ Omega[i][k]-Omega[j][k] : k in [1..N]]>;
		if not Type(RestrictionIdeal) eq BoolElt then
			gluingidealsarr[<i,j>] +:= RestrictionIdeal;
		end if;
		Groebner(gluingidealsarr[<i,j>]);
	end for;

	gluingideals := {@ gluingidealsarr[pair] : pair in Keys(gluingidealsarr) | gluingidealsarr[pair] ne ideal<P|1> @};

	//block structures on gluing ideals
	gluingblockstructures := AssociativeArray(gluingideals);
	for I in gluingideals do
		gluings := {@ p : p in Keys(gluingidealsarr) | gluingidealsarr[p] subset I @};
		blocks := {@ @};
		todo := {@i : i in [1..#Omega] @};
		while not IsEmpty(todo) do
			i := todo[1];
			iblock := {@i@};
			while exists(p){ p : p in gluings | p[1] in iblock or p[2] in iblock} do
				iblock join:={@p[1],p[2]@};
				gluings diff:={@p@};
			end while;
			todo diff:=iblock;
			blocks join:={@Sort(iblock)@};
		end while;
		gluingblockstructures[I] := blocks;
	end for;

	//take all sums of gluing ideals until we don't get anything new anymore (not: geometrically, we take intersections)
	stratification := gluingideals;
	stratificationinfo := AssociativeArray(); //will carry ideals involved in the sum
	for i:=1 to #gluingideals do
		stratificationinfo[gluingideals[i]] := {@i@};
	end for;
	lastideals := gluingideals;
	newideals := gluingideals;
	while not IsEmpty(newideals) do
		newideals := {@ @};
		for I in gluingideals do
			for J in lastideals do
				K := I + J; //sum!
				Groebner(K);
				if K ne ideal<P|1> and K notin stratification then
					stratification join:={@K@};
					stratificationinfo[K] := stratificationinfo[I] join stratificationinfo[J];
					newideals join:={@K@};
				end if;
			end for;
		end for;
		lastideals := newideals;
	end while;

	//block structure for each stratification ideal (obtained by gluing the covering gluingblockstructures)
	blockstructures := AssociativeArray(stratification);
	for I in stratification do
		coveringblockstructures := {gluingblockstructures[J] : J in gluingideals | J subset I};
		blockstructures[I] := GlueSets(coveringblockstructures);
	end for;

	print "Gluing ideals: ";
	IndentPush();
	for i:=1 to #gluingideals do
		print i;
		IndentPush();
		print gluingideals[i];
		IndentPop();
	end for;
	IndentPop();

	blockstructuresset := {@root@} join {@ blockstructures[I] : I in Keys(blockstructures) @};

	//sanity check (get different block structure on each stratification ideal as these are sums of gluing ideals (and these were all distinct) )
	assert #blockstructuresset eq #Keys(blockstructures) + 1; //+1 because of root

	//create graph
	str := "digraph G { \n";
	//root
	str *:= "0 [label=\"";
	str *:= "{";
	for j:=1 to #root do
		F := Sort(root[j]);
		str *:= "{";
		for k:=1 to #F do
			str *:= Sprint(F[k]);
			if k lt #F then
				str *:= ",";
			end if;
		end for;
		str *:= "}";
		if j lt #root then
			str *:= ",";
		end if;
	end for;
	str *:="}\"";
	str *:="];\n";

	for i:=1 to #stratification do
		I := stratification[i];
		B := blockstructures[I];
		str *:= Sprint(i)*" [label=\"";
		str *:= "{";
		for j:=1 to #B do
			F := Sort(B[j]);
			str *:= "{";
			for k:=1 to #F do
				str *:= Sprint(F[k]);
				if k lt #F then
					str *:= ",";
				end if;
			end for;
			str *:= "}";
			if j lt #B then
				str *:= ",";
			end if;
		end for;
		str *:="}";
		str *:="\\n";
		L := stratificationinfo[I];
		str *:= "";
		for j:=1 to #L do
			str *:= Sprint(L[j]);
			if j lt #L then
				str *:= "+";
			end if;
		end for;
		str *:= "";
		str *:= "\"";
		str *:="];\n";
	end for;

	//arrows

	//arrows from root to gluing ideals
	for i:=1 to #gluingideals do
		p := Position(stratification, gluingideals[i]);
		str *:= "0 -> "*Sprint(p)*";\n";
	end for;

	//all other arrows
	for i:=1 to #stratification do
		I := stratification[i];
		coverings := {j : j in {1..#stratification} | stratification[j] subset I and stratification[j] ne I};
		for j in coverings do
			str *:= Sprint(j)*" -> "*Sprint(i)*";\n";
		end for;
	end for;


/*
	//rank
	ranks := {#blockloci[B][2] : B in blockstructuresset};
	for r in ranks do
		samerank := {@i : i in [1..#blockstructuresset] | #blockloci[blockstructuresset[i]][2] eq r@};
		str *:= "{rank=same; ";
		for i:=1 to #samerank do
			str *:= Sprint(samerank[i]);
			if i lt #samerank then
				str *:= ",";
			end if;
		end for;
		str *:= "};\n";
	end for;
*/

	str *:= "}";

	file := CHAMP_GetDir()*"Tmp/"*Tempname("blockgraph_");
	Write(file*".dot", str);
	res := System("tred "*file*".dot | dot -T svg > "*file*".svg");
	res := System("open "*file*".svg");

	return blockstructures;

end intrinsic;
