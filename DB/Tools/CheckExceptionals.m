//Check if GAP3 data matches CHAMP data for ECRGs.

E := func<n | RootOfUnity(n)>;

for i:=4 to 37 do
	
	G:=ExceptionalComplexReflectionGroup(i); 
	CharacterTable(~G); 
	mynames:=G`CharacterNames; 
	for j:=1 to #mynames do; 
		mynames[j]:=Replace(mynames[j], "\\\\", ""); 
		mynames[j]:=Replace(mynames[j], "_", ""); 
	end for; 
	names:=eval Read("GAP3_Data/XN"*Sprint(i));
	
	print "###G"*Sprint(i);
	
	if names ne mynames then
		print "The following labels changed: \n";
		firstchange := true;
		for k:=1 to #names do
			if names[k] ne mynames[k] then
				if firstchange then
					firstchange := false;
					print "| CHAMP | gap3-jm5 |";
					print "| ------------- |:-------------:|";
				end if;
				if Replace(names[k], "'", "") ne Replace(mynames[k], "'", "") then
					error "Unallowed change!";
				end if;
				print "| "*mynames[k] * "\t | \t" * names[k];
			end if;
		end for;		
	else
		print "Nothing changed.";
	end if;
	
	gapclasses := eval Read("GAP3_Data/C"*Sprint(i));
	gapcharacters := eval Read("GAP3_Data/X"*Sprint(i));
	if #G`CharacterTable ne #gapcharacters then
		error "Number of characters not correct!";
	end if;
	for j:=1 to #G`CharacterTable do
		chi := G`CharacterTable[j];
		for k:=1 to #gapclasses do
			w := Reverse(gapclasses[k]);
			g := WordToElement(G,w);
			if chi(g) ne gapcharacters[j][k] then
				error "Character "*Sprint(j)*" not correct.";
			end if;
		end for;
	end for;

	R<q> := PolynomialRing(Integers());
	fakes := eval Read("GAP3_Data/fake"*Sprint(i));
	FakeDegrees(~G);
	assert fakes eq G`FakeDegrees;
	
	
	print "\n\n";
	
end for;