/*
    CHAMP (CHerednik Algebra Magma Package)
    Copyright (C) 2013, 2014 Ulrich Thiel
    Licensed under GNU GPLv3, see COPYING.
    thiel@mathematik.uni-stuttgart.de
*/
/*
    Testing complex reflection groups.
*/
print "Running self check \"CRG\"";
zeit := Cputime();

for i:=4 to 37 do
    G := ExceptionalComplexReflectionGroup(i);
    CharacterTable(~G:Check:=true);
        
    for j:=1 to #G`CharacterTable do
    	if not IsIrreducible(G`CharacterTable[j]) then
    		error "Non-irreducible character found.";
    	end if;
    end for;
    
    if #SequenceToSet(G`CharacterTable) ne #G`CharacterTable then
    	error "Isomorphic characters in character table.";
    end if;
    
    //check (d,b) labeling
    Degrees(~G);
    ReflectionLibrary(~G);
    FakeDegrees(~G : UseDB:=false, Method:="Formula");
    PhiNames(~G);
    unprimed := [ Replace(G`CharacterNames[i], "'", "") : i in [1..#G`CharacterNames]];
    assert unprimed eq G`PhiNames;  //checks labelings and fake degree computation
    
    //check database fake degrees
    dbfakes := G`FakeDegrees;
    G := ExceptionalComplexReflectionGroup(i);
    FakeDegrees(~G : UseDB := true);
    assert dbfakes eq G`FakeDegrees;
    
    //check direct computation of fake degrees
    //this is in particular a test of DecompositionInGradedGrothendieckgroup for ModGrGrp
    //we cannot check all groups as this gets really big
    if i in {4,5,6,7,8,9,10,11,12,13,14,15,16,17,20,21,22,23,24,25} then
        G:=ExceptionalComplexReflectionGroup(i);
        FakeDegrees(~G:Method:="Direct", UseDB:=false);
        assert G`FakeDegrees eq dbfakes;
    end if;    
    
    G:=ExceptionalComplexReflectionGroup(i);
    CharacterTable(~G);
    
    //check representations up to G32
    if i lt 32 then
        Representations(~G:Check:=true);
        for j:=1 to #G`CharacterTable do
        	if Character(G`Representations[0][j]) ne G`CharacterTable[j] then
        		error "Character Representation Problem with G"*Sprint(i);
        	end if;
        end for;
    end if;

    if i in {4,5,6,7,8,9,10,12,13,14,15,16,20,22,23,24} then
    	res := Gordon(G);
    end if;
    
    res := RouquierFamilies(G);

    printf "\b\b\b%o ", i;

end for;

//Checking realizations of symmetric reflection groups
for i:=2 to 22 do; 
	//print i;
	G:=SymmetricReflectionGroup(i:UseDB:=true); 
	H:=SymmetricReflectionGroup(i:UseDB:=false); 
	assert [G.j:j in [1..Ngens(G)]] eq [H.j:j in [1..Ngens(H)]]; 
end for;

IndentPush();
print "Time: "*Sprint(Cputime(zeit));
IndentPop();

