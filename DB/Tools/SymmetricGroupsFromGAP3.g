#Read data about symmetric complex reflection groups from gap3

for i in [2..22] do
    
    G:=ComplexReflectionGroup(1,1,i);
	PrintTo(ConcatenationString("GAP3_Data/G", String(i)), G.matgens); 
	PrintTo(ConcatenationString("GAP3_Data/C",String(i)), WordsClassRepresentatives(G)); 
	PrintTo(ConcatenationString("GAP3_Data/X", String(i)), CharTable(G).irreducibles); 	
	#PrintTo(ConcatenationString("GAP3_Data/var", String(i)), var(i));
    PrintTo(ConcatenationString("GAP3_Data/XN", String(i)), CharNames(CharTable(G)));
    
    q := X( Cyclotomics ); q.name := "q";
    PrintTo(ConcatenationString("GAP3_Data/fake", String(i)), FakeDegrees( G, q ));
    
	#for d in DivisorsInt(Size(G)) do
	#	if IsPrime(d) then
	#		PrintTo(ConcatenationString("GAP3_Data/hyperplanes", String(i), "_", String(d)), EssentialHyperplanes(G,d));
	#	fi;
	#od;
	#PrintTo(ConcatenationString("GAP3_Data/Rou", String(i)), DisplayAllBlocks(G));
	
	if i < 10 then
		PrintTo(ConcatenationString("GAP3_Data/R", String(i)), Representations(G));
	fi;
	
od;