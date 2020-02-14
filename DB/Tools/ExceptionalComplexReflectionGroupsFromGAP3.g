#Read data about exceptional complex reflection groups from gap3

for i in [4..37] do
    
    G:=ComplexReflectionGroup(i);
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
	
od;

#representations

#for i in [35] do 
#    G:=ComplexReflectionGroup(i);  
#    PrintTo(ConcatenationString("GAP3_Data/R", String(i)), Representations(G));
#od;