reps := true; //skip representations

load "../../CHAMP.m";


E := func<n | RootOfUnity(n)>;

groups := [];


for j in [<2,1,2>, <2,1,3>, <2,1,4>, <2,1,5>] do
    print j;
    i:= Sprint(j[1])*"_"*Sprint(j[2])*"_"*Sprint(j[3]);
    //Gi
    Gmat := eval Read("GAP3_Data/G"*Sprint(i));
    Gpre := MatrixGroup(Gmat);
    G := MatrixGroup<Dimension(Gpre), BaseRing(Gpre) | [ Transpose(Gpre.i) : i in [1..Ngens(Gpre)]]>;
    Append(~groups, G);
    str := "/*";
    str *:= "\n";
    str *:= "    Exceptional complex reflection group G("*Sprint(i[1])*","*Sprint(i[2])*","*Sprint(i[3])*") as obtained from the transposed of ComplexReflectionGroup("*Sprint(i[1])*","*Sprint(i[2])*","*Sprint(i[3])*") in GAP3v5.";
    str *:= "\n\n";
    str *:= "    Date: "*Date();
    str *:= "\n";
    str *:= "*/\n"; 
  
    str *:= Sprint(G, "Magma");
    
    CHAMP_SaveToDB(str, "GrpMat/G"*Sprint(i)*"_CHEVIE/GrpMat.m");
    
    //dual of Gi
    str := "/*";
    str *:= "\n";
    str *:= "    Dual of exceptional complex reflection group G("*Sprint(i[1])*","*Sprint(i[2])*","*Sprint(i[3])*") as obtained from the transposed of ComplexReflectionGroup("*Sprint(i[1])*","*Sprint(i[2])*","*Sprint(i[3])*") in GAP3v5.";
    str *:= "\n\n";
    str *:= "    Date: "*Date();
    str *:= "\n";
    str *:= "*/\n"; 
  
    str *:= Sprint(DualGroup(G), "Magma");
    
    CHAMP_SaveToDB(str, "GrpMat/G"*Sprint(i)*"_CHEVIE_Dual/GrpMat.m");
    
    //character data of Gi
    Cpre := eval Read("GAP3_Data/C"*Sprint(i));
    C := [ Reverse(Cpre[j]) : j in [1..#Cpre] ]; //classes
    X := eval Read("GAP3_Data/X"*Sprint(i)); //character table
    N := eval Read("GAP3_Data/XN"*Sprint(i)); //character names
    
    //rewrite N
    for i:=1 to #N do
        str := N[i];
        str := Replace(str, "{", "_{");
        str := Replace(str, "phi", "\\phi");
        N[i] := str;
    end for;
    
    str := "/*";
    str *:= "\n";
    str *:= "    Character data of exceptional complex reflection group G"*Sprint(i)*" as obtained by reversing the lists in WordsClassRepresentatives in GAP3v4 (Reversing because we transposed the group!), using CharTable, and using CharNames.";
    str *:= "\n\n";
    str *:= "    Date: "*Date();
    str *:= "\n";
    str *:= "*/\n";
    
    str *:= Sprint(<C,X,N>, "Magma");

    CHAMP_SaveToDB(str, "GrpMat/G"*Sprint(i)*"_CHEVIE/CharacterData.m");

        //Gi
        Rpre := eval Read("GAP3_Data/R"*Sprint(i));
        R := [* *];
        for j:=1 to #Rpre do
            entries := SequenceToSet(Flat(Rpre[j]));
            K := MinimalCyclotomicField(entries);
            mats := [ Transpose(Matrix(K, Rpre[j][k])) : k in [1..Ngens(G)] ];
            Append(~R, mats);
        end for;

        str := "/*";
        str *:= "\n";
        str *:= "    Transposed of irreducible representations of exceptional complex reflection group G"*Sprint(i)*" as obtained by Representations in GAP3v4.";
        str *:= "\n\n";
        str *:= "    Date: "*Date();
        str *:= "\n";
        str *:= "*/\n";
         
        str *:= Sprint(R, "Magma");    
    
        CHAMP_SaveToDB(str, "GrpMat/G"*Sprint(i)*"_CHEVIE/Representations_0.m");

end for;
