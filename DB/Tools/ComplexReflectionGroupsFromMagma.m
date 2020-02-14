load "../../CHAMP.m";

for i:=4 to 37 do
    print i;
    
    str := "/*";
    str *:= "\n";
    str *:= "    Exceptional complex reflection group G"*Sprint(i)*" as obtained from ShephardTodd("*Sprint(i)*") in Magma "*GetVersionString()*".";
    str *:= "\n\n";
    str *:= "    Date: "*Date();
    str *:= "\n";
    str *:= "*/\n"; 
    
    G := ShephardTodd(i);
    str *:= Sprint(G, "Magma");
    
    CHAMP_SaveToDB(str, "GrpMat/G"*Sprint(i)*"_Magma/GrpMat.m");
    
    
    str := "/*";
    str *:= "\n";
    str *:= "    Dual of exceptional complex reflection group G"*Sprint(i)*" as obtained from DualGroup(ShephardTodd("*Sprint(i)*")) in Magma "*GetVersionString()*".";
    str *:= "\n\n";
    str *:= "    Date: "*Date();
    str *:= "\n";
    str *:= "*/\n"; 
    
    G := DualGroup(ShephardTodd(i));
    str *:= Sprint(G, "Magma");
    
    CHAMP_SaveToDB(str, "GrpMat/G"*Sprint(i)*"_Magma_Dual/GrpMat.m");
    
end for;