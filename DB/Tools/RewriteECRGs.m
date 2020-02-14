for i:=4 to 37 do
	str := Read(CHAMP_GetDBDir()*"GrpMat/G"*Sprint(i)*"_CHEVIE/GrpMat.m");
	G := eval str;
	pos := Position(str, "CET 2014");
	assert pos ne 0;
	str := str[1..pos+7];
	str *:= "\n\n    Rewritten with Magma Sprint on "*Date()*"\n*/\n";
	str *:= Sprint(G, "Magma");
	print str;
end for;
