//Compute Calogero-Moser cellular characters for the exceptional complex
//reflection group G_n for the spetsial parameter. Writes output to a Markdown
//file.
n := 27;

//No change from here on
name := "G"*Sprint(n)*"_spetsial";

W:=ComplexReflectionGroup(n);

//Move group and representations to common overfield to avoid problems
//with coercion later.
Representations(~W);
K := CommonOverfield([* BaseRing(Codomain(W`Representations[0][i])) : i in [1..#W`Representations[0]] *]);
Wnew:=ChangeRing(W, K);
Wnew`DBDir := W`DBDir;
Representations(~Wnew);
LiftRepresentationsToCommonBaseField(~Wnew);

//Now compute cellular characters
c:=CherednikParameter(Wnew : Type:="Spetsial");
cellchar := CalogeroMoserCellularCharacters(Wnew,c);
eulerfams := EulerFamilies(Wnew,c);

//Output
str := "# Calogero-Moser cellular characters for G"*Sprint(n)*" for spetsial parameter\n\n";

str *:= "Cellular characters and their decomposition into irreducible characters are listed per (non-trivial) Euler family.\n\n";

for i:=1 to #eulerfams do
  //Skip singleton Euler family
  if #eulerfams[i][1] eq 1 then
    continue;
  end if;

  for j:=1 to #eulerfams[i][1] do
    str *:= "| "*Wnew`CharacterNames[eulerfams[i][1][j]];
  end for;
  str *:= " |";
  str *:= "\n";

  for j:=1 to #eulerfams[i][1] do
    str *:= "| ----";
  end for;
  str *:= " |";
  str *:= "\n";

  for k:=1 to Nrows(cellchar[i][2]) do
    for j:=1 to Ncols(cellchar[i][2]) do
      str *:= "| "*Sprint(cellchar[i][2][k][j]);
    end for;
    str *:= " |";
    str *:= "\n";
  end for;

  str *:= "\n";

end for;

Write(name*".md", str : Overwrite:=true);
