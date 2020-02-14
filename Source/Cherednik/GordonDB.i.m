/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

function LaurentSeriesRingEltToFractionFieldElt(K, f)

	x := K.1;
	coeffs, val, den := Coefficients(f);
	return ArraySum([ (K!coeffs[i])*x^(val-1+i) : i in [1..#coeffs] | coeffs[i] ne 0] : ZeroElement:=Zero(K));

end function;

//============================================================================
intrinsic WriteGordonRecord(G::GrpMat, res::Rec) -> MonStgElt
{}

	R := res`ParameterRing;
	CharacterTable(~G);

	//Header
	str := "/*\n";
    str *:= "\tData about the representation theory of restricted rational Cherednik algebra\n";
    str *:= "\tDate: "*Sprint(Date())*"\n";
    str *:= "\tVersion: "*Sprint(CHAMP_GetVersion())*"\n*/\n";

	str *:= "\n//Definition of record format\n";

	str *:= "RRCARec := recformat<ParameterRing, BlGen, DecGenStratification, Data>;\n";

    str *:= "RRCADataRec := recformat<SimpleDims, SimplePSeries, SimpleGModStruct, SimpleGradedGModStruct, VermaDecomposition, CMFamilies, CuspidalCMFamilies, VermaGradedDecomposition>;\n";

    str *:= "\n//Base ring for Poincare series\n";
    str *:= "P<q> := RationalFunctionField(Rationals());\n";

    P<q> := RationalFunctionField(Rationals());

    str *:= "\n//Graded Grothendieck group\n";
    str *:= "V := KSpace(P, "*Sprint(#G`CharacterTable)*");\n";

    str *:= "\n//Grothendieck group\n";
    str *:= "W := RModule(Integers(), "*Sprint(#G`CharacterTable)*");\n";

    str *:= "\n//Base ring for hyperplanes\n";
    str *:= "Q := PolynomialRing(Rationals(), "*Sprint(Ngens(R))*");\n";

    //clean old naming
    names := Names(R);
    for i:=1 to #names do
    	names[i] := Replace(names[i], "_{", "");
    	names[i] := Replace(names[i], ",", "_");
    	names[i] := Replace(names[i], "}", "");
    end for;

    AssignNames(~R, names);

    str *:= "AssignNames(~Q, [";
    for i:=1 to #names do
    	str *:= "\""*Sprint(names[i])*"\"";
    	if i lt #names then
    		str *:= ",";
    	end if;
    end for;
    str *:= "]);\n";

    for i:=1 to Ngens(R) do
        str *:= Sprint(names[i])*" := Q."*Sprint(i)*";\n";
    end for;

    str *:= "\n//Now, the results\n";

    str *:= "Results := rec<RRCARec|>;\n\n";

    str *:= "Results`ParameterRing := Q;\n\n";

    if assigned res`BlGen then
    	str *:= "//BlGen\n";
    	str *:= "Results`BlGen := [\n";
    	N := #res`BlGen;
    	count := 0;
    	for H in res`BlGen do
			count +:= 1;
			str *:= "\t"*Sprint(H);
			if count lt N then
				str *:= ",\n";
			end if;
    	end for;
    	str *:= "\n\t];\n\n";
    end if;

   	str *:= "Results`Data := AssociativeArray(PowerSet(Q));\n";

	currHcount := 0;
    for H in Keys(res`Data) do

    	currHcount +:= 1;

		Hstr := "{";
		N := #H;
		count := 0;
		for h in H do
			count +:=1;
			Hstr *:= Sprint(h);
			if count lt N then
				Hstr *:= ", ";
			end if;
		end for;
		Hstr *:= "}";

		str *:= "\n//Stratum "*Hstr*"\n";

		str *:= "Results`Data["*Hstr*"] := rec<RRCADataRec|>;\n";

		if assigned res`Data[H]`SimpleDims then

			str *:= "Results`Data["*Hstr*"]`SimpleDims := "*Sprint(res`Data[H]`SimpleDims)*";\n";

			str *:= "Results`Data["*Hstr*"]`SimplePSeries :=\n[\n";
			for i:=1 to #res`Data[H]`SimplePSeries do
				str *:= "\t"*Sprint(res`Data[H]`SimplePSeries[i]);
				if i lt #res`Data[H]`SimplePSeries then
					str *:= ",\n";
				end if;
			end for;
			str *:= "\n];\n";

		end if;

		if assigned res`Data[H]`SimpleGModStruct then

			str *:= "Results`Data["*Hstr*"]`SimpleGModStruct := \n[\n";
			for i:=1 to #res`Data[H]`SimpleGModStruct do
				str *:= "\tW!"*Sprint(Eltseq(res`Data[H]`SimpleGModStruct[i]))*"";
				if i lt #res`Data[H]`SimpleGModStruct then
					str *:= ",\n";
				end if;
			end for;
			str *:= "\n];\n";

		end if;

		if assigned res`Data[H]`SimpleGradedGModStruct then

			str *:= "Results`Data["*Hstr*"]`SimpleGradedGModStruct := \n[\n";
			for i:=1 to #res`Data[H]`SimpleGradedGModStruct do
				str *:= "\tV!"*Replace(Replace(Sprint(Eltseq(res`Data[H]`SimpleGradedGModStruct[i])), "\n", ""), " ", "")*"";
				if i lt #res`Data[H]`SimpleGradedGModStruct then
					str *:= ",\n";
				end if;
			end for;
			str *:= "\n];\n";

		end if;

		if assigned res`Data[H]`VermaDecomposition then

			str *:= "Results`Data["*Hstr*"]`VermaDecomposition := \n[\n";
			for i:=1 to #res`Data[H]`VermaDecomposition do
				str *:= "\tW!"*Sprint(Eltseq(res`Data[H]`VermaDecomposition[i]))*"";
				if i lt #res`Data[H]`VermaDecomposition then
					str *:= ",\n";
				end if;
			end for;
			str *:= "\n];\n";

		end if;

		if assigned res`Data[H]`VermaGradedDecomposition then
			str *:= "Results`Data["*Hstr*"]`VermaGradedDecomposition := \n[\n";
				for i:=1 to #res`Data[H]`VermaGradedDecomposition do
					str *:= "\tV!"*Replace(Replace(Sprint(Eltseq(res`Data[H]`VermaGradedDecomposition[i])), "\n", ""), " ", "")*"";
					if i lt #res`Data[H]`VermaGradedDecomposition then
						str *:= ",\n";
					end if;
				end for;
			str *:= "\n];\n";
		end if;

		if assigned res`Data[H]`CMFamilies then

			str *:= "Results`Data["*Hstr*"]`CMFamilies := \n{\n";
			fam := SetToSequence(res`Data[H]`CMFamilies);
			for i:=1 to #fam do
				str *:= "\t"*Sprint(fam[i]);
				if i lt #fam then
					str *:= ",\n";
				end if;
			end for;
			str *:= "\n};\n";

		end if;

	end for;

	print "";

	str *:= "return Results";

	return str;

    //CHAMP_SaveToDB(str, G`DBDir, "Cherednik/Gordon");

end intrinsic;

//=============================================================================
/*
intrinsic CleanGordonDBs()
{}

	for n in [4,5,6,7,8,9,10,12,13,14,15,16,20,22,23,24] do
		print n;
		IndentPush();
		G := ExceptionalComplexReflectionGroup(n);
		res := Gordon(G);

		for H in Keys(res) diff {1} do
			if not assigned res[H]`SimpleDims or res[H] eq res[1] then
				Remove(~res, H);
				print H;
			end if;
		end for;

		WriteGordonRecord(G,res);

		IndentPop();
		print "";
	end for;

end intrinsic;
*/

//============================================================================
function FixCharacterNamesHTML(str)

	str := Replace(str, "\\\\phi", "&#981;");
	str := Replace(str, "_{", "<sub>");
	str := Replace(str, "}", "</sub>");
	return str;

end function;

//============================================================================
function PrintGordonHTML(G, res, h)

	str := "";

	if assigned res[h]`CMFamilies then

		str *:= "<h3>Non-singleton Calogero&ndash;Moser k-families</h3>\n";

		nonsigneltonfams := ([fam : fam in res[h]`CMFamilies | #fam ne 1]);

		if #nonsigneltonfams eq 0 then
			str *:= "There are none.";
		end if;

		for j:=1 to #nonsigneltonfams do
			fam := SetToSequence(nonsigneltonfams[j]);
			if #fam eq 1 then
				continue;
			end if;
			famstr := "";
			for i:=1 to #fam do
				famstr *:= G`CharacterNames[fam[i]];
				if i lt #fam then
					famstr *:= ", &nbsp;";
				end if;
			end for;
			famstr := FixCharacterNamesHTML(famstr);
			str *:= "{"*famstr*"}";
			if j lt #nonsigneltonfams then
				str *:= ", &nbsp; ";
			end if;
		end for;
	end if;

	str *:= "<br><br>\n";


	rou := RouquierFamilies(G);
	P := Universe(Keys(rou));
	hsharp := P!NormalizeRationalHyperplaneEquation(G`MartinoSharp(h));
	if hsharp in Keys(rou) then
		roufams := rou[hsharp];
	else
		roufams := rou[1];
	end if;
	str *:= "<h3>Non-singleton Rouquier k<sup>#</sup>-families</h3>\n";

	if roufams eq res[h]`CMFamilies then
		str *:= "Same as Calogero&ndash;Moser k-families.";
	else
		nonsigneltonfams := ([fam : fam in roufams | #fam ne 1]);

		if #nonsigneltonfams eq 0 then
			str *:= "There are none.";
		end if;

		for j:=1 to #nonsigneltonfams do
			fam := SetToSequence(nonsigneltonfams[j]);
			if #fam eq 1 then
				continue;
			end if;
			famstr := "";
			for i:=1 to #fam do
				famstr *:= G`CharacterNames[fam[i]];
				if i lt #fam then
					famstr *:= ", &nbsp;";
				end if;
			end for;
			famstr := FixCharacterNamesHTML(famstr);
			str *:= "{"*famstr*"}";
			if j lt #nonsigneltonfams then
				str *:= ", &nbsp; ";
			end if;
		end for;
	end if;

	str *:= "<br><br>\n";

	if assigned res[h]`SimpleDims and assigned res[h]`VermaDecomposition and assigned res[h]`SimplePSeries then

		str *:= "<h3>Dimensions, Poincar&eacute; series and diagonal Verma multiplicities of the simple modules</h3>\n";

		str *:= "<table>\n";

		str *:= "<tr>\n";
		str *:= "<td> &#981; </td> <td>dim L(&#981;) </td> <td> [&Delta;(&#981;) : L(&#981;)] </td> <td> P<sub>L(&#981;)</sub> </td> </tr>\n";
		for i:=1 to #G`CharacterTable do
			str *:= "<tr><td>"*FixCharacterNamesHTML(G`CharacterNames[i])*"</td>";
			str *:= "<td>"*Sprint(res[h]`SimpleDims[i])*"</td>";
			str *:= "<td>"*Sprint(res[h]`VermaDecomposition[i][i])*"</td>";
			str *:= "<td>"*HTML(res[h]`SimplePSeries[i])*"</td>";
			str *:= "</tr>\n";
		end for;
		str *:= "</table>";

		str *:= "<br>\n";

	end if;

	if assigned res[h]`SimpleGModStruct then
		str *:= "<h3>Characters of the simple modules</h3>\n";

		str *:= "<table>\n";
		str *:= "<tr>\n";
		str *:= "<td> &#981; </td>\n";
		for i:=1 to #G`CharacterTable do
			str *:= "<td>L("*FixCharacterNamesHTML(G`CharacterNames[i])*")</td>\n";
		end for;
		for i:=1 to #G`CharacterTable do
			str *:= "<tr>";
			str *:= "<td>"*FixCharacterNamesHTML(G`CharacterNames[i])*"</td>\n";
			for j:=1 to #G`CharacterTable do
				str *:= "<td title=\"Multiplicity of "*FixCharacterNamesHTML(G`CharacterNames[i])*" in L("*FixCharacterNamesHTML(G`CharacterNames[j])*")\">"*Sprint(res[h]`SimpleGModStruct[j][i])*"</td>\n";
			end for;
			str *:= "</tr>";
		end for;
		str *:= "</table><br>\n";
	end if;

	if assigned res[h]`SimpleGradedGModStruct then
		str *:= "<h3>Graded characters of the simple modules</h3>\n";

		str *:= "<table>\n";
		str *:= "<tr>\n";
		str *:= "<td> &#981; </td>\n";
		for i:=1 to #G`CharacterTable do
			str *:= "<td>L("*FixCharacterNamesHTML(G`CharacterNames[i])*")</td>\n";
		end for;
		for i:=1 to #G`CharacterTable do
			str *:= "<tr>";
			str *:= "<td>"*FixCharacterNamesHTML(G`CharacterNames[i])*"</td>\n";
			for j:=1 to #G`CharacterTable do
				str *:= "<td title=\"Graded multiplicity of "*FixCharacterNamesHTML(G`CharacterNames[i])*" in L("*FixCharacterNamesHTML(G`CharacterNames[j])*")\">"*HTML(res[h]`SimpleGradedGModStruct[j][i])*"</td>\n";
			end for;
			str *:= "</tr>";
		end for;
		str *:= "</table><br>\n";
	end if;

	return str;

end function;

//============================================================================
intrinsic GordonDBAvailableData(W::GrpMat) -> Rec
{
	Analyzes the available data in a Gordon record.
}


	GordonDBDataAvailability := recformat<CMFamilies, CuspidalCMFamilies, SimpleDims, VermaDecomposition, SimplePSeries, SimpleGradedGModStruct, SimpleGModStruct, ExceptionalLocus>;

	avail := rec<GordonDBDataAvailability|>;
	avail`CMFamilies := "no";
	avail`ExceptionalLocus := "no";
	avail`CuspidalCMFamilies := "no";
	avail`SimpleDims := "no";
	avail`VermaDecomposition := "no";
	avail`SimplePSeries := "no";
	avail`SimpleGradedGModStruct := "no";
	avail`SimpleGModStruct := "no";


	if CHAMP_ExistsInDB(W`DBDir, "Cherednik/Gordon") then
		gordon := CHAMP_GetFromDB(W`DBDir, "Cherednik/Gordon");
	else
		return avail;
	end if;

	N := Rank(Codomain(CherednikParameter(W)));

	//CMFamilies
	for prop in Names(avail) do

		if N eq 1 then
			if prop in Names(gordon[1]) and assigned gordon[1]``prop then
				avail``prop := "all";
				if prop eq "CMFamilies" then
					avail`ExceptionalLocus := "yes";
				end if;
			end if;
		else
			speccount := 0;
			for f in Keys(gordon) diff {1} do
				if prop in Names(gordon[f]) and assigned gordon[f]``prop then
					speccount +:= 1;
				end if;
			end for;
			if prop in Names(gordon[1]) and assigned gordon[1]``prop then
				if Keys(gordon) eq {1} then
					avail``prop := "generic";
				else
					if speccount eq #(Keys(gordon) diff {1}) then
						avail``prop := "all";
						if prop eq "CMFamilies" then
							avail`ExceptionalLocus := "yes";
						end if;
					elif speccount gt 0 then
						avail``prop := "generic, some special";
					else
						avail``prop := "generic";
					end if;
				end if;
			elif speccount gt 0 then
				avail``prop := "some special";
			end if;
		end if;

	end for;

	return avail;

end intrinsic;

//============================================================================
intrinsic PrintGordonDBAvailableDataHTML()
{}

	file := CHAMP_GetDir()*"../Private/Results/Gordon/exceptional_summary.html";
	str := "<html>\n<head>\n";
	str *:= "<title>Summary of available data about the representation theory of the restricted rational Cherednik algebra for exceptional complex reflection groups";
	str *:="</title>\n";
	str *:= "<link href=\"style.css\" rel=\"stylesheet\" type=\"text/css\">\n";
	str *:= "</head>";
	str *:= "<body>\n";

	str *:= "<h1>Summary of available data about the representation theory of the restricted rational Cherednik algebra for exceptional complex reflection groups</h1>";

	str *:= "Computed by <a href=\"http://www.mathematik.uni-stuttgart.de/~thiel\">Ulrich Thiel</a>\n ";
	str *:= "using <a href=\"http://thielul.github.io/CHAMP\">CHAMP</a> (see LMS J. Comput. Math. 18 (2015), no. 1, 266&ndash;307). Some parts are joint work with <a href=\"http://www.math.univ-montp2.fr/~bonnafe/\">C&eacute;dric Bonnaf&eacute;</a>.<br>\n";
	str *:= "Last update on "*Date()*".\n<br><br>";

	str *:= "<table>";

	str *:= "<th>";
	str *:= "<td>Calogero&ndash;Moser<br>families</td>
		<td>Exceptional<br>locus</td>
		<td>Cuspidal<br>Calogero&ndash;Moser<br>families</td>
		<td>Dimensions of<br>simple modules</td>
		<td>Poincar&eacute; series of<br>simple modules</td>
		<td>W-module structure<br>of simple modules</td>
		<td>Graded W-module<br>structure of<br>simple modules</td>
		<td>Decomposition of<br>Verma modules</td>
		</th>\n";

	for n:=4 to 37 do
		W := ExceptionalComplexReflectionGroup(n);
		avail := GordonDBAvailableData(W);
		str *:= "<tr><td>"*Sprint(n)*"</td>";
		str *:= "<td>"*avail`CMFamilies*"</td>";
		str *:= "<td>"*avail`ExceptionalLocus*"</td>";
		str *:= "<td>"*avail`CuspidalCMFamilies*"</td>";
		str *:= "<td>"*avail`SimpleDims*"</td>";
		str *:= "<td>"*avail`SimplePSeries*"</td>";
		str *:= "<td>"*avail`SimpleGModStruct*"</td>";
		str *:= "<td>"*avail`SimpleGradedGModStruct*"</td>";
		str *:= "<td>"*avail`VermaDecomposition*"</td>";
		str *:= "</tr>";
	end for;

	str *:= "</table>";


	str *:= "</body>\n</html>";

	Write(file, str : Overwrite:=true);

end intrinsic;

//============================================================================
intrinsic PrintGordonHTMLForExceptionalGroup(n::RngIntElt)
{}

	G := ExceptionalComplexReflectionGroup(n);
	gordon := Gordon(G);
	N := Rank(Codomain(CherednikParameter(G)));

	CharacterTable(~G);
	file := CHAMP_GetDir()*"../Private/Results/Gordon/G"*Sprint(n)*".html";
	str := "<html>\n<head>\n";
	str *:= "<title>The representation theory of the restricted rational Cherednik algebra for G"*Sprint(n);
	if n eq 23 then
		str *:= "=H3";
	elif n eq 28 then
		str *:= "=F4";
	end if;
	str *:="</title>\n";
	str *:= "<link href=\"style.css\" rel=\"stylesheet\" type=\"text/css\">\n";
	str *:= "</head>";
	str *:= "<body>\n";
	str *:= "<h1>The representation theory of the restricted rational Cherednik algebra for G<sub>"*Sprint(n)*"</sub>";
	if n eq 23 then
		str *:= "=H<sub>3</sub>";
	elif n eq 28 then
		str *:= "=F<sub>4</sub>";
	end if;
	str *:= "</h1>\n";

	str *:= "Computed by <a href=\"http://www.mathematik.uni-stuttgart.de/~thiel\">Ulrich Thiel</a>\n ";
	str *:= "using <a href=\"http://thielul.github.io/CHAMP\">CHAMP</a> (see LMS J. Comput. Math. 18 (2015), no. 1, 266&ndash;307). Some parts are joint work with <a href=\"http://www.math.univ-montp2.fr/~bonnafe/\">C&eacute;dric Bonnaf&eacute;</a>.<br>\n";
	str *:= "Last update on "*Date()*".\n<br><br>";

	str *:= "<h2>Summary</h2>";

	//availability
	str *:= "<h3>Available data</h3>";
	str *:= "<table>\n";

	avail := GordonDBAvailableData(G);

	str *:= "<tr><td>Calogero&ndash;Moser families</td><td>"*Sprint(avail`CMFamilies)*"</td></tr>\n";
	str *:= "<tr><td>Exceptional locus</td><td>"*Sprint(avail`ExceptionalLocus)*"</td></tr>\n";
	str *:= "<tr><td>Cuspidal Calogero&ndash;Moser families</td><td>"*Sprint(avail`CuspidalCMFamilies)*"</td></tr>\n";
	str *:= "<tr><td>Dimensions of simple modules</td><td>"*Sprint(avail`SimpleDims)*"</td></tr>\n";
	str *:= "<tr><td>Poincar&eacute; series of simple modules</td><td>"*Sprint(avail`SimplePSeries)*"</td></tr>\n";
	str *:= "<tr><td>W-module structure of simple modules</td><td>"*Sprint(avail`SimpleGModStruct)*"</td></tr>\n";
	str *:= "<tr><td>Graded W-module structure of simple modules</td><td>"*Sprint(avail`SimpleGradedGModStruct)*"</td></tr>\n";
	str *:= "<tr><td>Decomposition of Verma modules</td><td>"*Sprint(avail`VermaDecomposition)*"</td></tr>\n";


	str *:= "</table>\n";

	//results
	str *:= "<h3>Results</h3>";

	rou := RouquierFamilies(G);

	if assigned gordon[1]`CMFamilies then
		str *:= "<ul type=\"disc\" style=\"margin-top:0px; margin-bottom:20px;\">\n";

		//cm families
		if rou[1] eq gordon[1]`CMFamilies then
			str *:= "<li>The generic Rouquier families equal the generic Calogero&ndash;Moser families&mdash;Martino's generic parameter conjecture holds.</li>";
		else
			roufams := rou[1];
			cmfams := gordon[1]`CMFamilies;
			roufinercm := true;
			for cmfam in cmfams do
				subroufams := { F : F in roufams | F subset cmfam };
				if not { i : i in F, F in subroufams} eq cmfam then
					roufinercm := false;
				end if;
			end for;
			if roufinercm then
				str *:= "<li>The generic Rouquier families refine, but do <b>not</b> equal, the generic the Calogero&ndash;Moser families&mdash;Martino's generic parameter conjecture does <b>not</b> hold.</li>";
			else
				str *:= "<li>The generic Rouquier families do <b>not</b> refine the generic Calogero&ndash;Moser families&mdash;Martino's generic parameter conjecture does <b>not</b> hold.</li>";
			end if;
		end if;

		allcms := false;
		if N gt 1 and #Keys(gordon) gt 1 then
			allcms := true;
		end if;

		if allcms then
			P := Universe(Keys(gordon));
			Q := Universe(Keys(rou));
			roufinercm := true;
			roueqcm := true;
			rouexeqcmgen := true;
			rouexcontainedincmgen := true;

			rouex := { P!NormalizeRationalHyperplaneEquation(P!G`MartinoSharp(H)) : H in Keys(rou) };
			if rouex ne Keys(gordon) then
				rouexeqcmgen := false;
			end if;
			if not rouex subset Keys(gordon) then
				rouexcontainedincmgen := false;
			end if;
			for H in Keys(rou) do
				roufams := rou[H];
				Hsharp := P!NormalizeRationalHyperplaneEquation(P!G`MartinoSharp(H));
				if Hsharp notin Keys(gordon) then
					cmfams := gordon[1]`CMFamilies;
					rouexcontainedincmgen := false;
					rouexeqcmgen := false;
				else
					cmfams := gordon[Hsharp]`CMFamilies;
				end if;
				if roufams ne cmfams then
					roueqcm := false;
					//check if each cm family is union of rouquier families (martino's conjecture)
					for cmfam in cmfams do
						subroufams := { F : F in roufams | F subset cmfam };
						if not { i : i in F, F in subroufams} eq cmfam then
							roufinercm := false;
						end if;
					end for;
				end if;
			end for;
			for H in Keys(gordon) do
				cmfams := gordon[H]`CMFamilies;
				Hsharp := NormalizeRationalHyperplaneEquation(Q!G`MartinoSharp(H));
				if Hsharp notin Keys(rou) then
					roufams := rou[1];
				else
					roufams := rou[Hsharp];
				end if;
				if roufams ne cmfams then
					roueqcm := false;
					//check if each cm family is union of rouquier families (martino's conjecture)
					for cmfam in cmfams do
						subroufams := { F : F in roufams | F subset cmfam };
						if not { i : i in F, F in subroufams} eq cmfam then
							roufinercm := false;
						end if;
					end for;
				end if;
			end for;


			if roueqcm then
				str *:= "<li>The Rouquier k<sup>#</sup>-families equal the Calogero&ndash;Moser k-families for all k&mdash;Martino's special parameter conjecture holds.</li>";
			else
				if roufinercm then
					str *:= "<li>The Rouquier k<sup>#</sup>-families refine, but in general do not equal, the Calogero&ndash;Moser k-families for all k&mdash;Martino's special parameter conjecture holds.</li>";
				else
					str *:= "<li>The Rouquier k<sup>#</sup>-families do <b>not</b> refine the Calogero&ndash;Moser k-families for some k&mdash;Martino's special parameter conjecture does <b>not</b> hold.</li>";
				end if;
			end if;

			if rouexeqcmgen then
				str *:= "<li>The image of the locus of essential hyperplanes for Rouquier families under # equals the exceptional locus for Calogero&ndash;Moser families.</li>";
			else
				if rouexcontainedincmgen then
					str *:= "<li>The image of the locus of essential hyperplanes for Rouquier families under # is contained in, but not equal to, the exceptional locus for Calogero&ndash;Moser families.</li>";
				else
					str *:= "<li>The image of the locus of essential hyperplanes for Rouquier families under # is <b> not </b>contained in the exceptional locus for Calogero&ndash;Moser families.</li>";
				end if;
			end if;

		end if;

	end if;

	str *:= "</ul>";

	//now, the data

	//str *:= "Note: In the larger tables each cell has a mouseover tooltip providing information about the cell.<br><br>\n";

	str *:= "Quick navigation: ";
	str *:= "<a href=\"#hyperplanes\">Exceptional hyperplanes</a><br><br>";

	str *:= "<h2><a href=\"#generic\" name=\"generic\">For generic parameters</a></h2>\n";

	str *:= PrintGordonHTML(G, gordon, 1);

	hyperplanes := Sort(SetToSequence(Keys(gordon) diff {1}));

	str *:= "<h2><a name=\"hyperplanes\" href=\"#hyperplanes\">Exceptional hyperplanes</a></h2>\n";

	if #hyperplanes eq 0 then
		if Ngens(Universe(Keys(gordon))) ne 1 then
			str *:= "Unknown<br>";
		else
			str *:= "There are none.<br>";
		end if;
	end if;

	for i:=1 to #hyperplanes do
		H := Sprint(hyperplanes[i]);
		H := Replace(H, "\\*", "");
		H := Replace(H, "_{", "<sub>");
		H := Replace(H, "}", "</sub>");
		H := Replace(H, "-", "&minus;");
		str *:= "<a href=\"#hyerplane-"*Sprint(i)*"\">"*H*"</a><br> \n";
	end for;

	str *:= "<br>\n";

	for i:=1 to #hyperplanes do
		h := hyperplanes[i];
		H := Sprint(hyperplanes[i]);
		H := Replace(H, "\\*", "");
		H := Replace(H, "_{", "<sub>");
		H := Replace(H, "}", "</sub>");
		H := Replace(H, "-", "&minus;");
		str *:= "<h2><a name=\"hyerplane-"*Sprint(i)*"\" href=\"#hyerplane-"*Sprint(i)*"\"> For the generic point of the hyperplane "*H*"</a></h2>\n";
		str *:= "Quick navigation: ";
		str *:= "<a href=\"#hyperplanes\">Exceptional hyperplanes</a>, ";
		str *:= "<a href=\"#generic\">For generic parameters</a><br><br>";
		str *:= PrintGordonHTML(G, gordon, h);
	end for;

	str *:= "</body>\n</html>";

	Write(file, str : Overwrite:=true);

end intrinsic;
