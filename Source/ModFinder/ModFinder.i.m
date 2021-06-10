/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Abstract subspaces and finding modules with prescribed abstract structure.

  History:
      Friday, September 20, 2013 11:49:14 : Initial.
*/

declare type SubAbs;

declare attributes SubAbs:
    BasespaceDimension,
    SubspaceDimension,
    Variables,
    RowVariablesMap,
		VariableRowsMap,
		VariableWeights,
		VariablesSortedByWeight,
		ColumnVariablesMap,
    KnownPart,
    UnknownPart
    ;

declare type ModFinder;

declare attributes ModFinder:
    AbstractStructure,
    RowDegrees,
		RowDegreeMap,
		DegreeRowsMap,
		DeterminedVariables,
		DeterminedVariablesMap,
		NonDetRows,
		SupportedColsPerDegree,
		NumberOfEquations,
		NumberOfAuxVars,
    Rounds,
    BaseModule
    ;

declare verbose ModFinder, 5;

vardatarec := recformat< dependencyblock, nondetdepvars, depvars, deprows, auxvars, totalgoodcolsinaux, totalgoodcolsinker, kerdata, auxdata>;

//====================================================================
intrinsic AbstractSubspace(U::SeqEnum, N::RngIntElt) -> SubAbs
/*
    History:
        Friday, September 20, 2013 11:54:50: Initial.
*/
{Construct the abstract structure of the subspace with basis U.}

    M := New(SubAbs);
    M`BasespaceDimension := N;
    M`SubspaceDimension := #U;

    M`KnownPart := SparseMatrix(M`SubspaceDimension, M`BasespaceDimension);
	M`UnknownPart := SparseMatrix(M`SubspaceDimension, M`BasespaceDimension);

    entries := [];
	for i:=1 to M`SubspaceDimension do
		b := U[i];
		supp := Sort(SetToSequence(Support(b))); //order to be safe
		M`KnownPart[i][supp[1]] := 1; //because of echelon form of U
		for j:=2 to #supp do
			k := Position(entries, b[supp[j]]);
			if k eq 0 then
				Append(~entries, b[supp[j]]);
				k := #entries;
			end if;
			M`UnknownPart[i][supp[j]] := k;
		end for;
	end for;


    M`Variables := {@ i : i in [1..#entries] @};

	rowvariables := [* *];
	for i:=1 to M`SubspaceDimension do
		Append(~rowvariables, {@ M`UnknownPart[i][j] : j in Support(M`UnknownPart,i) @});
	end for;
	M`RowVariablesMap := rowvariables;

	variablerows := [];
	for k:=1 to #M`Variables do
		Append(~variablerows, {@ i : i in [1..M`SubspaceDimension] | k in M`RowVariablesMap[i] @});
	end for;
	M`VariableRowsMap := variablerows;

	columnvariables := [* *];
	for j:=1 to M`BasespaceDimension do
		Append(~columnvariables, {@ M`UnknownPart[i][j] : i in Support(Column(M`UnknownPart,j)) @});
	end for;
	M`ColumnVariablesMap := columnvariables;

	M`VariablesSortedByWeight := {@ {@ i : i in [1..#entries] @} @};

    return M;

end intrinsic;

//====================================================================
intrinsic Print(M::SubAbs)
/*
    History:
        Friday, September 20, 2013 12:17:04: Initial.
*/
{}

    printf "Abstract subspace of dimension %o with %o variables in vector space of dimension %o.", M`SubspaceDimension, #M`Variables, M`BasespaceDimension;

end intrinsic;

//====================================================================
intrinsic AbstractSubspace(U::ModRng, V::ModRng) -> SubAbs
/*
    History:
        Friday, September 20, 2013 12:12:59: Initial.
*/
{}

    emb := Morphism(U,V);
	return AbstractSubspace([ emb(U.i) : i in [1..Dimension(U)] ], Dimension(V));

end intrinsic;

//====================================================================
intrinsic InitializeModuleFinder(A::SubAbs, V::ModGrOld) -> ModFinder
/*
    History:
        Saturday, September 21, 2013 20:08:35: Initial.
*/
{}

    F := New(ModFinder);
    F`AbstractStructure := A;
    SubAbsHelper_CheckIfGradedAndGetRowDegreeMap(V, ~F);
	K := BaseRing(V);

	F`DeterminedVariables := {@ @};
	F`DeterminedVariablesMap := [* false : i in [1..#F`AbstractStructure`Variables] *];

	F`NonDetRows := {@ c : c in [1..F`AbstractStructure`SubspaceDimension] @};
	SubAbsHelper_NonDetRows(~F);
	SubAbsHelper_SupportedColsPerDegree(~F);
	F`NumberOfEquations := 0;
	F`NumberOfAuxVars := 0;

    //new
    F`Rounds := 0;
    F`BaseModule := V;

	return F;

end intrinsic;

//====================================================================
intrinsic Print(F::ModFinder)
/*
    History:
        Saturday, September 21, 2013 20:26:46: Initial.
*/
{}

    printf "Finder for module with prescribed abstract structure.\n";
    printf "Current status: %o out of %o variables determined (%o%%); %o rounds of finder have been run; %o equations have been used.", #F`DeterminedVariables, #F`AbstractStructure`Variables, ChangePrecision((100.0*#F`DeterminedVariables)/#F`AbstractStructure`Variables, 4), F`Rounds, F`NumberOfEquations;

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_CheckIfGradedAndGetRowDegreeMap(V::ModGrOld, ~F::ModFinder)
/*
    History:
        Saturday, September 21, 2013 20:14:17: Initial; ported from old project.
*/
{}
	F`RowDegreeMap := [];

	for i:=1 to F`AbstractStructure`SubspaceDimension do
		supp := Support(F`AbstractStructure`KnownPart[i]) join Support(F`AbstractStructure`UnknownPart[i]);
		degrees := { Degree(V.j) : j in supp };
		if #degrees gt 2 then
			error "Abstract submodule won't be graded inside this module.";
		end if;
		Append(~F`RowDegreeMap, SetToSequence(degrees)[1]);
	end for;

	F`RowDegrees := {@ F`RowDegreeMap[c] : c in [1..F`AbstractStructure`SubspaceDimension] @};

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_NonDetRows(~F::ModFinder)
/*
	Update set of all rows in which a non-determined variable occurs.

    History:
        Saturday, September 21, 2013 20:18:27: Initial; ported from old project.
*/
{}

	F`NonDetRows := {@ row : row in F`NonDetRows | F`AbstractStructure`RowVariablesMap[row] subset F`DeterminedVariables eq false @};

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetVariableDependencies(F::ModFinder, varnum::RngIntElt) -> SetIndx, SetIndx
/*
	Return all rows which contain variable varnum and all variables occuring in
	these rows.

    History:
        Saturday, September 21, 2013 20:38:59: Initial; ported from old project.
*/
{}

	deprows := {@ r : r in [1..F`AbstractStructure`SubspaceDimension] | varnum in F`AbstractStructure`RowVariablesMap[r] and F`AbstractStructure`RowVariablesMap[r] subset F`DeterminedVariables eq false @};
	depvars := FlatFixed({@ F`AbstractStructure`RowVariablesMap[r] diff F`DeterminedVariables : r in deprows @});
	return depvars, deprows;

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetRowNonDetSupp(F::ModFinder, row::RngIntElt) -> SetIndx
/*
	Get all columns in row which contain a non-determined variable.

    History:
        Saturday, September 21, 2013 20:39:06: Initial; ported from old project.
*/
{}

	return {@ col : col in Support(F`AbstractStructure`UnknownPart[row]) | F`AbstractStructure`UnknownPart[row][col] notin F`DeterminedVariables @};

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetGoodSupportedColsForImageDegree(F::ModFinder, imagedeg::RngIntElt) -> SetIndx
/*
	Get all supported cols of imagedeg where all entries are determined.

    History:
        Saturday, September 21, 2013 20:39:13: Initial; ported from old project.
*/
{}

	goodcols := F`SupportedColsPerDegree[Position(F`RowDegrees, imagedeg)];

	for col in goodcols do
		colvars := F`AbstractStructure`ColumnVariablesMap[col];

		if colvars subset F`DeterminedVariables eq false then
			goodcols diff:={@ col @};
		end if;
	end for;

	return goodcols;

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetNonSupportedColsForImageDegree(F::ModFinder,imagedeg::RngIntElt) -> SetIndx
/*
	All columns of imagedeg which are non-supported.

    History:
        Saturday, September 21, 2013 20:39:20: Initial; ported from old project.
*/
{}

	if imagedeg notin F`RowDegrees then
		goodcols := {@ c : c in [1..F`AbstractStructure`BasespaceDimension] @};
	else
		goodcols := {@ c : c in [1..F`AbstractStructure`BasespaceDimension] | c notin F`SupportedColsPerDegree[Position(F`RowDegrees, imagedeg)] @};
	end if;

	return goodcols;

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetGoodSupportedColsForRowGenPair(F::ModFinder, row::RngIntElt, gendeg::RngIntElt, genmat::MtrxSprs) -> SetIndx
/*
	All GoodSupportedColsForImageDegree for imagedeg equal to rowdegree+gendegree
	such that genmat has support in these cols.

    History:
        Saturday, September 21, 2013 20:39:26: Initial; ported from old project.
*/
{}

	return SubAbsHelper_GetGoodSupportedColsForImageDegree(F,F`RowDegreeMap[row]+gendeg) meet FlatFixed( {@ SetToIndexedSet(Support(genmat[c])) : c in SubAbsHelper_GetRowNonDetSupp(F,row) @} );

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetGoodSupportedColsForRowGenPair(F::ModFinder, row::RngIntElt, gendeg::RngIntElt, genmat::Mtrx) -> SetIndx
/*
	All GoodSupportedColsForImageDegree for imagedeg equal to rowdegree+gendegree
	such that genmat has support in these cols.

    History:
        Saturday, September 21, 2013 20:39:26: Initial; ported from old project.
*/
{}

	return SubAbsHelper_GetGoodSupportedColsForImageDegree(F,F`RowDegreeMap[row]+gendeg) meet FlatFixed( {@ SetToIndexedSet(Support(genmat[c])) : c in SubAbsHelper_GetRowNonDetSupp(F,row) @} );

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetGoodNonSupportedColsForRowGenPair(F::ModFinder, row::RngIntElt, gendeg::RngIntElt, genmat::MtrxSprs) -> SetIndx
/*
	All GoodNonSupportedColsForImageDegree for imagedeg equal to rowdegree+gendegree
	such that genmat has support in these cols.

    History:
        Saturday, September 21, 2013 20:40:10: Initial; ported from old project.
*/
{}

	return SubAbsHelper_GetNonSupportedColsForImageDegree(F,F`RowDegreeMap[row]+gendeg) meet FlatFixed( {@ SetToIndexedSet(Support(genmat[c])) : c in SubAbsHelper_GetRowNonDetSupp(F,row) @} );

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetGoodNonSupportedColsForRowGenPair(F::ModFinder, row::RngIntElt, gendeg::RngIntElt, genmat::Mtrx) -> SetIndx
/*
	All GoodNonSupportedColsForImageDegree for imagedeg equal to rowdegree+gendegree
	such that genmat has support in these cols.

    History:
        Saturday, September 21, 2013 20:40:10: Initial; ported from old project.
*/
{}

	return SubAbsHelper_GetNonSupportedColsForImageDegree(F,F`RowDegreeMap[row]+gendeg) meet FlatFixed( {@ SetToIndexedSet(Support(genmat[c])) : c in SubAbsHelper_GetRowNonDetSupp(F,row) @} );

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetGoodImageRowsForCols(F::ModFinder, imagedeg::RngIntElt, goodcols::SetIndx) -> SetIndx
/*
	All rows of imagedeg such that known part or unknown part in this row have goodcols.

    History:
        Saturday, September 21, 2013 20:41:13: Initial; ported from old project.
*/
{}

	return {@ row : row in [1..F`AbstractStructure`SubspaceDimension] | F`RowDegreeMap[row] eq imagedeg and #( (Support(F`AbstractStructure`KnownPart[row]) join Support(F`AbstractStructure`UnknownPart[row])) meet goodcols) ne 0 @};

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_SupportedColsPerDegree(~F::ModFinder)
/*
	For each degree all columns in which known part or unknown part are supported.

    History:
        Saturday, September 21, 2013 20:38:52: Initial; ported from old project.
*/
{}

	F`SupportedColsPerDegree := {@ @};
	for deg in F`RowDegrees do
		rowsofdeg := {@ r : r in [1..F`AbstractStructure`SubspaceDimension] | F`RowDegreeMap[r] eq deg @};
		supprows := FlatFixed({@ SetToIndexedSet(Support(F`AbstractStructure`KnownPart[r])) join SetToIndexedSet(Support(F`AbstractStructure`UnknownPart[r])) : r in rowsofdeg @});
		F`SupportedColsPerDegree join:= {@ supprows @};
	end for;

end intrinsic;


//====================================================================
intrinsic SubAbsHelper_GetVariableData(F::ModFinder, levelgens::SetIndx, var::RngIntElt) -> Rec
/*
    History:
        Wednesday, January 15, 2014 at 16:52:17: Changes Exclude to Diff.
        Tuesday, October 08, 2013 19:26:18: Initial.
*/
{}

    depvars, deprows := SubAbsHelper_GetVariableDependencies(F,var);

    nondetdepvars := Diff(depvars, IndexedSetToSet(F`DeterminedVariables));

    kerdata := {@ @};
    for row in deprows do
        for gen in levelgens do
            goodcols := SubAbsHelper_GetGoodNonSupportedColsForRowGenPair(F,row,F`BaseModule`MatrixDegrees[gen],F`BaseModule`Matrices[gen]);

            if #goodcols ne 0 then
                kerdata join:={@ <row,gen,goodcols> @};
            end if;
        end for;
    end for;

    auxdata := {@ @};
    for row in deprows do
        for gen in levelgens do
            imagedeg := F`RowDegreeMap[row]+F`BaseModule`MatrixDegrees[gen];
            if imagedeg notin F`RowDegrees then
                continue;
            end if;
            goodcols := SubAbsHelper_GetGoodSupportedColsForRowGenPair(F,row,F`BaseModule`MatrixDegrees[gen],F`BaseModule`Matrices[gen]);

            if #goodcols ne 0 then
                goodimagerows := SubAbsHelper_GetGoodImageRowsForCols(F, imagedeg, goodcols);

                if #goodimagerows ne 0 then
                    auxdata join:={@ <row,gen,goodcols,goodimagerows> @};
                end if;
            end if;
        end for;
    end for;

    auxvars := ArraySum([ #data[4] : data in auxdata ]);
    totalgoodcolsinaux := ArraySum([ #data[3] : data in auxdata]);
    totalgoodcolsinker:= ArraySum([ #data[3] : data in kerdata]);

    dependencyblock := SubAbsHelper_GetVariableDependenciesComplete(F,var);

    return rec<vardatarec | dependencyblock:=dependencyblock, nondetdepvars:=nondetdepvars, depvars:=depvars, deprows:=deprows, auxvars:=auxvars, totalgoodcolsinaux:=totalgoodcolsinaux, totalgoodcolsinker:=totalgoodcolsinker, kerdata:=kerdata, auxdata:=auxdata>;
    //LHS of system will be of size (#nondetdepvars+auxvars) times (totalgoodcolsinaux+totalgoodcolsinker)

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetVariableDependenciesComplete(F::ModFinder, var::RngIntElt) -> SetIndx
/*
    History:
        Wednesday, October 09, 2013 14:29:41: Initial.
*/
{}
    depvarstotal, deprows := SubAbsHelper_GetVariableDependencies(F,var);
    depvarstotal join:={@var@};
    checked := {@var@};

    while checked ne depvarstotal and not IsEmpty(depvarstotal) do
        var := (depvarstotal diff checked)[1];
        checked join:={@var@};
        depvars, deprows := SubAbsHelper_GetVariableDependencies(F,var);
        depvars join:={@var@};
        depvarstotal join:=depvars;
    end while;

    return depvarstotal;

end intrinsic;

//====================================================================
intrinsic SubAbsHelper_GetSystem(F::ModFinder, vardata::Rec : InvolveAuxialiary:=true) -> Mtrx, Mtrx
/*
    History:
        Saturday, October 12, 2013 13:43:11: Initial.
*/
{}

    K := BaseRing(F`BaseModule);

    //shortcut for data of workdatum
    depvars := vardata`depvars;
    nondetdepvars := vardata`nondetdepvars;
    totalgoodcolsinker := vardata`totalgoodcolsinker;
    totalgoodcolsinaux := vardata`totalgoodcolsinaux;
    auxvars := vardata`auxvars;
    kerdata := vardata`kerdata;
    auxdata := vardata`auxdata;
    deprows := vardata`deprows;

    if InvolveAuxialiary then
        lhs := ZeroMatrix(K, #nondetdepvars+auxvars, totalgoodcolsinaux+totalgoodcolsinker);
        rhs := Zero(VectorSpace(K, totalgoodcolsinaux+totalgoodcolsinker));
    else
        lhs := ZeroMatrix(K, #nondetdepvars, totalgoodcolsinker);
        rhs := Zero(VectorSpace(K, totalgoodcolsinker));
    end if;

    eqnnum := 0;
    for i:=1 to #kerdata do
        data := kerdata[i];
        row := data[1];
        gen := data[2];
        goodcols := data[3];
        rowknownsupp := SetToIndexedSet(Support(F`AbstractStructure`KnownPart[row]));
        rowunknownsupp := SetToIndexedSet(Support(F`AbstractStructure`UnknownPart[row]));
        rownondetsupp := SubAbsHelper_GetRowNonDetSupp(F,row);
        for j:=1 to #goodcols do
            eqnnum +:= 1;
            goodcol := goodcols[j];

        //RHS
            for c in rowknownsupp do
                rhs[eqnnum] -:= F`AbstractStructure`KnownPart[row][c]*F`BaseModule`Matrices[gen][c][goodcol];
            end for;

            rowunknownbutdetsupp := rowunknownsupp diff rownondetsupp;
            for c in rowunknownbutdetsupp do
                rhs[eqnnum] -:= F`DeterminedVariablesMap[F`AbstractStructure`UnknownPart[row][c]]*F`BaseModule`Matrices[gen][c][goodcol];
            end for;

        //LHS
            for c in rownondetsupp do
                varnum := F`AbstractStructure`UnknownPart[row][c];
                curvarnum := Position(nondetdepvars, varnum);
                lhs[curvarnum][eqnnum] +:= F`BaseModule`Matrices[gen][c][goodcol];
            end for;

        end for;

    end for;

    if InvolveAuxialiary then
        auxvarblockcounter := 0;
        for i:=1 to #auxdata do
            data := auxdata[i];
            row := data[1];
            gen := data[2];
            goodcols := data[3];
            goodrows := data[4];
            rowknownsupp := SetToIndexedSet(Support(F`AbstractStructure`KnownPart[row]));
            rowunknownsupp := SetToIndexedSet(Support(F`AbstractStructure`UnknownPart[row]));
            rownondetsupp := SubAbsHelper_GetRowNonDetSupp(F,row);
            for j:=1 to #goodcols do
                eqnnum +:= 1;
                goodcol := goodcols[j];

            //RHS
                for c in rowknownsupp do
                    rhs[eqnnum] -:= F`AbstractStructure`KnownPart[row][c]*F`BaseModule`Matrices[gen][c][goodcol];
                end for;

                rowunknownbutdetsupp := rowunknownsupp diff rownondetsupp;
                for c in rowunknownbutdetsupp do
                    rhs[eqnnum] -:= F`DeterminedVariablesMap[F`AbstractStructure`UnknownPart[row][c]]*F`BaseModule`Matrices[gen][c][goodcol];
                end for;

            //LHS
                for c in rownondetsupp do
                    varnum := F`AbstractStructure`UnknownPart[row][c];
                    curvarnum := Position(nondetdepvars, varnum);
                    lhs[curvarnum][eqnnum] +:= F`BaseModule`Matrices[gen][c][goodcol];
                end for;

                for k:=1 to #goodrows do
                    goodrow := goodrows[k];
                    //print [#nondetdepvars,auxvarblockcounter,k];
                    if goodcol in Support(F`AbstractStructure`UnknownPart[goodrow]) then
                        varnum := F`AbstractStructure`UnknownPart[goodrow][goodcol];
                        lhs[#nondetdepvars+auxvarblockcounter+k][eqnnum] -:= F`DeterminedVariablesMap[varnum];
                    else
                        lhs[#nondetdepvars+auxvarblockcounter+k][eqnnum] -:= F`AbstractStructure`KnownPart[goodrow][goodcol];
                    end if;
                end for;

            end for;

            auxvarblockcounter +:= #goodrows;

        end for;
    end if;

    return lhs,rhs;

end intrinsic;

//====================================================================
intrinsic FindSubmoduleWithAbstractStructure(~F::ModFinder :
    gens:={},
    solver:="Magma",
    sortdatums:=true,
    removeredundant:=false,
    updateinterval:="Automatic",
    saveuseless:=true
)
/*
    Here it is, the finder for modules with prescribed abstract structure.

    History:
        Saturday, September 21, 2013 20:47:03: Initial; ported from old project.

        Tuesday, October 08, 2013 20:56:24: New system selection.

        Wednesday, October 09, 2013 17:50:04: Again, new selection.
*/
{}

	if F`AbstractStructure`SubspaceDimension eq 0 then
		return;
	end if;

    //INIT DATA
	K := BaseRing(F`BaseModule);
	maxrowdegree := Maximum(F`RowDegrees);
	minrowdegree := Minimum(F`RowDegrees);
	SubAbsHelper_NonDetRows(~F);
	maxauxvars := 0;
	maxeqnnum := 0;
	eqnsizes := [];
    uselesssystems := {};

	if #F`DeterminedVariables eq #F`AbstractStructure`Variables then
		print "All variables are already determined.";
		return;
	end if;

	vprint ModFinder, 3: "Abstract submodule has "*Sprint(#F`AbstractStructure`Variables)*" variables.";

	if Type(updateinterval) eq MonStgElt and updateinterval eq "Automatic" then
        updateinterval := Round(0.20*#F`AbstractStructure`Variables);
    end if;

    vprint ModFinder, 3: "Update interval "*Sprint(updateinterval)*".";

	if gens eq {} then
        gens := {i : i in [1..Ngens(F`BaseModule)]};
        vprint ModFinder, 3: "Selected generators: All.";
    else
        vprint ModFinder, 3: "Selected generators: "*Sprint(gens)*".";
    end if;

    //collect all possible datums and data
    vprint ModFinder,5: "Collecting datums and data.";
    counter:=0;
    datums := { <var,SetToIndexedSet(gens)> : var in FlatFixed(F`AbstractStructure`VariablesSortedByWeight) | var notin F`DeterminedVariables };
    vardata := AssociativeArray(Universe({<1,{1}>}));
    for x in datums do
        counter+:=1;
        PrintPercentage(counter, #datums);
        vardata[x] := SubAbsHelper_GetVariableData(F,SetToIndexedSet(x[2]),x[1]);
    end for;

     //remove reduntant variables
    counter:=0;
    if removeredundant then
        vprint ModFinder,5: "Removing redundant datums.";
        nonredundantdatums := {};
        for x in datums do
            counter+:=1;
            PrintPercentage(counter, #datums);
            if exists{y : y in nonredundantdatums | x[2] eq y[2] and IsEqual(vardata[y], vardata[x])} then
                continue;
            else
                nonredundantdatums join:={x};
            end if;
        end for;
    else
        nonredundantdatums := datums;
    end if;

    trieddatums := {};
    nondetvars := { x[1] : x in datums };
    newdetvars := {};
    systemscounter := 0;
    usefulsystemscounter := 0;
    InvolveAuxialiary := false;
    progress := true;
    update := false;

    //Now, the fun starts
	zeit:=Cputime();
	while progress and trieddatums ne nonredundantdatums and #F`DeterminedVariables ne #F`AbstractStructure`Variables do

        //collect variable data
        if update then
            vprint ModFinder,5: "Updating datums.";
            counter:=0;
            for x in nonredundantdatums do
                counter+:=1;
                PrintPercentage(counter, #nonredundantdatums);
                updated := false;
                if IsDefined(vardata, x) then
                    oldvardata := vardata[x];
                    updated := true;
                end if;
                vardata[x] := SubAbsHelper_GetVariableData(F,x[2],x[1]);
                if updated then
                    if not IsEqual(oldvardata, vardata[x]) then
                        trieddatums diff:={x};
                    end if;
                end if;
            end for;
        end if;
        update := false;

        //SORTING PROCESS
        if sortdatums then
            //collect those variables which _probably_ cannot be determined as there do not have sufficiently many equations

            //bad datums without auxiliary
            if not InvolveAuxialiary then
                baddatums := { x : x in nonredundantdatums | ((vardata[x]`totalgoodcolsinker) eq 0 or (#vardata[x]`nondetdepvars)/(vardata[x]`totalgoodcolsinker) gt 1) };
            else
            //bad datums with auxiliary
                baddatums := { x : x in nonredundantdatums | ((vardata[x]`totalgoodcolsinaux+vardata[x]`totalgoodcolsinker) eq 0 or (#vardata[x]`nondetdepvars+vardata[x]`auxvars)/(vardata[x]`totalgoodcolsinaux+vardata[x]`totalgoodcolsinker) gt 1) };
            end if;

            gooddatums := nonredundantdatums diff baddatums;

            //sort good datums
            selectdata := [];
            for x in gooddatums do
                lhs,rhs := SubAbsHelper_GetSystem(F, vardata[x]:InvolveAuxialiary:=InvolveAuxialiary);
                val := <NumberOfNonZeroEntries(lhs), Ncols(lhs), Nrows(lhs), &+[NumberOfNonZeroEntries(F`BaseModule`Matrices[i]) : i in x[2]], x[1], SetToSequence(x[2])>;
                Append(~selectdata, val);
            end for;

            Sort(~selectdata);

            selection := [ < selectdata[i][#selectdata[i]-1], SequenceToIndexedSet(selectdata[i][#selectdata[i]])> : i in [1..#selectdata] ];

            //sort bad datums
            selectdata := [];
            for x in baddatums do
                if not InvolveAuxialiary then
                    val := <vardata[x]`totalgoodcolsinker, #vardata[x]`nondetdepvars, &+[NumberOfNonZeroEntries(F`BaseModule`Matrices[i]) : i in x[2]], x[1], SetToSequence(x[2])>;
                else
                    val := <vardata[x]`totalgoodcolsinaux+vardata[x]`totalgoodcolsinker, #vardata[x]`nondetdepvars+vardata[x]`auxvars, &+[NumberOfNonZeroEntries(F`BaseModule`Matrices[i]) : i in x[2]],x[1],SetToSequence(x[2])>;
                end if;
                Append(~selectdata,val);
            end for;

            Sort(~selectdata);

            selection := selection cat [ < selectdata[i][#selectdata[i]-1], SequenceToIndexedSet(selectdata[i][#selectdata[i]])> : i in [1..#selectdata] ];

            vprint ModFinder, 5: "There are "*Sprint(#gooddatums)*" good datums and "*Sprint(#baddatums)*" bad datums.";
        else
            selection := [ x : x in nonredundantdatums ];
        end if;

        //NOW WE START SOLVING THE SYSTEMS IN THE DETERMINED SEQUENCE
        progress := false;
        //print selection;
        while not IsEmpty(selection) and #F`DeterminedVariables ne #F`AbstractStructure`Variables  do

            workdatum := < selection[1][#selection[1]-1], selection[1][#selection[1]]> ;

            vprint ModFinder, 5: "        Trying datum "*Sprint(workdatum)*".";

            vardata[workdatum] := SubAbsHelper_GetVariableData(F,workdatum[2],workdatum[1]);

            Remove(~selection, 1);
            trieddatums join:={workdatum};

            //shortcut for data of workdatum
            depvars := vardata[workdatum]`depvars;
            nondetdepvars := vardata[workdatum]`nondetdepvars;
            totalgoodcolsinker := vardata[workdatum]`totalgoodcolsinker;
            totalgoodcolsinaux := vardata[workdatum]`totalgoodcolsinaux;
            auxvars := vardata[workdatum]`auxvars;
            kerdata := vardata[workdatum]`kerdata;
            auxdata := vardata[workdatum]`auxdata;
            deprows := vardata[workdatum]`deprows;

            if not InvolveAuxialiary and totalgoodcolsinker eq 0 then
                vprint ModFinder, 5: "            System is useless.";
                continue;
            elif InvolveAuxialiary and totalgoodcolsinaux+totalgoodcolsinker eq 0 then
                vprint ModFinder, 5: "            System is useless.";
                continue;
            end if;

            //vprint ModFinder, 5: "            Setting up system.";
            lhs, rhs := SubAbsHelper_GetSystem(F, vardata[workdatum]:InvolveAuxialiary:=InvolveAuxialiary);


            if MD5(lhs) in uselesssystems then
                vprint ModFinder, 5: "            System is useless.";
                continue;
            end if;

            systemscounter +:= 1;

            vprint ModFinder, 5: "            Variable has "*Sprint(#nondetdepvars)*" dependent non-determined variables: "*Sprint(nondetdepvars)*".";
            vprint ModFinder, 5: "            System has "*Sprint(Ncols(lhs))*" equations involving "*Sprint(Nrows(lhs))*" variables.";

            eqnnum := Ncols(lhs);
            F`NumberOfEquations +:= eqnnum;
            if eqnnum gt maxeqnnum then
                maxeqnnum := eqnnum;
            end if;
            if auxvars gt maxauxvars then
                maxauxvars := auxvars;
            end if;
            Append(~eqnsizes, <Nrows(lhs), Ncols(lhs)>);

            vprint ModFinder, 5: "            System density: "*Sprint(ChangePrecision(Density(lhs),2))*".";
            vprint ModFinder, 5: "            Checking consistency of system.";

            //we don't want to solve the whole system and check consistency. it's enough to know what the solution looks like if there is one.
            tsys := Cputime();

            if solver eq "Mine" then

                E,T := EchelonForm(Transpose(lhs));
                u := rhs*Transpose(T);
                newdetvars := {@ @};
                for i:=1 to #nondetdepvars do
                    supp := Support(E[i]);
                    if supp subset {1..#nondetdepvars} and #supp eq 1 then
                        k := SetToSequence(supp)[1];
                        sol := u[k];
                        newvar := nondetdepvars[k];
                        newdetvars join:={@newvar@};
                        F`DeterminedVariables join:={@ newvar @};
                        F`DeterminedVariablesMap[newvar] := sol;
                    end if;
                end for;

            elif solver eq "Magma" then

                tsys := Cputime();
                t,sol,ker := IsConsistent(lhs,rhs);

                if t eq false then
                    print "System is inconsistent: There is no submodule with this abstract structure.";
                    return;
                end if;

                newdetvars := {@@};
                for i:=1 to #nondetdepvars do
                    if exists{ b : b in Basis(ker) | i in Support(b)} then
                        continue;
                    end if;
                    newdetvars join:={@ nondetdepvars[i] @};
                    F`DeterminedVariables join:={@ nondetdepvars[i] @};
                    F`DeterminedVariablesMap[nondetdepvars[i]] := sol[i];
                end for;

            end if;

            vprint ModFinder, 5:  "            Finished. Time: "*Sprint(Cputime(tsys))*".";

            vprint ModFinder, 5: "            System is consistent. Analyzing kernel.";

            if #newdetvars ne 0 then
                progress := true;
                update := true;

                SubAbsHelper_NonDetRows(~F);

                //remove all newly determined variables from datums
                while exists(x){x : x in nonredundantdatums | x[1] in newdetvars} do
                    nonredundantdatums diff:={x};
                end while;

                vprint ModFinder, 5: "            *** Could determine "*Sprint(#newdetvars)*" variables: "*Sprint(newdetvars)*" ***.";
                usefulsystemscounter +:= 1;

                vprint ModFinder, 5: "            Status: "*Sprint(ChangePrecision(100.0*#F`DeterminedVariables/#F`AbstractStructure`Variables,4))*"%.";
                PrintPercentage(#F`DeterminedVariables, #F`AbstractStructure`Variables);



                if updateinterval ne 0 and systemscounter mod updateinterval eq 0 then
                    break;
                end if;
            else

                if saveuseless then
                    uselesssystems join:={MD5(lhs)};
                end if;

                vprint ModFinder, 5: "            System is useless.";

                vprint ModFinder, 5: "            Status: "*Sprint(ChangePrecision(100.0*#F`DeterminedVariables/#F`AbstractStructure`Variables,4))*"%.";
                PrintPercentage(#F`DeterminedVariables, #F`AbstractStructure`Variables);
            end if;

        end while;

        if not progress and not InvolveAuxialiary then
            vprint ModFinder, 5: "Switching on auxiliary variables.";
            progress := true;
            InvolveAuxialiary := true;
            trieddatums:={};
        end if;

    end while;

    if #F`DeterminedVariables ne #F`AbstractStructure`Variables then
        print "CANNOT determine submodule.";
        return;
    end if;

	if GetVerbose("ModFinder") ge 3 then
		print "";
		print "Could determine all variables.";
		print "    Time: "*Sprint(Cputime(zeit));
        print "    Number of systems: "*Sprint(systemscounter);
        print "    Number of useful systems: "*Sprint(usefulsystemscounter);
		print "    System efficiency: "*Sprint( ChangePrecision(1.0*usefulsystemscounter/systemscounter, 4));
		maxsize := Maximum( [ x[1]*x[2] : x in eqnsizes] );
		t := exists(y){ x : x in eqnsizes | x[1]*x[2] eq maxsize};
		maxsize := y;
		print "    Maximal system size: "*Sprint(y);
		avgsize := Sqrt(ArraySum([ x[1]*x[2] : x in eqnsizes])/#eqnsizes);
		print "    Average system size: "*Sprint(ChangePrecision(avgsize,4));
        print "    Number of equations: "*Sprint(F`NumberOfEquations);
        print "    Equation efficiency: "*Sprint(ChangePrecision(1.0*#F`AbstractStructure`Variables/F`NumberOfEquations,4));
	end if;

end intrinsic;

//====================================================================
intrinsic Concretize(F::ModFinder) -> List
/*
    History:
        Saturday, September 21, 2013 21:18:52: Initial.
*/
{If all variables of F have been determined, instantiate F as a concrete subspace of its base module.}

	Ubasis := [**];
    K:=BaseRing(F`BaseModule);
	for r in [1..F`AbstractStructure`SubspaceDimension] do
		v := Zero(F`BaseModule`VectorSpace);
		for c in Support(F`AbstractStructure`KnownPart[r]) do
			v[c] := K!F`AbstractStructure`KnownPart[r][c];
		end for;
		for c in Support(F`AbstractStructure`UnknownPart[r]) do
			varnum := F`AbstractStructure`UnknownPart[r][c];
			v[c] := F`DeterminedVariablesMap[varnum];
		end for;
        z := New(ModGrOldElt);
        z`VectorSpaceElement := v;
        z`Parent := F`BaseModule;
		Append(~Ubasis, z);
	end for;

    return Ubasis;

end intrinsic;
