/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Intrinsics around reflections.
*/


//============================================================================
/*
    New attributes.
*/

declare attributes GrpMat:
    IsReflectionGroup,
    IsDiagonalizableReflectionGroup,
    Reflections,
    ReflectionLibrary,
    ReflectionLibraryFlat,
    ReflectionLibraryClasses,
    ReflectionClasses,
    NumberOfReflectionClasses,
    NumberOfReflections,
    RootOperators,
    CorootOperators;

//============================================================================
/*
    New structures
*/
    /*
        The reflection data record; the entries should be self-explanatory.
        We have one for matrices and one for matrix group elements. The latter contains a bit more data.
    */
    algmatrefldata := recformat<IsReflection:BoolElt, Element:AlgMatElt, Hyperplane:ModTupFld, Root:ModTupFldElt, Coroot:ModTupFldElt, Order:RngIntElt, IsDiagonalizable:BoolElt, Eigenvalue:FldElt, DualElement:AlgMatElt>;
    grpmatrefldata := recformat<IsReflection:BoolElt, Element:GrpElt, Word:SeqEnum, Class:RngIntElt, ReflectionClass:RngIntElt, Hyperplane:ModTupFld, Root:ModTupFldElt, Coroot:ModTupFldElt, Order:RngIntElt, IsDiagonalizable:BoolElt, Eigenvalue:FldElt, ID:Tup, ReflectionNumber:RngIntElt, DualElement:GrpMatElt>;

//============================================================================
intrinsic FixedSpace(M::AlgMatElt) -> ModTupFld
/*
    History:
        Friday, August 09, 2013 16:32:16: Initial.
*/
{The fixed space of M.}

    if not IsQuadratic(M) then
        error "Matrix has to be quadratic.";
    end if;
    K := BaseRing(M);
	n := Ncols(M);
    return Kernel(IdentityMatrix(K,n) - Matrix(M));

end intrinsic;

//============================================================================
intrinsic FixedSpace(M::GrpMatElt) -> ModTupFld
/*
    History:
        Friday, August 09, 2013 16:34:10: Initial.
*/
{The fixed space of M.}

    return FixedSpace(Matrix(M));

end intrinsic;

//============================================================================
intrinsic CoFixedSpace(M::AlgMatElt) -> ModTupFld
/*
    History:
        Friday, August 09, 2013 16:35:54: Initial.
*/
{The co-fixed space, Im(id-M).}

    if not IsQuadratic(M) then
        error "Matrix has to be quadratic.";
    end if;
    K := BaseRing(M);
	n := Ncols(M);
    return Image( IdentityMatrix(K,n) - Matrix(M));

end intrinsic;

//============================================================================
intrinsic CoFixedSpace(M::GrpMatElt) -> ModTupFld
/*
    History:
        Friday, August 09, 2013 16:36:50: Initial.
*/
{The co-fixed space, Im(id-M).}

    return CoFixedSpace(Matrix(M));

end intrinsic;

//============================================================================
intrinsic IsReflection(M::AlgMatElt : InvertibleCheck:=true) -> BoolElt
/*
    History:
        Friday, August 09, 2013 15:21:00: Initial.
*/
{True iff M is a reflection}

    if InvertibleCheck and not IsInvertible(M) then
        return false;
    end if;
    return Dimension(CoFixedSpace(M)) eq 1;

end intrinsic;



//============================================================================
intrinsic IsReflection(G::GrpMatElt) -> BoolElt
/*
    History:
        Friday, August 09, 2013 15:25:52: Initial.
*/
{True iff M is a reflection}

    return IsReflection(Matrix(G) : InvertibleCheck:=false);

end intrinsic;

//============================================================================
intrinsic ReflectionData(M::AlgMatElt) -> Rec
/*
    History:
        Friday, May 17, 2013 1:08:43 PM: Initial. Ported from old CHAMP.

        Friday, August 09, 2013 15:09:54: New record format.
*/
{If M is a reflection, then the hyperplane, a root, and admissible corresponding coroot (as in my thesis), the order, true/false depending on whether M is diagonalizable and in case it is diagonalizable, its non-trivial eigenvalue are returned. This works in any characteristic!}

	K := BaseRing(M);
	n := Ncols(M);
	result := rec<algmatrefldata|Element:=M>;

    result`IsReflection := IsReflection(M);
    if result`IsReflection eq false then
        return result;
    end if;

    result`Hyperplane := FixedSpace(M);
    result`Root := Vector(K,n,Eltseq(CoFixedSpace(M).1)); //to get the right parent
    result`Order := Order(M);

    evals := Eigenvalues(M);
	if IsDiagonalizable(M:Evals:=evals) eq false then
        result`IsDiagonalizable := false;
    else
		a := [ eig[1] : eig in evals | eig[1] ne 1 ][1];
		result`Eigenvalue := a;
        result`IsDiagonalizable := true;
	end if;

    //compute coroot
	crt := [];
	N := 1;
	f := IdentityMatrix(K,n) - Matrix(M);
	while result`Root[N] eq 0 do
		N +:= 1;
	end while;
	for i:=1 to n do
		row := f[i];
		scalar := row[N]/result`Root[N];
		Append(~crt, scalar);
	end for;
	result`Coroot := Vector(K,n,crt);

	return result;

end intrinsic;


//============================================================================
intrinsic ReflectionData(M::GrpMatElt) -> Rec
/*
    History:
        Friday, May 17, 2013 1:08:43 PM: Initial. Ported from old CHAMP.

        Friday, August 09, 2013 15:10:05: New record format.
*/
{Same as for AlgMatElt.}

	K := BaseRing(M);
	n := Ncols(M);
	result := rec<grpmatrefldata|Element:=M>;

    result`IsReflection := IsReflection(M);
    if result`IsReflection eq false then
        return result;
    end if;

    result`Hyperplane := FixedSpace(M);
    result`Root := Vector(K,n,Eltseq(CoFixedSpace(M).1)); //to get the right parent
    result`Order := Order(M);

    evals := Eigenvalues(M);
	if IsDiagonalizable(M:Evals:=evals) eq false then
        result`IsDiagonalizable := false;
    else
		a := [ eig[1] : eig in evals | eig[1] ne 1 ][1];
		result`Eigenvalue := a;
        result`IsDiagonalizable := true;
	end if;

    //compute coroot
	crt := [];
	N := 1;
	f := IdentityMatrix(K,n) - Matrix(M);
	while result`Root[N] eq 0 do
		N +:= 1;
	end while;
	for i:=1 to n do
		row := f[i];
		scalar := row[N]/result`Root[N];
		Append(~crt, scalar);
	end for;
	result`Coroot := Vector(K,n,crt);

    G := Parent(M);
    if assigned G`Classes then
    	if not assigned G`ClassMap then
    		G`ClassMap := ClassMap(G);
    	end if;
        result`Class := G`ClassMap(M);
    end if;

    result`DualElement := Transpose(M^-1);

	return result;

end intrinsic;

//============================================================================
intrinsic IsReflectionGroup(~G::GrpMat)
/*
    History:
        Thursday, August 08, 2013 11:46:30: Initial.

        Friday, August 09, 2013 16:46:00: Implemented.
*/
{True iff G is generated by reflections.}

    if assigned G`IsReflectionGroup then
        return;
    end if;

    ReflectionClasses(~G);
    refls := {@@};
    for i in G`ReflectionClasses do
        refls join:= Class(G, ClassRepresentative(G, i));
    end for;
    G`Reflections := SetToSequence(refls);
    G`IsReflectionGroup := sub<G|refls> eq G;

end intrinsic;

//============================================================================
intrinsic IsReflectionGroup(G::GrpMat) -> BoolElt
/*
    History:
        Friday, August 09, 2013 16:46:19: Initial.
*/
{True iff G is generated by reflections.}

    IsReflectionGroup(~G);
    return G`IsReflectionGroup;

end intrinsic;

//============================================================================
intrinsic IsDiagonalizable(M::AlgMatElt : Evals:=false) -> BoolElt
/*
    History:
        Friday, May 17, 2013 1:08:43 PM: Initial. Ported from old CHAMP.
*/
{Checks if a matrix M is diagonalizable. If the eigenvalues of M are already known, they can optionally be given as Evals.}

	if Type(Evals) eq BoolElt then
		Evals := Eigenvalues(M);
	end if;

	sum := 0;
	for eig in Evals do
		sum +:= Dimension(Eigenspace(M,eig[1]));
	end for;
	if sum eq Ncols(M) then
		return true;
	else
		return false;
	end if;

end intrinsic;

//============================================================================
intrinsic IsDiagonalizable(M::GrpMatElt : Evals:=false ) -> BoolElt
/*
    History:
        Friday, June 21, 2013 11:08:57: Initial.
*/
{Same as for AlgMatElt.}

    return IsDiagonalizable(Matrix(M) : Evals:=Evals );

end intrinsic;

//============================================================================
intrinsic IsDiagonalizableReflectionGroup(~G::GrpMat)
/*
    History:
        Wednesday, July 03, 2013 16:38:44: Initial.

        Friday, August 09, 2013 16:50:13: Updated.
*/
{True iff G is a reflection group (i.e., generated by reflections) and every reflection is diagonalizable.}

    ReflectionClasses(~G);
    if not IsReflectionGroup(G) then
        G`IsDiagonalizableReflectionGroup := false;
        return;
    end if;
    refls := {@@};
    for s in G`ReflectionClasses do
        if not IsDiagonalizable(s) then
            G`IsDiagonalizableReflectionGroup := false;
            return;
        end if;
    end for;

    G`IsDiagonalizableReflectionGroup := true;

end intrinsic;

//============================================================================
intrinsic IsDiagonalizableReflectionGroup(G::GrpMat) -> BoolElt
/*
    History:
        Friday, August 09, 2013 16:50:37: Initial.
*/
{True iff G is a reflection group (i.e., generated by reflections) and every reflection is diagonalizable.}

    IsDiagonalizableReflectionGroup(~G);
    return G`IsDiagonalizableReflectionGroup;

end intrinsic;

//============================================================================
intrinsic ReflectionClasses(~G::GrpMat)
/*
    History:
        Friday, August 09, 2013 15:13:18: Initial.
*/
{The classes of reflections in the same numbering as G`Classes.}

    if assigned G`ReflectionClasses then
        return;
    end if;

    ReflectionLibrary(~G);  //we call reflection library as this sorts the reflections

end intrinsic;

//============================================================================
intrinsic ReflectionClasses(G::GrpMat) -> SeqEnum
/*
    History:
        Friday, August 09, 2013 15:15:53: Initial.
*/
{}

    ReflectionClasses(~G);
    return G`ReflectionClasses;

end intrinsic;


//============================================================================
intrinsic Reflections(~G::GrpMat)
/*
    History:
        Friday, August 09, 2013 16:54:13: Initial.
*/
{The set of reflections of G.}

    if assigned G`Reflections then
        return;
    end if;

    ReflectionLibrary(~G);  //we call reflection library as this sorts the reflections

end intrinsic;

//============================================================================
intrinsic Reflections(G::GrpMat) -> SetEnum
/*
    History:
        Friday, August 09, 2013 17:04:20: Initial.
*/
{The set of reflections of G.}

    Reflections(~G);
    return G`Reflections;

end intrinsic;

//============================================================================
intrinsic ReflectionLibrary(~G::GrpMat : Sorting:="Quick", Verbose:=false)
{Reflections with all data sorted hierachically by hyperplane orbit, hyperplanes, reflections.}

    if assigned G`ReflectionLibrary then
        return;
    end if;

    if Verbose then
    	print "Building reflection library.";
    end if;

    //classes of reflections
    Classes(~G);
    G`ReflectionClasses := [ ClassRepresentative(G,i) : i in [1..#G`Classes] | IsReflection(ClassRepresentative(G,i)) ];
    G`NumberOfReflectionClasses := #G`ReflectionClasses;

    //reflections (this list will be reordered)
    G`Reflections := [];
    for i in G`ReflectionClasses do
        G`Reflections cat:=SetToSequence(IndexedSetToSet(Class(G, ClassRepresentative(G, i))));
    end for;

    Reflections(~G);

    reflswithdata := [* ReflectionData(g) : g in G`Reflections *];
    reflswithdatarefls := [ s`Element : s in reflswithdata ];

    if Sorting eq "RWSGroup" then
        RWSGroup(~G);
    elif Sorting eq "FPGroup" then
        FPGroup(~G);
    end if;

    //reflssorted is the a list of numbers corresponding to G`Reflections, resorted in some or the other way.
    if Sorting eq "None" then
        reflssorted := {@ i : i in [1..#G`Reflections] @};
    elif Sorting eq "Quick" then
        reflssorted := {@ @};
        for i:=1 to Ngens(G) do
            for j:=1 to Order(G.i) do
                k := Position(reflswithdatarefls, G.i^j);
                if k eq 0 then
                    continue;
                end if;
                reflssorted join:={@k@};
                reflswithdata[k]`Word := [ i : l in [1..j] ];
            end for;
        end for;
        Diff(~reflssorted, {@ 0 @});
        reflssorted join:={@ i : i in [1..#G`Reflections] @};
    else
        for i:=1 to #G`Reflections do
            reflswithdata[i]`Word := ElementToWord(reflswithdata[i]`Element : Method:=Sorting);
        end for;
        words := {@ reflswithdata[i]`Word : i in [1..#G`Reflections] @};
        //now we prepare the words for sorting; lexicographically by word length, sign list, absolute value list, and word itself.
        wordsprep := [ <#w, [ -Sign(w[i]) : i in [1..#w] ], [Abs(w[i]) : i in [1..#w]], w > : w in words ];
        wordsprepsorted := Sort(wordsprep);
        reflssorted := {@ Position(words, wordsprepsorted[i][4]) : i in [1..#words] @};
    end if;

    G`ReflectionLibrary := [ ];

    Omegacount := 0;
    while not IsEmpty(reflssorted) do
        //select Omega occording to sorting
        //Omega is hyperplane orbit of first reflection in new order
        Omega := Orbit(G, reflswithdata[reflssorted[1]]`Hyperplane);
        Omegacount +:= 1;
        //find all reflections whose hyperplanes lies in Omega
        Omegarefls := {@ i : i in reflssorted | reflswithdata[i]`Hyperplane in Omega @};
        Diff(~reflssorted, Omegarefls);

        //print "Omegarefls: "*Sprint(Omegarefls);

        OmegaHyperplanes := [ ]; //will be the hierarchy of hyperplanes in a fixed Omega

        Hcount := 0;
        while not IsEmpty(Omegarefls) do
            //select hyperplane in Omega according to sorting
            H := reflswithdata[Omegarefls[1]]`Hyperplane;
            Hcount +:= 1;
            Hrefls := {@ i : i in Omegarefls | reflswithdata[i]`Hyperplane eq H @};
            Diff(~Omegarefls, Hrefls);

            //print "Hrefls: "*Sprint(Hrefls);

            //now resort the reflections with hyperplane H
            Hreflssorted := [ ];
            orders := Sort({@ reflswithdata[i]`Order : i in Hrefls @});
            //print "Orders: "*Sprint(orders);
            reflcount := 0;
            while not IsEmpty(orders) do
                o := orders[1];
                Diff(~orders, {@o@});
                orefls := [* i : i in Hrefls | reflswithdata[i]`Order eq o *];
                //print "Orefls: "*Sprint(orefls);
                hasroot, ozeta := HasRootOfUnity(o, BaseRing(G));
                //reflections are diagonalizable iff ozeta in K
                if not hasroot then
                    for i in orefls do
                        reflcount +:=1;
                        reflswithdata[i]`ID := <Omegacount, Hcount, reflcount>;
                        Hreflssorted cat:=[ reflswithdata[i] ];
                    end for;
                else
                    for r:=1 to o-1 do
                        for j in [ i : i in orefls | reflswithdata[i]`Eigenvalue eq ozeta^r ] do
                            reflcount +:=1;
                            reflswithdata[j]`ID := <Omegacount, Hcount, reflcount>;
                            Hreflssorted cat:= [ reflswithdata[j] ]; //all reflections with eigenvalue ozeta^r
                        end for;
                    end for;
                    //print Hreflssorted;
                end if;
            end while;
            OmegaHyperplanes cat:=[ Hreflssorted ];
        end while;
        G`ReflectionLibrary cat:=[ OmegaHyperplanes ];

    end while;

    //compatibly set
    G`ReflectionClasses := [];
    G`ReflectionLibraryClasses := [];
    cl := {};
    reflclnumber := [];
    for i:=1 to #G`ReflectionLibrary do
        for j:=1 to #G`ReflectionLibrary[i] do
            for k:=1 to #G`ReflectionLibrary[i][j] do
                s := G`ReflectionLibrary[i][j][k];
                if s`Class notin cl then
                    Append(~G`ReflectionClasses, s`Element);
                    Append(~G`ReflectionLibraryClasses, s);
                    Append(~reflclnumber, s`Class);
                    cl join:={s`Class};
                end if;
                G`ReflectionLibrary[i][j][k]`ReflectionClass := Position(reflclnumber, s`Class);
            end for;
        end for;
    end for;

    G`ReflectionLibraryFlat := FlatFixed(G`ReflectionLibrary);
    G`NumberOfReflections := #G`ReflectionLibraryFlat;
    G`Reflections := [ G`ReflectionLibraryFlat[i]`Element : i in [1..#G`ReflectionLibraryFlat] ];
    for i:=1 to #G`ReflectionLibraryFlat do
    	G`ReflectionLibraryFlat[i]`ReflectionNumber := i;
    end for;


end intrinsic;

//============================================================================
intrinsic ReflectionLibrary(G::GrpMat) -> SeqEnum
/*
    History:
        Monday, August 12, 2013 10:43:55: Initial.
*/
{}

    ReflectionLibrary(~G);
    return G`ReflectionLibrary;

end intrinsic;

//============================================================================
intrinsic ReflectionID(g::GrpMatElt) -> Tup
/*
    History:
        Saturday, August 10, 2013 16:12:12: Initial
*/
{}
    G := Parent(g);

    res := exists(pos){i : i in [1..#G`Reflections] | G`ReflectionLibraryFlat[i]`Element eq g};
    if res ne false then
        return G`ReflectionLibraryFlat[pos]`ID;
    else
        error "Element not in reflections library.";
    end if;

end intrinsic;

//============================================================================
intrinsic MinimalGeneratingSystemsOfReflections(G::GrpMat) -> SetEnum
/*
    History:
        Monday, September 23, 2013 19:39:29: Initial.
*/
{}

    Reflections(~G);
    gens := {};
    d := 0;
    for i:=1 to #G do
        for S in CartesianProduct([ G`Reflections : j in [1..i] ]) do
            Sset := {S[j] : j in [1..i]};
            if sub<G|Sset> eq G then
                d := i;
                break i;
            end if;
        end for;
    end for;


    for S in CartesianProduct([ G`Reflections : j in [1..d] ]) do
        Sset := {S[j] : j in [1..d]};
        if sub<G|Sset> eq G then
            gens join:={Sset};
        end if;
    end for;

    return gens;

end intrinsic;

//==============================================================================
intrinsic IsRegular(G::GrpMat, v::ModTupFldElt) -> BoolElt
/*
    Intrinsic: IsRegular

    Description:
        Here, +v+ is an element of the vector space +V+ on which the matrix group +G+ acts. This intrinsic returns true if and only if +v+ is a regular vector with respect to this action, i.e., the stabilizer of +v+ is trivial.

    History:
        * Monday, April 14, 2014 at 12:27:10: Initial
*/
{}

    return #Stabilizer(G, v) eq 1;

end intrinsic;



//==============================================================================
intrinsic CorootOperator(G::GrpMat, i::RngIntElt, p::RngMPolElt : SaveToTable:=true) -> RngMPolElt
/*
    Description:
        The coroot operator as defined in Gordon, Baby Verma Modules, 3.5, for the reflection of +G+ with index +i+ in ReflectionLibraryFlat.
*/
{}

    P := Parent(p);
    result := Zero(P);

    for t in Terms(p) do

        tmon := Monomials(t)[1];
        tcoeff := Coefficients(t)[1];
        if tmon eq 1 then
            continue;   //operator is zero on scalars
        end if;

        //find in table
        /*find, op := IsDefined(G`CorootOperators[i], tmon);
        if find then
            result := result + tcoeff*op;
            continue;
        end if;*/

        //compute
        exp := Exponents(tmon);

        //if element of V
        if Degree(tmon) eq 1 then
            j := GetVariableNumber(tmon);
            op := G`ReflectionLibraryFlat[i]`Coroot[j];
            /*if SaveToTable then
                Set(~G`CorootOperators[i], tmon, op);
            end if;*/
            result := result + tcoeff*op;
            continue;
        end if;

        //use formula
        //split tmon into leftpart*rightpart with rightpart of degree 1
        //find right most non-zero position in exponent sequence
        n:=Rank(P);
        while exp[n] eq 0 do
            n := n - 1;
        end while;

        leftpartexp := exp;
        leftpartexp[n] := leftpartexp[n] - 1;
        leftpart := Monomial(P, leftpartexp);
        rightpart := P.n;

        op := CorootOperator(G,i,leftpart)*rightpart + leftpart*CorootOperator(G,i,rightpart) - CorootOperator(G,i,leftpart)*CorootOperator(G,i,rightpart);

        /*if SaveToTable then
            Set(~G`CorootOperators[i], tmon, op);
        end if;*/

        result := result + tcoeff*op;

    end for;

    return result;

end intrinsic;

//==============================================================================
intrinsic RootOperator(G::GrpMat, i::RngIntElt, p::RngMPolElt : SaveToTable:=true) -> RngMPolElt
/*
    Description:
        The root operator as defined in Gordon, Baby Verma Modules, 3.5, for the reflection of +G+ with index +i+ in ReflectionLibraryFlat.
*/
{}

    P := Parent(p);
    result := Zero(P);

    for t in Terms(p) do

        tmon := Monomials(t)[1];
        tcoeff := Coefficients(t)[1];
        if tmon eq 1 then
            continue;   //operator is zero on scalars
        end if;

        //find in table
        /*find, op := IsDefined(G`RootOperators[i], tmon);
        if find then
            result := result + tcoeff*op;
            continue;
        end if;*/

        //compute
        exp := Exponents(tmon);

        //if element of V
        if Degree(tmon) eq 1 then
            j := GetVariableNumber(tmon);
            op := G`ReflectionLibraryFlat[i]`Root[j];
            /*if SaveToTable then
                Set(~G`RootOperators[i], tmon, op);
            end if;*/
            result := result + tcoeff*op;
            continue;
        end if;

        //use formula
        //split tmon into leftpart*rightpart with rightpart of degree 1
        //find right most non-zero position in exponent sequence
        n:=Rank(P);
        while exp[n] eq 0 do
            n := n - 1;
        end while;

        leftpartexp := exp;
        leftpartexp[n] := leftpartexp[n] - 1;
        leftpart := Monomial(P, leftpartexp);
        rightpart := P.n;

        op := RootOperator(G,i,leftpart)*rightpart + leftpart*RootOperator(G,i,rightpart) - RootOperator(G,i,leftpart)*RootOperator(G,i,rightpart);

        /*if SaveToTable then
            Set(~G`RootOperators[i], tmon, op);
        end if;*/

        result := result + tcoeff*op;

    end for;

    return result;

end intrinsic;
