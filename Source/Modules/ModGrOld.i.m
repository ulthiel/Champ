/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
    An old graded module structure (used by ModFinder)
*/


//==============================================================================
intrinsic GradedModuleOld(Q::SeqEnum, Qdegs::SeqEnum, Algdegs::SeqEnum) -> ModGrOld
/*
    Intrinsic: GradedModule

    Creates a graded module

    Declaration:
        :intrinsic GradedModule(Q::SeqEnum, Qdegs::SeqEnum, Algdegs::SeqEnum) -> ModGrOld

    Parameters:
       Q - a list of matrices defining the action of the generators in +Algdegs+
       Qdegs - the row degrees of the matrices in +Q+
       Algdegs - degrees of the generators

    Description:
        Create a graded module defined by the matrices in +Q+ with row degrees as defined in +Qdegs+ and with algebra generator degrees as defined in +Algdegs+.
*/
{Create a graded module defined by the matrices in Q with row degrees as defined in Qdegs and with algebra generator degrees as defined in Algdegs.}

    M := New(ModGrOld);
    M`Matrices := Q;
    M`RowDegrees := Qdegs;
    M`MatrixDegrees := Algdegs;
    M`Dimension := Ncols(Q[1]);
    M`BaseRing := BaseRing(Q[1]);
    M`VectorSpace := KSpace(M`BaseRing, M`Dimension);
    if Type(Q[1]) eq MtrxSprs then
        M`Rep := "Sparse";
    else
        M`Rep := "Dense";
    end if;
    return M;

end intrinsic;

//==============================================================================
intrinsic Print(M::ModGrOld)
/*
    Intrinsic: Print

    Prints information about a graded module

    Declaration:
        :intrinsic Print(M::ModGrOld)

    Parameters:
       M - a graded module
*/
{}

    printf "Graded module of dimension "*Sprint(M`Dimension)*" over an algebra with generator degrees "*Sprint(M`MatrixDegrees)*" over "*Sprint(M`BaseRing, "Minimal")*".";

end intrinsic;

//==============================================================================
intrinsic Ngens(M::ModGrOld) -> RngIntElt
/*
    Intrinsic: Ngens

    Number of generators of a graded module

    Declaration:
        :intrinsic Ngens(M::ModGrOld) -> RngIntElt

    Parameters:
       M - a graded module
*/
{}

    return #M`Matrices;

end intrinsic;

//==============================================================================
intrinsic BaseRing(M::ModGrOld) -> Rng
/*
    Intrinsic: BaseRing

    Base ring of a graded module

    Declaration:
        :intrinsic BaseRing(M::ModGrOld) -> Rng

    Parameters:
       M - a graded module
*/
{}

    return M`BaseRing;

end intrinsic;

//==============================================================================
intrinsic Dimension(M::ModGrOld) -> Rng
/*
    Intrinsic: Dimension

    Dimension of a graded module

    Declaration:
        :intrinsic Dimension(M::ModGrOld) -> Rng

    Parameters:
       M - a graded module
*/
{}

    return M`Dimension;

end intrinsic;

//==============================================================================
intrinsic VectorSpace(M::ModGrOld) -> ModTupFld
/*
    Intrinsic: VectorSpace

    The underlying vector space of a graded module

    Declaration:
        :intrinsic VectorSpace(M::ModGrOld) -> ModTupFld

    Parameters:
       M - a graded module
*/
{}

    return M`VectorSpace;

end intrinsic;

//==============================================================================
/*
    Intrinsic: RModule

    The underlying (ungraded) RModule of a graded module

    Declaration:
        :intrinsic RModule(~M::ModGrOld)
        :intrinsic RModule(M::ModGrOld) -> ModRng

    Parameters:
       M - a graded module
*/
intrinsic RModule(~M::ModGrOld)
{}

    if M`Rep eq "Dense" then
        M`RModule := RModule(M`Matrices);
    else
        M`RModule := RModule([ Matrix(M`Matrices[i]) : i in [1..#M`Matrices] ]);
    end if;

end intrinsic;

//==============================================================================
intrinsic RModule(M::ModGrOld) -> ModRng
{}

    RModule(~M);
    return M`RModule;

end intrinsic;


//==============================================================================
intrinsic PoincareSeries(M::ModGrOld) -> FldFunRatUElt
/*
    Intrinsic: PoincareSeries

    The Poincare series of a graded module

    Declaration:
        :intrinsic PoincareSeries(M::ModGrOld) -> RngSerLaur

    Parameters:
       M - a graded module
*/
{}

    R<q> := RationalFunctionField(Rationals());
    P := Zero(R);
    for i:=1 to #M`RowDegrees do
        P +:= q^M`RowDegrees[i];
    end for;

    return P;

end intrinsic;

//==============================================================================
intrinsic Density(M::ModGrOld) -> FldReElt
/*
    Intrinsic: Density

    The total density of the action matrices of a graded module

    Declaration:
        :intrinsic Density(M::ModGrOld) -> FldReElt

    Parameters:
       M - a graded module
*/
{}

    nonzeros := &+[ NumberOfNonZeroEntries(M`Matrices[i]) : i in [1..Ngens(M)] ];
    total := Ngens(M)*Dimension(M)^2;
    return RealField()!(nonzeros/total);

end intrinsic;

//==============================================================================
intrinsic NumberOfNonZeroEntries(M::ModGrOld) -> RngIntElt
/*
    Intrinsic: NumberOfNonZeroEntries

    The total number of non-zero entries of the action matrices of a graded module

    Declaration:
        :intrinsic NumberOfNonZeroEntries(M::ModGrOld) -> RngIntElt

    Parameters:
       M - a graded module
*/
{}

    return &+[ NumberOfNonZeroEntries(M`Matrices[i]) : i in [1..Ngens(M)] ];

end intrinsic;


//==============================================================================
intrinsic DimensionOfHomogeneousComponent(M::ModGrOld, d::RngIntElt) -> RngIntElt
/*
    Intrinsic: DimensionOfHomogeneousComponent

    The dimension of a homogeneous component of a graded module

    Declaration:
        :intrinsic DimensionOfHomogeneousComponent(M::ModGrOld, d::RngIntElt) -> RngIntElt

    Parameters:
       M - a graded module
       d - an integer

    Description:
    	The dimension of the homogeneous component of +M+ of degree +d+.
*/
{The dimension of the homogeneous component of M of degree d.}

    return #{ i : i in [1..Dimension(M)] | M`RowDegrees[i] eq d };

end intrinsic;

//==============================================================================
intrinsic HomogeneousData(~M::ModGrOld)
/*
    Intrinsic: HomogeneousData

    Attaches homogeneous data for a graded module

    Declaration:
        :intrinsic HomogeneousData(~M::ModGrOld)

    Parameters:
       M - a graded module

    Description:
    	Attach projections, embeddings and maps of and between homogeneous components of +M+.
*/
{Attach projections, embeddings and maps of and between homogeneous components of M.}

    if assigned M`HomogeneousComponents then
        return;
    end if;

    M`HomogeneousComponents := AssociativeArray(Universe({1}));
    M`HomogeneousComponentProjections := AssociativeArray(Universe({1}));
    M`HomogeneousComponentEmbeddings := AssociativeArray(Universe({1}));

    degrees := SequenceToSet(M`RowDegrees);
    for d in degrees do
        rowsofdegreed := [ i : i in [1..Dimension(M)] | M`RowDegrees[i] eq d ];
        M`HomogeneousComponents[d] := KSpace(BaseRing(M), #rowsofdegreed);

        //compute projection
        A := ZeroMatrix(BaseRing(M), Dimension(M), #rowsofdegreed);
        for i:=1 to Dimension(M) do
            p := Position(rowsofdegreed, i);
            if p eq 0 then
                continue;
            else
                A[i][p] := 1;
            end if;
        end for;
        M`HomogeneousComponentProjections[d] := hom<VectorSpace(M)->M`HomogeneousComponents[d]| A >;

        //compute embedding
        A := ZeroMatrix(BaseRing(M), #rowsofdegreed, Dimension(M));
        for i:=1 to #rowsofdegreed do
            A[i][rowsofdegreed[i]] := 1;
        end for;
        M`HomogeneousComponentEmbeddings[d] := hom<M`HomogeneousComponents[d] -> VectorSpace(M) | A>;
    end for;

    M`HomogeneousComponentMatrices := AssociativeArray(Universe({<1,1>}));
    degrees := SequenceToSet(M`RowDegrees);
    for i:=1 to Ngens(M) do
        for d in degrees do
            targetdegree := d + M`MatrixDegrees[i];
            rowsofdegreed := [ j : j in [1..Dimension(M)] | M`RowDegrees[j] eq d];
            rowsoftargetdegree :=  [ j : j in [1..Dimension(M)] | M`RowDegrees[j] eq targetdegree];
            degreeddimension := DimensionOfHomogeneousComponent(M,d);
            targetdegreedimension := DimensionOfHomogeneousComponent(M,targetdegree);
            if M`Rep eq "Sparse" then
                A := SparseMatrix(BaseRing(M), degreeddimension, targetdegreedimension);
            else
                A := ZeroMatrix(BaseRing(M), degreeddimension, targetdegreedimension);
            end if;
            for j in rowsofdegreed do
                for k in Support(M`Matrices[i][j]) do
                    A[Position(rowsofdegreed,j)][Position(rowsoftargetdegree, k)] := M`Matrices[i][j][k];
                end for;
            end for;
            M`HomogeneousComponentMatrices[<i,d>] := A;
        end for;
    end for;

end intrinsic;

//==============================================================================
intrinsic Print(x::ModGrOldElt)
/*
    Intrinsic: Print

    Prints a vector of a graded module

    Declaration:
        :intrinsic Print(x::ModGrOldElt)

    Parameters:
       x - a vector of a graded module
*/
{}

    printf "%o", x`VectorSpaceElement;

end intrinsic;

//==============================================================================
intrinsic '.'(M::ModGrOld, i::RngIntElt) -> ModGrOldElt
/*
    Intrinsic: '.'

    Returns a basis vector of a graded module

    Declaration:
        :intrinsic '.'(M::ModGrOld, i::RngIntElt) -> ModGrOldElt

    Parameters:
       i - an integer

    Description:
    	Returns the +i+-th basis vector of the underlying vector space of +M+
*/
{}
    x := New(ModGrOldElt);
    x`Parent := M;
    x`Degree := M`RowDegrees[i];
    x`VectorSpaceElement := M`VectorSpace.i;
    return x;

end intrinsic;

//==============================================================================
intrinsic Basis(M::ModGrOld) -> List
{}

    return [M.i : i in [1..Dimension(M)]];

end intrinsic;

//==============================================================================
intrinsic '+'(x::ModGrOldElt, y::ModGrOldElt) -> ModGrOldElt
/*
    Intrinsic: '+'

    Addition of vectors of a graded module

    Declaration:
        :intrinsic '+'(x::ModGrOldElt, y::ModGrOldElt) -> ModGrOldElt

    Parameters:
       x - a vector of a graded module
       y - a vector of a graded module
*/
{}

    z := New(ModGrOldElt);
    z`VectorSpaceElement := x`VectorSpaceElement + y`VectorSpaceElement;
    z`Parent := x`Parent;
    if assigned x`Degree and assigned y`Degree and x`Degree eq y`Degree then
        z`Degree := x`Degree;
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic '*'(a::RngElt, x::ModGrOldElt) -> ModGrOldElt
/*
    Intrinsic: '*'

    Scalar multiplication of a vector of a graded module

    Declaration:
        :intrinsic '*'(a::RngElt, x::ModGrOldElt) -> ModGrOldElt

    Parameters:
    	a - an element of the base ring of the graded module of +x+
    	x - a vector of a graded module
*/
{}

    z := New(ModGrOldElt);
    z`VectorSpaceElement := a*x`VectorSpaceElement;
    z`Parent := x`Parent;
    if assigned x`Degree then
        z`Degree := x`Degree;
    end if;

    return z;

end intrinsic;

//==============================================================================
intrinsic 'eq'(x::ModGrOldElt, y::ModGrOldElt) -> BoolElt
/*
    Intrinsic: 'eq'

    Equality of vectors of a graded module

    Declaration:
        :intrinsic 'eq'(x::ModGrOldElt, y::ModGrOldElt) -> BoolElt

    Parameters:
       x - a vector of a graded module
       y - a vector of a graded module
*/
{}

    return x`VectorSpaceElement eq y`VectorSpaceElement;

end intrinsic;

//==============================================================================
intrinsic Parent(x::ModGrOldElt) -> ModGrOld
/*
    Intrinsic: Parent

    The graded module of a vector of a graded module

    Declaration:
        :intrinsic Parent(x::ModGrOldElt) -> ModGrOld

    Parameters:
       x - a vector of a graded module
*/
{}

    return x`Parent;

end intrinsic;


//==============================================================================
intrinsic Support(x::ModGrOldElt) -> SetEnum
/*
    Intrinsic: Support

    The non-zero indices in the basis representation of a vector

    Declaration:
        :intrinsic Support(x::ModGrOldElt) -> SetEnum

    Parameters:
       x - a vector of a graded module
*/
{}

    return Support(x`VectorSpaceElement);

end intrinsic;

//==============================================================================
intrinsic IsHomogeneous(x::ModGrOldElt) -> BoolElt, RngIntElt
/*
    Intrinsic: IsHomogeneous

    Checks if a vector of a graded module is homogeneous

    Declaration:
        :intrinsic IsHomogeneous(x::ModGrOldElt) -> BoolElt, RngIntElt

    Parameters:
       x - a vector of a graded module
*/
{}
    if assigned x`Degree then
        return true, x`Degree;
    end if;

    if IsZero(x`VectorSpaceElement) then
        return false, _;
    end if;

    M := Parent(x);
    degrees := { };
    for i in Support(x`VectorSpaceElement) do
        d := M`RowDegrees[i];
        if not IsEmpty(degrees) and d notin degrees then
            return false,_;
        else
            degrees join:={d};
        end if;
    end for;

    x`Degree := d;
    return true, d;

end intrinsic;

//==============================================================================
intrinsic Degree(x::ModGrOldElt) -> RngIntElt
/*
    Intrinsic: Degree

    The degree of a homogeneous vector of a graded module

    Declaration:
        :intrinsic Degree(x::ModGrOldElt) -> RngIntElt

    Parameters:
       x - a vector of a graded module
*/
{}

    if assigned x`Degree then
        return x`Degree;
    else
        res, d := IsHomogeneous(x);
        if res then
            return d;
        else
            error "Element is not homogeneous.";
        end if;
    end if;

end intrinsic;

//==============================================================================
intrinsic ActionGenerator(M::ModGrOld, i::RngIntElt) -> Mtrx
/*
    Intrinsic: ActionGenerator

    Returns an action generator of a graded module

    Declaration:
        :intrinsic ActionGenerator(M::ModGrOld, i::RngIntElt) -> Mtrx

    Parameters:
       M - a graded module
       i - an integer
*/
{}

    return M`Matrices[i];

end intrinsic;


//==============================================================================
intrinsic '^'(x::ModGrOldElt, i::RngIntElt) -> ModGrOldElt
/*
    Intrinsic: '^'

    Action of an action generator on a vector of a graded module

    Declaration:
        :intrinsic '^'(x::ModGrOldElt, i::RngIntElt) -> ModGrOldElt

    Parameters:
       M - a graded module
       i - an integer

    Description:
    	The action of the +i+-th generator of the parent of +x+ on +x+.
*/
{The action of the i-th generator of the parent of x on x.}

    M := Parent(x);
    z := New(ModGrOldElt);
    z`VectorSpaceElement := x`VectorSpaceElement*M`Matrices[i];
    z`Parent := M;
    if assigned x`Degree and not IsZero(z`VectorSpaceElement) then
        z`Degree := x`Degree + M`MatrixDegrees[i];
    end if;
    return z;

end intrinsic;

//==============================================================================
intrinsic Spin(M::ModGrOld, U::List) -> List, Assoc
/*
    Intrinsic: Spin

    Spins a list homogenous vectors of a graded module

    Declaration:
        :intrinsic Spin(M::ModGrOld, U::List) -> List, Assoc
        :intrinsic Spin(V::ModGrOld, x::ModGrOldElt) -> ModGrOld

    Parameters:
       M - a graded module
       U - a list

    Description:
    	Spin the list +U+ of homogenous vectors of +M+ to the submodule of +M+ they generate and return the basis and the homogeneous components of this submodule.
*/
{Spin the subset U of M consisting of homogeneous vectors to a submodule of M and return the basis and the homogeneous components of this submodule.}

    for x in U do
        if not IsHomogeneous(x) then
            error "Only works for homogeneous submodules.";
        end if;
    end for;

    HomogeneousData(~M);
    auxspaces := AssociativeArray(Integers()); //this will be the homogeneous components of the submodule generated by Q
    dimcounter := 0;
    for d in SequenceToSet(M`RowDegrees) do
        Ud := [ M`HomogeneousComponentProjections[d](x`VectorSpaceElement) : x in U | Degree(x) eq d ];
        auxspaces[d] := sub<M`HomogeneousComponents[d]|Ud>;
        dimcounter +:= Dimension(auxspaces[d]);
    end for;

    tospin := [* *];
    for x in U do
        d := x`Degree;
        Append(~tospin, <d, M`HomogeneousComponentProjections[d](x`VectorSpaceElement)>);
    end for;

    while not IsEmpty(tospin) do
        //i := Random([1..#tospin]);

        //take the "shortest" to spin instead of random one
        lengths := [ #Sprint(tospin[i][2]) : i in [1..#tospin] ];
        minlength := Minimum(lengths);
        minpos := Position(lengths, minlength);
        i := minpos;

        y := tospin[i]; //y[1] is the degree of y[2]
        Remove(~tospin, i);
        for i:=1 to Ngens(M) do
            if y[1] + M`MatrixDegrees[i] notin M`RowDegrees then //it's zero. we don't have to compute this!
                continue;
            end if;
            z := <y[1]+M`MatrixDegrees[i], y[2]*M`HomogeneousComponentMatrices[<i,y[1]>]>; //image of y[1] under generator i
            if z[2] notin auxspaces[z[1]] then
                auxspaces[z[1]] +:= sub<M`HomogeneousComponents[z[1]]|z[2]>;
                Append(~tospin, z);
                dimcounter +:=1;
            end if;

            if GetVerbose("ModGrOld") ge 5 then
                PrintAndDelete("Current dimension: "*Sprint(dimcounter)*" ("*Sprint(#tospin)*" to spin)                  ");
            end if;
        end for;
    end while;

    //construct basis
    basis := [* *];
    for d in SequenceToSet(M`RowDegrees) do
        for b in Basis(auxspaces[d]) do
            z := New(ModGrOldElt);
            z`Parent := M;
            z`VectorSpaceElement := M`HomogeneousComponentEmbeddings[d](b);
            z`Degree := d;
            Append(~basis, z);
        end for;
    end for;

    if GetVerbose("ModGrOld") ge 5 then
        printf "\n";
    end if;

    return basis, auxspaces;

end intrinsic;

//==============================================================================
intrinsic Spin(V::ModGrOld, x::ModGrOldElt) -> ModGrOld
{Spin the homogeneous vector x.}

    return Spin(V, [* x *]);

end intrinsic;


//==============================================================================
intrinsic Quotient(M::ModGrOld, U::List : Rep:="Dense", PerformSpinning:=true, Alg:="Own") -> ModGrOld, List
/*
    Intrinsic: Quotient

    The quotient of a graded module by a submodule

    Declaration:
        :intrinsic Quotient(M::ModGrOld, U::List : Rep:="Dense", PerformSpinning:=true, Alg:="Own") -> ModGrOld, List

    Parameters:
       M - a graded module
       U - a list

    Options:
    	Rep - "Dense" or "Sparse"
    	PerformSpinning - true or false
    	Alg - "Own" or "Magma"

    Description:
    	The quotient of +M+ by the submodule generated by +U+.
*/
{The quotient of M by the submodule generated by U.}

    if IsEmpty(U) then
        return M,U;
    end if;

    HomogeneousData(~M);

    vprint ModGrOld, 5: "Computing quotient.";

    if GetVerbose("ModGrOld") ge 5 then
        IndentPush();
    end if;

    if Alg eq "Own" then
        if PerformSpinning then
            vprint ModGrOld, 5: "Performing spinning.";
            U, Ucpt := Spin(M,U);
            Udim := #U;
        else
            //In this case U is assumed to be a basis for the submodule generated by U.
            //The user should know that this is true.
            Ucpt := AssociativeArray(Integers());
            Udim := 0;
            for d in SequenceToSet(M`RowDegrees) do
                Ud := [ M`HomogeneousComponentProjections[d](x`VectorSpaceElement) : x in U | Degree(x) eq d ];
                Ucpt[d] := sub<M`HomogeneousComponents[d]|Ud>;
                Udim +:= Dimension(Ucpt[d]);
            end for;
        end if;

        quotientspaces := AssociativeArray(Integers());
        quotientmaps := AssociativeArray(Integers());
        sections := AssociativeArray(Integers());

        nops := #SequenceToSet(M`RowDegrees);
        nopcounter := 0;
        for d in SequenceToSet(M`RowDegrees) do
            nopcounter +:= 1;
            if GetVerbose("ModGrOld") ge 5 then
                PrintAndDelete("Computing vector space quotients: "*Sprint( ChangePrecision( 100.0*nopcounter/nops, 4))*"      ");
            end if;
            quotientspaces[d], quotientmaps[d], sections[d] := Quotient(M`HomogeneousComponents[d], Ucpt[d]);
        end for;

        quotientdegrees := []; //will be the row degrees of the quotient
        Qdim:=0; //will be the dimension of the quotient
        for d in SequenceToSet(M`RowDegrees) do
            for i:=1 to Dimension(quotientspaces[d]) do
                Append(~quotientdegrees, d);
                Qdim +:= 1;
            end for;
        end for;

        if GetVerbose("ModGrOld") ge 5 then
            printf "\n";
        end if;

        quotientactions := AssociativeArray(Universe({<1,1>}));
        nops := Ngens(M)*#SequenceToSet(quotientdegrees);
        nopcounter := 0;
        for i:=1 to Ngens(M) do
            for d in SequenceToSet(quotientdegrees) do
                nopcounter +:= 1;
                if GetVerbose("ModGrOld") ge 5 then
                    PrintAndDelete("Computing quotient actions: "*Sprint( ChangePrecision( 100.0*nopcounter/nops, 4))*"      ");
                end if;
                if not IsDefined(quotientspaces, d+M`MatrixDegrees[i]) then //otherwise action of x_i on Q_d is zero!
                    continue;
                end if;
                A := ZeroMatrix(BaseRing(M), Dimension(quotientspaces[d]), Dimension(quotientspaces[d+M`MatrixDegrees[i]]));
                for j:=1 to Dimension(quotientspaces[d]) do
                    v := sections[d](quotientspaces[d].j);
                    w := v*M`HomogeneousComponentMatrices[<i,d>];
                    if not IsZero(w) then
                        A[j] := quotientmaps[d+M`MatrixDegrees[i]](w);
                    end if; //otherwise A[j] remains zero
                end for;
                quotientactions[<i,d>] := A;
            end for;
        end for;

        if GetVerbose("ModGrOld") ge 5 then
            printf "\n";
        end if;

        //now construct the module from this homogeneous data
        Qmatrices := [ ZeroMatrix(BaseRing(M), Qdim, Qdim) : i in [1..Ngens(M)] ];
        nops := Ngens(M)*#quotientdegrees;
        nopcounter := 0;
        for i:=1 to Ngens(M) do
            for d in quotientdegrees do
                nopcounter +:= 1;
                if GetVerbose("ModGrOld") ge 5 then
                    PrintAndDelete("Combining homogeneous data: "*Sprint( ChangePrecision( 100.0*nopcounter/nops, 4))*"      ");
                end if;
                drows := [ j : j in [1..#quotientdegrees] | quotientdegrees[j] eq d ];
                targetdegree := d + M`MatrixDegrees[i];
                targetdegreerows := [ j : j in [1..#quotientdegrees] | quotientdegrees[j] eq targetdegree ];
                //print d,drows,targetdegree,targetdegreerows;
                if IsEmpty(targetdegreerows) then //action is zero in this case
                    continue;
                end if;

                //by construction above we can work with blocks
                dstart := Minimum(drows);
                dend := Maximum(drows);
                targetdegreestart := Minimum(targetdegreerows);
                targetdegreeend := Maximum(targetdegreerows);

                InsertBlock(~Qmatrices[i], quotientactions[<i,d>], dstart, targetdegreestart);

            end for;
        end for;

        if GetVerbose("ModGrOld") ge 5 then
            printf "\n";
        end if;

        if GetVerbose("ModGrOld") ge 5 then
            IndentPop();
        end if;

        return GradedModuleOld(Qmatrices, quotientdegrees, M`MatrixDegrees), U;

    elif Alg eq "Magma" then

        RModule(~M);
        U,iota := sub<M`RModule| [ x`VectorSpaceElement : x in U ]>;
        Q,q := quo<M`RModule|U>;
        Usubspace := sub<VectorSpace(M`RModule) | [ VectorSpace(M`RModule)!iota(VectorSpace(U).i) : i in [1..Dimension(U)]] >;
        Qspace, Qspaceproj, Qspacesection := Quotient(VectorSpace(M`RModule), Usubspace);
        Qdegrees := [];
        for i:=1 to Dimension(Qspace) do
            v := Qspacesection(Qspace.i);
            x := New(ModGrOldElt);
            x`Parent := M;
            x`VectorSpaceElement := v;
            Append(~Qdegrees, Degree(x));
        end for;

        Ubasis := [* *];
        for i:=1 to Dimension(Usubspace) do
            x := New(ModGrOldElt);
            x`Parent := M;
            x`VectorSpaceElement := Usubspace.i;
            Append(~Ubasis, x);
        end for;

        return GradedModuleOld(ActionGenerators(Q), Qdegrees, M`MatrixDegrees), Ubasis;

    end if;

end intrinsic;

//==============================================================================
intrinsic ChangeRing(M::ModGrOld, S::Rng) -> ModGrOld
/*
    Intrinsic: ChangeRing

    Change the base ring of a graded module

    Declaration:
        :intrinsic ChangeRing(M::ModGrOld, S::Rng) -> ModGrOld

    Parameters:
       M - a graded module
       S - a ring
*/
{}

    mats := [ ChangeRing(M`Matrices[i], S) : i in [1..#M`Matrices] ];
    return GradedModuleOld(mats, M`RowDegrees, M`MatrixDegrees);

end intrinsic;



//==============================================================================
intrinsic DecompositionInGradedGrothendieckGroup(M::ModGrOld, G::Grp, groupgens::SeqEnum : UseCharacters:=true) -> ModTupRngElt
/*
	Intrinsic: DecompositionInGradedGrothendieckGroup

    Declaration:
        :DecompositionInGradedGrothendieckGroup(M::ModGrOld, G::Grp, groupgens::SeqEnum) -> ModTupRngElt

    Parameters:
    	M - a graded module
    	G - a group
    	groupgens - a sequence

    Options:
    	UseCharacters - true or false

    Description:
        Suppose that +M+ is a graded module such that the action matrices of +M+ indexed by +groupgens+ define a +G+-module. This means that the algebra over which +M+ is defined contains the group algebra of +G+ and we restrict +M+ to the group algebra. This yields a graded +G+-module and this intrinsic returns the decomposition of this module in the graded Grothendieck group.

*/
{}
    P<t> := PolynomialRing(Integers());
    p := Characteristic(BaseRing(M));
    Modules(~G,p);
    HomogeneousData(~M);
    V := RModule(P, #G`Modules[p]);

    dec := Zero(V);
    for d in SequenceToSet(M`RowDegrees) do
        Md := GModule(G, [ Matrix(M`HomogeneousComponentMatrices[<i,d>]) : i in groupgens]);
        ddec := DecompositionInGrothendieckGroup(Md : UseCharacters:=UseCharacters);
        for j:=1 to #G`Modules[p] do
            dec[j] +:= ddec[j]*t^d;
        end for;
    end for;

    return dec;

end intrinsic;

//=============================================================================
/*
    Namespace: ModGrOld

    Additions to the category +ModGrOld+.
*/
declare type ModGrOld[ModGrOldElt];

declare attributes ModGrOld:
    Rep,
    /*
    	Attribute: Rep

    	Representation type of the matrices (Dense or Sparse)
    */
    Matrices,
    /*
    	Attribute: Matrices

    	The action matrices
    */
    MatrixDegrees,
    /*
    	Attribute: MatrixDegrees

    	Degrees of the generators
    */
    RowDegrees,
    /*
    	Attribute: RowDegrees

    	Degrees of the rows of the action matrices
    */
    Dimension,
    /*
    	Attribute: Dimension

    	Dimension of the module
    */
    BaseRing,
    /*
    	Attribute: BaseRing

    	Base ring of the module
    */
    VectorSpace,
    /*
    	Attribute: VectorSpace

    	The underlying vector space of the module
    */
    RModule,
    /*
    	Attribute: RModule

    	The corresponding (ungraded) RModule
    */
    HomogeneousComponentMatrices,
    /*
    	Attribute: HomogeneousComponentMatrices

    	Matrices describing the actions on homogeneous components
    */
    HomogeneousComponents,
    /*
    	Attribute: HomogeneousComponents

    	Vector spaces of the homogeneous components
    */
    HomogeneousComponentProjections,
    /*
    	Attribute: HomogeneousComponentProjections

    	Projections to the homogeneous components
    */
    HomogeneousComponentEmbeddings;
    /*
    	Attribute: HomogeneousComponentEmbeddings

    	Embeddings of the homogeneous components
    */


/*
    Namespace: ModGrOldElt

    Additions to the category +ModGrOldElt+.
*/
declare attributes ModGrOldElt:
    Parent,
    /*
    	Attribute: Parent

    	The parent module of a vector
    */
    Degree,
    /*
    	Attribute: Degree

    	The degree of a vector
    */
    VectorSpaceElement;
    /*
    	Attribute: VectorSpaceElement

    	The corresponding vector space element
    */

declare verbose ModGrOld, 5;
