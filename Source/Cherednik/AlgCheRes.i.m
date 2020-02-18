/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
  Intrinsics for restricted rational Cherednik algebras.
  The implementation of the multiplication is in AlgCheResMult.
*/


//========================================================================

declare type AlgCheRes[AlgCheResElt];

declare attributes AlgCheRes:
    BaseRing,	//base ring of the Cherednik algebras
    Group,		//the reflection group of the Cherednik algebra
    cParameter,	//the c-Parameter
    GroupDimension,	//the dimension of the vector space on which the group acts
    NumberOfGroupGenerators,	//number of generators of the reflection group
    NumberOfGenerators,	//number of generators of the Cherednik algebra
    Generators,	//the generators of the Cherednik algebra
    GeneratorDegrees,	//the degrees of the generators of the Cherednik algebra
    GroupAlgebra,	//the underlying group algebra of the Cherednik algebra H (H is isomorphic to it as a module)
    yxAlgebra,	//R[V]_G \otimes R[V^*]_G
    xAlgebra,	//the algebra R[V]_G
    yAlgebra,   //the algebra R[V^*]_G
    xEmbedding,	//the embedding R[V]_G -> R[V]_G \otimes R[V^*]_G
    yEmbedding, //the embedding of R[V^*]_G into R[V]_G \otimes R[V^*]_G
    ySeq, //the sequence [1..n]
    xSeq, //the sequence [n+1..2*n]
    CoinvariantAlgebraToyxAlgebra, //the embedding K[V]_G -> R[V]_G \otimes R[V^*]_G
    DualCoinvariantAlgebraToyxAlgebra, //the embedding K[V^*]_G -> R[V]_G \otimes R[V^*]_G
    CoinvariantAlgebraToxAlgebra, //the embedding K[V]_G -> R[V]_G
    DualCoinvariantAlgebraToyAlgebra, //the embedding K[V^*]_G -> R[V^*]_G
    UseProductTable,	//if ProductTable should be used
    ProductTable,	//saves products of monomials in x with monomials in y
		UseCommutatorsTable,
		CommutatorsTable, // a nested array such that [i][k][mu] (note the order!) is [y_i, x^\mu]_{s_k} in H_c. Only used if UseTable is true. The numbering of the reflections is as in W`ReflectionLibrary and W`Reflections.
		ZeroSequenceOfLengthDim, //this will simply carry the sequence [ 0 : r in [1..Dim(W)]] to avoid repeated creation of this function in CherednikCommutator
		Zero,
		Basis, //a basis of the restricted rational Cherednik algebra
		MatrixAlgebra,
		VectorSpace, //underlying vector space
		MatrixAlgebraOfDegreeZeroPart,
		qField, //rational function field K(q) for characters etc.
		Initialized;

declare attributes AlgCheResElt:
    Parent,
    Element;


//=========================================================================
intrinsic RestrictedRationalCherednikAlgebra(G::GrpMat, c::Map : UseProductTable:=true, UseCommutatorsTable:=true, Verbose:=false, Init:=false) -> AlgCheRes
{}
    H := New(AlgCheRes);
    DualGroup(~G);
    SymplecticDoubling(~G);
    CherednikCoefficients(~G);

    H`Group := G;
    H`cParameter := c;
    H`BaseRing := Codomain(H`cParameter);
    H`GroupDimension := Dimension(G);
    H`NumberOfGroupGenerators := Ngens(G);
    ReflectionLibrary(~G);

    H`GeneratorDegrees:=[];
    for i in Reverse([1..#Generators(G)]) do
        Append(~H`GeneratorDegrees, 0);
    end for;
    for i in Reverse([1..H`GroupDimension]) do
        Append(~H`GeneratorDegrees, -1);
    end for;
    for i in Reverse([1..H`GroupDimension]) do
        Append(~H`GeneratorDegrees, 1);
    end for;

    H`NumberOfGenerators := #H`GeneratorDegrees;


    H`ySeq := [1..H`GroupDimension];
    H`xSeq := [H`GroupDimension+1..2*H`GroupDimension];

		H`UseProductTable := UseProductTable;


		H`UseCommutatorsTable := UseCommutatorsTable;

		//initialize commutators table
		if H`UseCommutatorsTable then
			ReflectionLibrary(~G);
			CoordinateAlgebra(~G);
			N := #G`ReflectionLibraryFlat;
			d := G`Dimension;

			// just a safe check
			for k:=1 to N do
				assert G`ReflectionLibraryFlat[k]`ReflectionNumber eq k;
			end for;

			H`CommutatorsTable := < < AssociativeArray(Universe({[1,1]})) : k in [1..N]> : i in [1..d] >;
		end if;

		H`ZeroSequenceOfLengthDim := [ 0 : r in [1..H`GroupDimension] ];


		R<q> := RationalFunctionField(BaseRing(G));
		H`qField := R;

		if Init then
			Initialize(~H : Verbose:=Verbose);
		end if;

    return H;

end intrinsic;

//=========================================================================
intrinsic RestrictedRationalCherednikAlgebra(G::GrpMat : Type:="GGOR",   UseProductTable:=true, UseCommutatorsTable:=true, Verbose:=false, Init:=false) -> AlgCheRes
{}
    param := CherednikParameter(G:Type:=Type,Rational:=true);
    return RestrictedRationalCherednikAlgebra(G, param:  UseProductTable:=UseProductTable, UseCommutatorsTable:=UseCommutatorsTable, Verbose:=Verbose, Init:=Init);

end intrinsic;

intrinsic Initialize(~H::AlgCheRes : Verbose:=false)
{}

	if assigned H`Initialized then
		return;
	end if;

	G := H`Group;

	//construct R[V]/m0, R[V^*]/m0 and R[V^*+V]/m0, where m0 is the respective origin
	//we do a trick here by assigning different weights to the y's and x's. this makes the Groebner basis stuff much more efficient
	H`yxAlgebra := ChangeOrder(PolynomialRing(H`BaseRing, 2*H`GroupDimension), <"grevlexw", [1: i in [1..H`GroupDimension]] cat [2 : i in [1..H`GroupDimension]]>);
	AssignNames(~H`yxAlgebra, (["y"*Sprint(i): i in [1..H`GroupDimension]]) cat (["x"*Sprint(i): i in [1..H`GroupDimension]]));
	Generators(~H`yxAlgebra);
	H`yxAlgebra`Zero := Zero(H`yxAlgebra);
	H`yxAlgebra`One := One(H`yxAlgebra);

	CoinvariantAlgebra(~G : Verbose:=Verbose);
	DualCoinvariantAlgebra(~G : Verbose:=Verbose);

	phi := hom<G`CoordinateAlgebra -> H`yxAlgebra | [ H`yxAlgebra.(H`GroupDimension+i) : i in [1..H`GroupDimension]]>;
	psi := hom<G`DualCoordinateAlgebra -> H`yxAlgebra | [ H`yxAlgebra.i : i in [1..H`GroupDimension]]>; //for coercion of hilbert ideal
	I := ideal<H`yxAlgebra|[psi(f) : f in Basis(G`DualHilbertIdeal)] cat [phi(f) : f in Basis(G`HilbertIdeal)]>;
	H`yxAlgebra := H`yxAlgebra/I;
	Generators(~H`yxAlgebra);
	H`yxAlgebra`Zero := Zero(H`yxAlgebra);
	H`yxAlgebra`One := One(H`yxAlgebra);

	H`xAlgebra := PolynomialRing(H`BaseRing, H`GroupDimension);
	AssignNames(~H`xAlgebra, ["x"*Sprint(i): i in [1..H`GroupDimension]]);
	Generators(~H`xAlgebra);
	H`xAlgebra`Zero := Zero(H`xAlgebra);
	H`xAlgebra`One := One(H`xAlgebra);
	phi := hom<G`CoordinateAlgebra -> H`xAlgebra | [ H`xAlgebra.i : i in [1..H`GroupDimension]]>;
	I := ideal<H`xAlgebra | [phi(f) : f in Basis(G`HilbertIdeal)]>;
	H`xAlgebra := H`xAlgebra/I;
	Generators(~H`xAlgebra);
	H`xAlgebra`Zero := Zero(H`xAlgebra);
	H`xAlgebra`One := One(H`xAlgebra);

	H`yAlgebra := PolynomialRing(H`BaseRing, H`GroupDimension);
	AssignNames(~H`yAlgebra, ["y"*Sprint(i): i in [1..H`GroupDimension]]);
	Generators(~H`yAlgebra);
	H`yAlgebra`Zero := Zero(H`yAlgebra);
	H`yAlgebra`One := One(H`yAlgebra);
	psi := hom<G`DualGroup`CoordinateAlgebra -> H`yAlgebra | [ H`yAlgebra.i : i in [1..H`GroupDimension]]>;
	I := ideal<H`yAlgebra | [psi(f) : f in Basis(G`DualHilbertIdeal)]>;
	H`yAlgebra := H`yAlgebra/I;
	Generators(~H`yAlgebra);
	H`yAlgebra`Zero := Zero(H`yAlgebra);
	H`yAlgebra`One := One(H`yAlgebra);

	H`yxAlgebra`xPart := H`xAlgebra;
	H`yxAlgebra`yPart := H`yAlgebra;

	H`yxAlgebra`xEmbedding := hom<H`xAlgebra -> H`yxAlgebra | [H`yxAlgebra.i : i in [H`GroupDimension+1..2*H`GroupDimension]]>;
	H`yxAlgebra`yEmbedding := hom<H`yAlgebra -> H`yxAlgebra | [H`yxAlgebra.i : i in [1..H`GroupDimension]]>;

	H`CoinvariantAlgebraToyxAlgebra := hom<G`CoinvariantAlgebra-> H`yxAlgebra | [H`yxAlgebra`Generators[i] : i in [H`GroupDimension+1..2*H`GroupDimension]]>;
	H`DualCoinvariantAlgebraToyxAlgebra := hom<G`DualCoinvariantAlgebra-> H`yxAlgebra | [H`yxAlgebra`Generators[i] : i in [1..H`GroupDimension]]>;

	//set bases
	Basis(~G`CoinvariantAlgebra : SortByDegree:=true);
	Basis(~G`DualCoinvariantAlgebra : SortByDegree:=true);
	phi := hom<G`CoinvariantAlgebra -> H`xAlgebra | [ H`xAlgebra.i : i in [1..H`GroupDimension]]>;
	H`CoinvariantAlgebraToxAlgebra := phi;
	H`xAlgebra`Basis := {@ phi(G`CoinvariantAlgebra`Basis[i]) : i in [1..#G`CoinvariantAlgebra`Basis] @};
	SetBasis(~H`xAlgebra, H`xAlgebra`Basis);
	phi := hom<G`DualCoinvariantAlgebra -> H`yAlgebra | [ H`yAlgebra.i : i in [1..H`GroupDimension]]>;
	H`DualCoinvariantAlgebraToyAlgebra := phi;
	H`yAlgebra`Basis := {@ phi(G`DualCoinvariantAlgebra`Basis[i]) : i in [1..#G`DualCoinvariantAlgebra`Basis] @};
	SetBasis(~H`yAlgebra, H`yAlgebra`Basis);

	H`GroupAlgebra := GroupAlgebra(H`yxAlgebra, H`Group : Rep:="Terms");

	//attach zero element
	H`Zero := New(AlgCheResElt);
	H`Zero`Parent := H;
	H`Zero`Element := Zero(H`GroupAlgebra);

	if H`UseProductTable then
			H`ProductTable := [ AssociativeArray({[1]}) : i in [1..H`GroupDimension] ];
		//initialize codomain
		for i:=1 to H`GroupDimension do
			H`ProductTable[i][[0 : j in [1..H`GroupDimension]]] := H.(Ngens(H`Group)+i); //this is x^0*y_i = 1*y_i = y_i
		end for;
	end if;

	H`Generators := [ H.i : i in [1..H`NumberOfGroupGenerators + 2*H`GroupDimension]];

	H`Initialized := true;

end intrinsic;

//=========================================================================
intrinsic 'eq'(H1::AlgCheRes, H2::AlgCheRes) -> BoolElt
{}

    return H1`Group eq H2`Group and H1`BaseRing eq H2`BaseRing and H1`cParameter eq H2`cParameter;

end intrinsic;

//=========================================================================
intrinsic 'eq'(h1::AlgCheResElt, h2::AlgCheResElt) -> BoolElt
{}

    return h1`Element eq h2`Element;

end intrinsic;

//=========================================================================
intrinsic Ngens(H::AlgCheRes) -> RngIntElt
{}

    return H`NumberOfGenerators;

end intrinsic;

//=========================================================================
intrinsic One(H::AlgCheRes) -> AlgCheResElt
{}

    one := New(AlgCheResElt);
    one`Parent := H;
    one`Element := One(H`GroupAlgebra);
    return one;

end intrinsic;

//=========================================================================
intrinsic Zero(H::AlgCheRes) -> AlgCheResElt
{}

    zero := New(AlgCheResElt);
    zero`Parent := H;
    zero`Element := Zero(H`GroupAlgebra);
    return zero;

end intrinsic;

//=========================================================================
intrinsic '.'(H::AlgCheRes, i::RngIntElt) -> AlgCheResElt
{}

    h := New(AlgCheResElt);
    h`Parent := H;
    G := H`Group;

    if i in [1..H`NumberOfGroupGenerators] then
        h`Element := H`GroupAlgebra!(H`Group.(i));
    elif i in [H`NumberOfGroupGenerators+1..H`NumberOfGroupGenerators+H`GroupDimension] then
        h`Element := H`GroupAlgebra!(H`yxAlgebra.(i-H`NumberOfGroupGenerators));
    elif i in [H`GroupDimension+H`NumberOfGroupGenerators..2*H`GroupDimension+H`NumberOfGroupGenerators] then
        h`Element := H`GroupAlgebra!(H`yxAlgebra.(i-H`NumberOfGroupGenerators));
    end if;

    return h;

end intrinsic;

//=========================================================================
intrinsic '+'(h1::AlgCheResElt, h2::AlgCheResElt) -> AlgCheResElt
{}

    sum := New(AlgCheResElt);
    sum`Parent := h1`Parent;
    sum`Element := h1`Element + h2`Element;

    return sum;

end intrinsic;

//=========================================================================
intrinsic '+'(a::RngElt, h::AlgCheResElt) -> AlgCheResElt
{}

    res := New(AlgCheResElt);
    res`Parent := h`Parent;
    res`Element := a + h`Element;

    return res;

end intrinsic;

//=========================================================================
intrinsic '+'(h::AlgCheResElt, a::RngElt) -> AlgCheResElt
{}

    res := New(AlgCheResElt);
    res`Parent := h`Parent;
    res`Element := h`Element + a;

    return res;

end intrinsic;

//=========================================================================
intrinsic '-'(a::RngElt, h::AlgCheResElt) -> AlgCheResElt
{}

    res := New(AlgCheResElt);
    res`Parent := h`Parent;
    res`Element := a - h`Element;

    return res;

end intrinsic;

//=========================================================================
intrinsic '-'(h::AlgCheResElt, a::RngElt) -> AlgCheResElt
{}

    res := New(AlgCheResElt);
    res`Parent := h`Parent;
    res`Element := h`Element - a;

    return res;

end intrinsic;

//=========================================================================
intrinsic '-'(h1::AlgCheResElt, h2::AlgCheResElt) -> AlgCheResElt
{}

    sum := New(AlgCheResElt);
    sum`Parent := h1`Parent;
    sum`Element := h1`Element - h2`Element;

    return sum;

end intrinsic;

//=========================================================================
intrinsic '-'(h::AlgCheResElt) -> AlgCheResElt
{}

    res := New(AlgCheResElt);
    res`Parent := h`Parent;
    res`Element := -1*h`Element;

    return res;

end intrinsic;

//=========================================================================
intrinsic '*'(a::RngElt, h::AlgCheResElt) -> AlgCheResElt
{}

    ah := New(AlgCheResElt);
    ah`Parent := h`Parent;
    ah`Element := a*h`Element;

    return ah;

end intrinsic;

//=========================================================================
intrinsic '*'(h::AlgCheResElt, a::RngElt) -> AlgCheResElt
{}

    ah := New(AlgCheResElt);
    ah`Parent := h`Parent;
    ah`Element := h`Element*a;

    return ah;

end intrinsic;

//=========================================================================
intrinsic '^'(h::AlgCheResElt, n::RngIntElt) -> AlgCheResElt
{}

    pow := One(h`Parent);
    for i:=1 to n do
        pow := pow*h;
    end for;

    return pow;

end intrinsic;



//=========================================================================
intrinsic Terms(h::AlgCheResElt) -> List
{}
    H := h`Parent;
    terms := [**];
    for g in Support(h`Element) do
        t := New(AlgCheResElt);
        t`Parent := H;
        t`Element := Coefficient(h`Element, g)*H`GroupAlgebra!g;
        Append(~terms, t);
    end for;

    return terms;

end intrinsic;


//=========================================================================
intrinsic IsCoercible(H::AlgCheRes, g::GrpMatElt) -> BoolElt
{}
    if H`Group ne Parent(g) then
        return false;
    end if;

    gH := New(AlgCheResElt);
    gH`Parent := H;
    gH`Element := H`GroupAlgebra!g;

    return true,gH;

end intrinsic;

//=========================================================================
intrinsic IsCoercible(H::AlgCheRes, f::RngMPolResElt) -> BoolElt
{
	Coercion of an element of K[V+V^*].
}

    fH := New(AlgCheResElt);
    fH`Parent := H;

    phi := hom<Parent(f)->H`yxAlgebra | [H`yxAlgebra.i : i in [1..Ngens(Parent(f))]]>;

    fH`Element := phi(f)*One(H`GroupAlgebra);

    return true, fH;

end intrinsic;


//=========================================================================
intrinsic IsCentral(h::AlgCheResElt) -> BoolElt
{True iff h is central in the Cherednik algebra.}
    H := Parent(h);

    for i:=1 to #H`Generators do
        if H.i*h ne h*H.i then
            return false;
        end if;
    end for;
    return true;

end intrinsic;

//============================================================================
intrinsic IsCoercible(H::AlgCheRes, h::AlgCheResElt) -> BoolElt, AlgCheResElt
{}

	if H`Group ne h`Parent`Group then
		return false, _;
	end if;
	if H`Group`ReflectionClasses ne h`Parent`Group`ReflectionClasses then
		return false, _;
	end if;
	if Type(h`Parent`BaseRing) ne Type(H`BaseRing) then
		return false, _;
	end if;
	W := H`Group;
	cparam := [ H`BaseRing!h`Parent`cParameter(i) : i in [1..#W`ReflectionClasses] ];
	if not cparam eq [ H`cParameter(i) : i in [1..#W`ReflectionClasses]] then
		return false, _;
	end if;
	newh := Zero(H);
	for w in Support(h`Element) do
		newh +:= (H`yxAlgebra!Coefficient(h`Element,w))*H`GroupAlgebra!w;
	end for;

	return true, newh;

end intrinsic;


//=========================================================================
intrinsic Print(H::AlgCheRes)
/*
    History:
        * Thursday, April 17, 2014 at 13:41:41: Initial.
*/
{}

    printf "Restricted rational Cherednik algebra\n", H`BaseRing;
    printf "Generators:\n";
    IndentPush();
    str := "";
    for i:=1 to H`NumberOfGroupGenerators do
        str *:= "w"*Sprint(i)*", ";
    end for;
    for i:=1 to H`GroupDimension do
        str *:= "y"*Sprint(i)*", ";
    end for;
    for i:=1 to H`GroupDimension do
        str *:= "x"*Sprint(i);
        if i lt H`GroupDimension then
            str *:= ", ";
        end if;
    end for;
    printf "%o\n", str;
    IndentPop();

    printf "Generator degrees:\n";
    IndentPush();
    str := "";
    for i:=1 to Ngens(H) do
        str *:= Sprint(H`GeneratorDegrees[i]);
        if i lt Ngens(H) then
            str *:= ", ";
        end if;
    end for;
    printf "%o\n", str;
    IndentPop();
    printf "Base ring:\n";
    IndentPush();
    printf "%o\n", H`BaseRing;
    IndentPop();
    printf "Group:\n";
    IndentPush();
    printf "%o\n", H`Group;
    IndentPop();
    printf "c-parameter:\n";
    IndentPush();
    printf "%o", H`cParameter;
    IndentPop();

end intrinsic;

//=========================================================================
intrinsic Parent(x::AlgCheResElt) -> AlgCheRes
{}

    return x`Parent;

end intrinsic;

//=========================================================================
intrinsic Print(x::AlgCheResElt)
/*
    History:
        Friday, September 20, 2013 16:27:50: Initial.
*/
{}

    if x`Element eq 0 then
        printf "0";
    else
        supp := Support(x`Element);
        H := Parent(x);
        N := #supp;
        i:=1;
        for g in supp do
            printf "%o*(%o)", g, Coefficient(x`Element, g);
            if i lt N then
                printf "\n+\n";
            end if;
            i +:= 1;
        end for;
    end if;


end intrinsic;

//============================================================================
intrinsic Basis(~H::AlgCheRes)
{Assigns a basis of H.}

	if assigned H`Basis then
		return;
	end if;

	W := H`Group;
	n := H`GroupDimension;
	Basis(~W`CoinvariantAlgebra);
	Basis(~W`DualCoinvariantAlgebra);
	NumberingMap(~W);

	phi := H`CoinvariantAlgebraToyxAlgebra;
	psi := H`DualCoinvariantAlgebraToyxAlgebra;

	H`Basis := {@ @};
	for i:=1 to #W do
		for j:=1 to #W`DualCoinvariantAlgebra`Basis do
			for k:=1 to #W`CoinvariantAlgebra`Basis do
				b := New(AlgCheResElt);
				b`Parent := H;
				b`Element := psi(W`DualCoinvariantAlgebra`Basis[j])*phi(W`CoinvariantAlgebra`Basis[k])*H`GroupAlgebra!(W`InverseNumberingMap(i));
				H`Basis join:={@b@};
			end for;
		end for;
	end for;

	H`VectorSpace := KSpace(H`BaseRing, #H`Basis);

	//self check
	for b:=1 to #H`Basis do
		assert Support(VectorSpaceElement(H`Basis[b])) eq {b};
	end for;

end intrinsic;

//============================================================================
intrinsic VectorSpaceElement(h::AlgCheResElt) -> ModTupFldElt
{Returns the corresponding vector space element of h in the assigned Basis.}

	H := h`Parent;
	Basis(~H);
	v:=Zero(H`VectorSpace);
	n:=H`GroupDimension;
	W:=H`Group;

	for w in Support(h`Element) do
		wcoeff := Coefficient(h`Element,w); //an element of R[V^*+V]
		i := W`NumberingMap(w);
		for mon in Monomials(wcoeff) do //mon is a monomial of R[V^*+V]
			coeff := MonomialCoefficient(wcoeff,mon); //an element of R
			exp := Exponents(mon);
			monx := Monomial(W`CoinvariantAlgebra, exp[n+1..2*n]);
			mony := Monomial(W`DualCoinvariantAlgebra, exp[1..n]);
			monxvector := W`CoinvariantAlgebra`VectorSpaceMap(monx);
			monyvector := W`DualCoinvariantAlgebra`VectorSpaceMap(mony);
			for j in Support(monyvector) do
				for k in Support(monxvector) do
					//construct basis element corresponding to indices i,j,k
					b := New(AlgCheResElt);
					b`Parent := H;
					b`Element := (H`DualCoinvariantAlgebraToyxAlgebra(W`DualCoinvariantAlgebra`Basis[j])*H`CoinvariantAlgebraToyxAlgebra(W`CoinvariantAlgebra`Basis[k]))*H`GroupAlgebra!w;
					p := Position(H`Basis, b);
					if p eq 0 then
						error "Something wrong.";
					end if;
					v[p] +:= coeff;
				end for;
			end for;
		end for;
	end for;

	return v;

end intrinsic;

//============================================================================
intrinsic Basis(H::AlgCheRes) -> SetIndx
{A basis of H.}

	Basis(~H);
	return H`Basis;

end intrinsic;

//============================================================================
intrinsic StructureConstants(H::AlgCheRes) -> SetEnum
{}

	Basis(~H);
	structs := {};
	for i:=1 to #H`Basis do
		PrintPercentage(i, #H`Basis);
		for j:=i to #H`Basis do
			f := (H`Basis[i]*H`Basis[j])`Element;
			for g in Support(f) do
				c := Coefficient(f,g);
				structs join:=SequenceToSet(Coefficients(c));
			end for;
		end for;
	end for;

	return structs;

end intrinsic;

//============================================================================
intrinsic MatrixAlgebra(~H::AlgCheRes : Verbose:=true)
{}

	if assigned H`MatrixAlgebra then
		return;
	end if;

	W := H`Group;
	n := H`GroupDimension;
	Basis(~H);
	opmats := [ ZeroMatrix(H`BaseRing, #H`Basis, #H`Basis) : i in [1..#H`Generators] ];

	//compute action matrices
	for i:=1 to #H`Generators do
		h := H.i;
		for b:=1 to #H`Basis do
			bh := H`Basis[b]*h;
			opmats[i][b] := VectorSpaceElement(bh);
		end for;
	end for;

	H`MatrixAlgebra := MatrixAlgebra< H`BaseRing, #H`Basis | opmats>;

end intrinsic;

//============================================================================
intrinsic MatrixAlgebra(H::AlgCheRes : Verbose:=true) -> AlgMat
{}

	MatrixAlgebra(~H);
	return H`MatrixAlgebra;

end intrinsic;

//============================================================================
intrinsic MatrixAlgebraOfDegreeZeroPart(~H::AlgCheRes : Verbose:=true)
{}

	Basis(~H);
	deg0basis := [ i : i in [1..#H`Basis] | Bidegree(Coefficient(H`Basis[i]`Element, SetToSequence(Support(H`Basis[i]`Element))[1]))[1] eq Bidegree(Coefficient(H`Basis[i]`Element, SetToSequence(Support(H`Basis[i]`Element))[1]))[2] ]; //ugly!
	U := sub<H`VectorSpace | [H`VectorSpace.i : i in deg0basis]>;
	Uabs := KSpace(H`BaseRing, #deg0basis);

	W := H`Group;
	n := H`GroupDimension;

	//the following are algebra generators of the degree 0 part (using PBW theorem)
	deg0gens := [ H!W.i : i in [1..Ngens(W)] ] cat [ H!(H`yxAlgebra.i*H`yxAlgebra.j) : i in [1..n], j in [n+1..2*n] ];

	opmats := [ ZeroMatrix(H`BaseRing, #deg0basis, #deg0basis) : i in [1..#deg0gens] ];

	//compute action matrices
	for i:=1 to #deg0gens do
		h := deg0gens[i];
		for b:=1 to #deg0basis do
			if Verbose then
				PrintPercentage( (i-1)*#deg0basis + b, #deg0gens*#deg0basis);
			end if;
			bh := H`Basis[deg0basis[b]]*h;
			v := VectorSpaceElement(bh);
			u := Uabs!Coordinates(U,v);
			opmats[i][b] := u;
		end for;
	end for;

	H`MatrixAlgebraOfDegreeZeroPart := MatrixAlgebra< H`BaseRing, #deg0basis | opmats>;

end intrinsic;

//============================================================================
intrinsic MatrixAlgebraOfDegreeZeroPart(H::AlgCheRes : Verbose:=true) -> AlgMat
{}

	MatrixAlgebraOfDegreeZeroPart(~H);

	return H`MatrixAlgebraOfDegreeZeroPart;

end intrinsic;
