/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Intrinsics for rational Cherednik algebras.
  The implementation of the multiplication is in AlgCheMult.
*/

//========================================================================

declare type AlgChe[AlgCheElt];

declare attributes AlgChe:
    BaseRing,	//base ring of the Cherednik algebras
    Group,		//the reflection group of the Cherednik algebra
    cParameter,	//the c-Parameter
    tParameter,	//the t-Parameter
    GroupDimension,	//the dimension of the vector space on which the group acts
    NumberOfGroupGenerators,	//number of generators of the reflection group
    NumberOfGenerators,	//number of generators of the Cherednik algebra
    Generators,	//the generators of the Cherednik algebra
    GeneratorDegrees,	//the degrees of the generators of the Cherednik algebra
    GroupAlgebra,	//the underlying group algebra of the Cherednik algebra H (H is isomorphic to it as a module)
    yxAlgebra,	//k[V \oplus V^*]
    xAlgebra,	//the algebra R[V]
    yAlgebra,   //the algebra R[V^*] (but actually we take it to be R[V] since we can work with dual action; this distinction becomes more relevant for AlgCheRes).
    xEmbedding,	//the embedding R[V] -> R[V^* + V]
    yEmbedding, //the embedding of R[V^*] into R[V^* + V] with image R[V^*]
    ySeq, //the sequence [1..n]
    xSeq, //the sequence [n+1..2*n]
    UseProductTable,	//if ProductTable should be used
    ProductTable,	//saves products of monomials in x with monomials in y
    Poisson,	//if the additional Cherednik algebra for Poisson bracket computation should be attached
    PoissonAlgebra,	//the additional Cherednik algebra for Poisson bracket computation
    PoissonBaseRingEmbedding,	//the embedding from the t=0 basering to the t=t/t^2 basering
    PoissonBaseRingProjection, //the projection from the t/t^2 basering to the t=0 basering
    PoissonModtSquare,	//Whether to compute Poisson brackets mod t^2 or not
	UseCommutatorsTable,
	CommutatorsTable, // a nested array such that [i][k][mu] (note the order!) is [y_i, x^\mu]_{s_k} in H_c. Only used if UseTable is true. The numbering of the reflections is as in W`ReflectionLibrary and W`Reflections.
	ZeroSequenceOfLengthDim, //this will simply carry the sequence [ 0 : r in [1..Dim(W)]] to avoid repeated creation of this function in CherednikCommutator
	DBDir,
	Zero;

declare attributes AlgCheElt:
    Parent,
    Element;

//=========================================================================
intrinsic RationalCherednikAlgebra(G::GrpMat, param::Tup : UseProductTable:=true, Poisson:=true, UseCommutatorsTable:=true, PoissonModtSquare:=true, GroupAlgebraRep:="") -> AlgChe
{}
    H := New(AlgChe);
    DualGroup(~G);
    SymplecticDoubling(~G);
    CherednikCoefficients(~G);

    H`Group := G;
    H`tParameter := param[1];
    H`cParameter := param[2];
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

	//the algebra R[V^* + V] where R is the base ring of H
    H`yxAlgebra := PolynomialRing(H`BaseRing, 2*H`GroupDimension);
    AssignNames(~H`yxAlgebra, (["y"*Sprint(i): i in [1..H`GroupDimension]]) cat (["x"*Sprint(i): i in [1..H`GroupDimension]]));
    Generators(~H`yxAlgebra);
    H`yxAlgebra`Zero := Zero(H`yxAlgebra);
    H`yxAlgebra`One := One(H`yxAlgebra);

    //the algebra R[V]
    H`xAlgebra := PolynomialRing(H`BaseRing, H`GroupDimension);
    AssignNames(~H`xAlgebra, ["x"*Sprint(i): i in [1..H`GroupDimension]]);
    Generators(~H`xAlgebra);
    H`xAlgebra`Zero := Zero(H`xAlgebra);
    H`xAlgebra`One := One(H`xAlgebra);

   	H`yAlgebra := PolynomialRing(H`BaseRing, H`GroupDimension);
    AssignNames(~H`yAlgebra, ["y"*Sprint(i): i in [1..H`GroupDimension]]);
    Generators(~H`yAlgebra);
    H`yAlgebra`Zero := Zero(H`yAlgebra);
    H`yAlgebra`One := One(H`yAlgebra);

    H`yxAlgebra`xPart := H`xAlgebra;
    H`yxAlgebra`yPart := H`yAlgebra;

    H`yxAlgebra`xEmbedding := hom<H`xAlgebra -> H`yxAlgebra | [H`yxAlgebra.i : i in [H`GroupDimension+1..2*H`GroupDimension]]>;
    H`yxAlgebra`yEmbedding := hom<H`yAlgebra -> H`yxAlgebra | [H`yxAlgebra.i : i in [1..H`GroupDimension]]>;

    H`GroupAlgebra := GroupAlgebra(H`yxAlgebra, H`Group : Rep:="Terms");

    H`Generators := [ H.i : i in [1..H`NumberOfGroupGenerators + 2*H`GroupDimension]];

    H`ySeq := [1..H`GroupDimension];
    H`xSeq := [H`GroupDimension+1..2*H`GroupDimension];

	H`UseProductTable := UseProductTable;
	if UseProductTable then
    	H`ProductTable := [ AssociativeArray({[1]}) : i in [1..H`GroupDimension] ];
		//initialize codomain
		for i:=1 to H`GroupDimension do
			H`ProductTable[i][[0 : j in [1..H`GroupDimension]]] := H.(Ngens(H`Group)+i); //this is x^0*y_i = 1*y_i = y_i
		end for;
	end if;

	//attach Poisson algebra
    H`Poisson := Poisson;
    if Poisson and H`tParameter eq 0 then
        R := Codomain(H`cParameter);
        S := PolynomialRing(R);
        AssignNames(~S, ["t"]);
        H`PoissonModtSquare := PoissonModtSquare;
        if PoissonModtSquare eq true then
       		T := quo<S|S.1^2>;
        	H`PoissonBaseRingEmbedding := func< r | T!r >;
        	H`PoissonBaseRingProjection := hom<T->R | [One(R)]>; //extract t-part
        	cPoisson := map<Domain(H`cParameter)->T|[<s,H`PoissonBaseRingEmbedding(H`cParameter(s))> : s in Domain(H`cParameter)]>;
        	H`PoissonAlgebra := RationalCherednikAlgebra(G,<T.1,cPoisson>);
        else
        	H`PoissonBaseRingEmbedding := hom<R->S | [S!(R.i) : i in [1..Ngens(R)]] >;
        	H`PoissonBaseRingProjection := func< f | R!MonomialCoefficient(f,S.1)>;
        	cPoisson := map<Domain(H`cParameter)->S|[<s,H`PoissonBaseRingEmbedding(H`cParameter(s))> : s in Domain(H`cParameter)]>;
        	H`PoissonAlgebra := RationalCherednikAlgebra(G,<S.1,cPoisson>);
        end if;
    end if;

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

	//attach zero element
	H`Zero := New(AlgCheElt);
    H`Zero`Parent := H;
    H`Zero`Element := Zero(H`GroupAlgebra);

    return H;

end intrinsic;

//=========================================================================
intrinsic RationalCherednikAlgebra(G::GrpMat, c::Map : UseProductTable:=true, UseCommutatorsTable:=true, Poisson:=true, PoissonModtSquare:=true, GroupAlgebraRep:="") -> AlgChe
{}

    return RationalCherednikAlgebra(G, <Zero(Codomain(c)),c> :   UseProductTable:=UseProductTable, UseCommutatorsTable:=UseCommutatorsTable, Poisson:=Poisson, PoissonModtSquare:=PoissonModtSquare);

end intrinsic;

//=========================================================================
intrinsic RationalCherednikAlgebra(G::GrpMat : Type:="GGOR", Rational:=false,  UseProductTable:=true, UseCommutatorsTable:=true, PoissonModtSquare:=true, GroupAlgebraRep:="") -> AlgChe
/*
    History:
        Sunday, October 27, 2013 16:51:05: Initial.
*/
{}
    param := FullCherednikParameter(G:Type:=Type,Rational:=Rational);
    return RationalCherednikAlgebra(G, param:  UseProductTable:=UseProductTable, UseCommutatorsTable:=UseCommutatorsTable, PoissonModtSquare:=PoissonModtSquare);

end intrinsic;

//=========================================================================
intrinsic RationalCherednikAlgebra(G::GrpMat, t::RngElt : Type:="GGOR", Rational:=false, UseProductTable:=true,UseCommutatorsTable:=true, Poisson:=true, PoissonModtSquare:=true, GroupAlgebraRep:="") -> AlgChe
/*
    History:
        Sunday, October 27, 2013 16:51:05: Initial.
*/
{}
    c:=CherednikParameter(G: Type:=Type, Rational:=Rational);
    R:=Codomain(c);
    H := RationalCherednikAlgebra(G, <R!t,c> :    UseProductTable:=UseProductTable,UseCommutatorsTable:=UseCommutatorsTable,Poisson:=Poisson,PoissonModtSquare:=PoissonModtSquare);
    if assigned G`DBDir then
    	H`DBDir := G`DBDir*"Cherednik/"*Type*"/Generic";
    end if;
    return H;

end intrinsic;


//=========================================================================
intrinsic 'eq'(H1::AlgChe, H2::AlgChe) -> BoolElt
{}

    return H1`Group eq H2`Group and H1`BaseRing eq H2`BaseRing and H1`tParameter eq H2`tParameter and H1`cParameter eq H2`cParameter;

end intrinsic;

//=========================================================================
intrinsic 'eq'(h1::AlgCheElt, h2::AlgCheElt) -> BoolElt
{}

    return h1`Element eq h2`Element;

end intrinsic;

//=========================================================================
intrinsic Ngens(H::AlgChe) -> RngIntElt
{}

    return H`NumberOfGenerators;

end intrinsic;

//=========================================================================
intrinsic One(H::AlgChe) -> AlgCheElt
{}

    one := New(AlgCheElt);
    one`Parent := H;
    one`Element := One(H`GroupAlgebra);
    return one;

end intrinsic;

//=========================================================================
intrinsic Zero(H::AlgChe) -> AlgCheElt
{}

    zero := New(AlgCheElt);
    zero`Parent := H;
    zero`Element := Zero(H`GroupAlgebra);
    return zero;

end intrinsic;

//=========================================================================
intrinsic '.'(H::AlgChe, i::RngIntElt) -> AlgCheElt
{}

    h := New(AlgCheElt);
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
intrinsic '+'(h1::AlgCheElt, h2::AlgCheElt) -> AlgCheElt
{}

    sum := New(AlgCheElt);
    sum`Parent := h1`Parent;
    sum`Element := h1`Element + h2`Element;

    return sum;

end intrinsic;

//=========================================================================
intrinsic '+'(a::RngElt, h::AlgCheElt) -> AlgCheElt
{}

    res := New(AlgCheElt);
    res`Parent := h`Parent;
    res`Element := a + h`Element;

    return res;

end intrinsic;

//=========================================================================
intrinsic '+'(h::AlgCheElt, a::RngElt) -> AlgCheElt
{}

    res := New(AlgCheElt);
    res`Parent := h`Parent;
    res`Element := h`Element + a;

    return res;

end intrinsic;

//=========================================================================
intrinsic '-'(a::RngElt, h::AlgCheElt) -> AlgCheElt
{}

    res := New(AlgCheElt);
    res`Parent := h`Parent;
    res`Element := a - h`Element;

    return res;

end intrinsic;

//=========================================================================
intrinsic '-'(h::AlgCheElt, a::RngElt) -> AlgCheElt
{}

    res := New(AlgCheElt);
    res`Parent := h`Parent;
    res`Element := h`Element - a;

    return res;

end intrinsic;

//=========================================================================
intrinsic '-'(h1::AlgCheElt, h2::AlgCheElt) -> AlgCheElt
{}

    sum := New(AlgCheElt);
    sum`Parent := h1`Parent;
    sum`Element := h1`Element - h2`Element;

    return sum;

end intrinsic;

//=========================================================================
intrinsic '-'(h::AlgCheElt) -> AlgCheElt
{}

    res := New(AlgCheElt);
    res`Parent := h`Parent;
    res`Element := -1*h`Element;

    return res;

end intrinsic;

//=========================================================================
intrinsic '*'(a::RngElt, h::AlgCheElt) -> AlgCheElt
{}

    ah := New(AlgCheElt);
    ah`Parent := h`Parent;
    ah`Element := a*h`Element;

    return ah;

end intrinsic;

//=========================================================================
intrinsic '*'(h::AlgCheElt, a::RngElt) -> AlgCheElt
{}

    ah := New(AlgCheElt);
    ah`Parent := h`Parent;
    ah`Element := h`Element*a;

    return ah;

end intrinsic;

//=========================================================================
intrinsic '^'(h::AlgCheElt, n::RngIntElt) -> AlgCheElt
{}

    pow := One(h`Parent);
    for i:=1 to n do
        pow := pow*h;
    end for;

    return pow;

end intrinsic;



//=========================================================================
intrinsic Terms(h::AlgCheElt) -> List
{}
    H := h`Parent;
    terms := [**];
    for g in Support(h`Element) do
        t := New(AlgCheElt);
        t`Parent := H;
        t`Element := Coefficient(h`Element, g)*H`GroupAlgebra!g;
        Append(~terms, t);
    end for;

    return terms;

end intrinsic;


//=========================================================================
intrinsic IsCoercible(H::AlgChe, g::GrpMatElt) -> BoolElt
{}
    if H`Group ne Parent(g) then
        return false;
    end if;

    gH := New(AlgCheElt);
    gH`Parent := H;
    gH`Element := H`GroupAlgebra!g;

    return true,gH;

end intrinsic;

//=========================================================================
intrinsic IsCoercible(H::AlgChe, f::RngMPolElt) -> BoolElt
{
	Coercion of an element of K[V+V^*].
}

    fH := New(AlgCheElt);
    fH`Parent := H;

    phi := hom<Parent(f)->H`yxAlgebra | [H`yxAlgebra.i : i in [1..Ngens(Parent(f))]]>;

    fH`Element := phi(f)*One(H`GroupAlgebra);

    return true, fH;

end intrinsic;


//=========================================================================
intrinsic EulerElement(H::AlgChe) -> AlgCheElt
/*
    As in [BR13], 4.4
*/
{
	Returns the Euler element of H.
}
    d := Dimension(H`Group);

    eu := Zero(H);

    eu`Element := eu`Element + H`GroupAlgebra!H`yxAlgebra!(d*H`tParameter);

    for i:=1 to d do
        eu`Element := eu`Element+H`GroupAlgebra!(H`yxAlgebra.(i+d)*H`yxAlgebra.i);
    end for;

    for s in H`Group`ReflectionLibraryFlat do
        eu`Element := eu`Element+H`GroupAlgebra!(H`yxAlgebra!(H`BaseRing!(Determinant(s`Element)/(Determinant(s`Element)-1))*H`cParameter(s`ReflectionClass)))*s`Element;
    end for;

    return eu;

end intrinsic;

//=========================================================================
intrinsic IsCentral(h::AlgCheElt) -> BoolElt
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
intrinsic GroupAlgebraPart(h::AlgCheElt) -> AlgGrpElt
{The map Omega from Bonnafe-Rouquier.}

	H := h`Parent;
	W := H`Group;
	R := GroupAlgebra(H`BaseRing, W);
	grppart := Zero(R);

	for g in Support(h`Element) do
		coeff := Coefficient(h`Element, g);
		coeffkilled := Evaluate(coeff, [ Zero(H`BaseRing) : i in [1..2*H`GroupDimension]]);
		grppart +:= (H`BaseRing!coeffkilled)*R!g;
	end for;

	return grppart;

end intrinsic;

//============================================================================
intrinsic IsGeneric(H::AlgChe) -> BoolElt, MonStgElt
{
	True if H is the generic rational Cherednik algebra at t=0. If so, it also returns the parameter type.
}

	R := H`BaseRing;
	W := H`Group;
	if Type(R) eq RngMPol and Ngens(R) eq #W`ReflectionClasses then
		//determine type
		c := CherednikParameter(W:Type:="EG", Rational:=false);
		if [H`cParameter(i) : i in [1..#W`ReflectionClasses]] eq [ H`BaseRing!c(i) : i in [1..#W`ReflectionClasses] ] then
			return true, "EG";
		end if;
		c := CherednikParameter(W:Type:="GGOR", Rational:=false);
		if [H`cParameter(i) : i in [1..#W`ReflectionClasses]] eq [ H`BaseRing!c(i) : i in [1..#W`ReflectionClasses] ] then
			return true, "GGOR";
		end if;
		return true,_; //cannot determine parameter type
	else
		return false,_;
	end if;

end intrinsic;

//============================================================================
intrinsic IsCoercible(H::AlgChe, h::AlgCheElt) -> BoolElt, AlgCheElt
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
	tparam := H`BaseRing!h`Parent`tParameter;
	W := H`Group;
	cparam := [ H`BaseRing!h`Parent`cParameter(i) : i in [1..#W`ReflectionClasses] ];
	if not (tparam eq H`tParameter and cparam eq [ H`cParameter(i) : i in [1..#W`ReflectionClasses]]) then
		return false, _;
	end if;
	newh := Zero(H);
	for w in Support(h`Element) do
		newh +:= (H`yxAlgebra!Coefficient(h`Element,w))*H`GroupAlgebra!w;
	end for;

	return true, newh;

end intrinsic;


//=========================================================================
intrinsic Print(H::AlgChe)
/*
    History:
        * Thursday, April 17, 2014 at 13:41:41: Initial.
*/
{}

    printf "Rational Cherednik algebra\n", H`BaseRing;
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
    printf "t-parameter:\n";
    IndentPush();
    printf "%o\n", H`tParameter;
    IndentPop();
    printf "c-parameter:\n";
    IndentPush();
    printf "%o", H`cParameter;
    IndentPop();

end intrinsic;

//=========================================================================
intrinsic Parent(x::AlgCheElt) -> AlgChe
{}

    return x`Parent;

end intrinsic;

//=========================================================================
intrinsic Print(x::AlgCheElt)
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
intrinsic Rprint(h::AlgCheElt) -> MonStgElt
{
	Reversible print of h.
}

	SetColumns(0); //idiotic!

	W := h`Parent`Group;
	K := BaseRing(W);
	H := h`Parent;

	str := "/*\n\tCode for a Cherednik algebra element\n";
	str *:= "\tVersion: "*CHAMP_GetVersion()*"\n";
	str *:= "\tDate: "*Date()*"\n*/\n";
	str *:= "//base ring of the group\n";
	str *:= "K := ";
	fieldstr := "";
	if Type(K) eq FldRat then
		fieldstr := "RationalField()";
		str *:= fieldstr*";\n";
	elif Type(K) eq FldCyc then
		fieldstr := "CyclotomicField("*Sprint(CyclotomicOrder(K))*")";
		str *:= fieldstr*";\n";
		str *:= Sprint(K.1)*" := RootOfUnity("*Sprint(CyclotomicOrder(K))*");\n";
	else
		error "Not yet implemented for this type of base ring.";
	end if;
	str *:= "//the group\n";
	str *:= "W := MatrixGroup<"*Sprint(Dimension(W))*", "*fieldstr*" | ";
	for i:=1 to Ngens(W) do
		str *:= Replace(Sprint(Eltseq(W.i)), "\n", "");
		if i lt Ngens(W) then
			str *:= " , ";
		end if;
	end for;
	str *:= ">;\n";

	str *:= "//the parameters of the Cherednik algebra\n";
	R := H`BaseRing;
	if Type(R) eq RngMPol then
		if BaseRing(R) ne K then
			error "Not yet implemented for this type of base ring.";
		end if;
		rstr := "PolynomialRing("*fieldstr*", "*Sprint(Ngens(R))*")";
		str *:= "R := "*rstr*";\n";
		for i:=1 to Ngens(R) do
			str *:= Sprint(R.i)*" := R."*Sprint(i)*";\n";
		end for;
	elif Type(R) eq Fld then
		if R ne K then
			error "Not yet implemented for this type of base ring.";
		end if;
		rstr := fieldstr;
		str *:= "R := "*rstr*";\n";
	end if;
	str *:= "t := "*Sprint(H`tParameter)*";\n";
	str *:= "c := map<"*Sprint(Domain(H`cParameter))*"-> R |[";
	for i:=1 to #W`ReflectionClasses do
		str *:= "<"*Sprint(i)*","*Sprint(H`cParameter(i))*">";
		if i lt #W`ReflectionClasses then
			str *:= ",";
		end if;
	end for;
	str *:= "]>;\n";

	str *:= "//the Cherednik algebra\n";
	str *:= "H := RationalCherednikAlgebra(W,<t,c>);\n";
	str *:= "A := H`GroupAlgebra;\n";
	for i:=1 to 2*Dimension(W) do
		str *:= Sprint(H`yxAlgebra.i)*" := H`yxAlgebra."*Sprint(i)*";\n";
	end for;

	str *:= "//the support of the Cherednik algebra element (as words in the generators)\n";
	str *:= "supp := [";
	supp := SetToSequence(Support(h`Element));
	for i:=1 to #supp do
		w := supp[i];
		str *:= Sprint(ElementToWord(w : Method:="FPGroup"));
		if i lt #supp then
			str *:= ",";
		end if;
	end for;
	str *:= "];\n";

	str *:= "//the coefficients of the Cherednik algebra element (these are elements of S)\n";
	for i:=1 to #supp do
		str *:= "coeff"*Sprint(i)*" := "*Sprint(Coefficient(h`Element, supp[i]))*";\n";
	end for;
	//str *:= "//construct element\n";
	str *:= "hA := Zero(A);\n";
	for i:=1 to #supp do
		str *:= "hA +:= coeff"*Sprint(i)*"*(A!WordToElement(W,supp["*Sprint(i)*"]));\n";
	end for;
	str *:= "h := Zero(H); h`Element := hA;\n";
	str *:= "return h";

	return str;

end intrinsic;

//============================================================================
intrinsic Bidegree(h::AlgCheElt) -> MonStgElt
{The bidegree of a homogeneous element h. It is not checked whether h is in fact homogeneous.}

	for w in Support(h`Element) do
		coeff := Coefficient(h`Element, w);
		return Bidegree(coeff);
	end for;

end intrinsic;
