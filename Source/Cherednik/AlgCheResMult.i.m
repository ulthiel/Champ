/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Implements the multiplication in the rational Cherednik algebra.
*/


//=========================================================================
function MultiplyWithG(H, h, g)
/*
    Multiplies Cherednik algebra element h from right with group element g.
*/

    if IsIdentity(g) then
        return h;
    end if;

    prod := Zero(H);
    for i in Support(h`Element) do
        prod`Element +:= SymplecticDoublingAction(Coefficient(h`Element, i), g )*H`GroupAlgebra!i;
    end for;
    prod`Element := prod`Element*g;

    return prod;

end function;

//============================================================================
function CommutatorOfYVariableWithMonomialInX(H, i, mu, r)
/*
	The coefficient [y_i, x^\mu]_{s_r} of the reflection s_k in the commutator [y_i, x^\mu] as in [Thi15], Lemma 1.15. This is an element of R[V]_G.
*/

	W := H`Group;
	A := W`CoinvariantAlgebra;
	//A:=H`xAlgebra; //computation over W`CoinvariantAlgebra is much quicker
	comm := A`Zero;

	if IsEmpty(mu) then
		return comm;
	end if;

	//if UseCommutatorsTable, try to find commutator in the table
	if H`UseCommutatorsTable then
		if IsDefined(H`CommutatorsTable[i][r], mu) then
			return H`CommutatorsTable[i][r][mu];
		end if;
	end if;

	//otherwise, compute it
	n := H`GroupDimension;
	s := W`ReflectionLibraryFlat[r];

	//with the following code I try to prevent repeated creation of the zero sequence [ 0 : r in [1..n] ] used in this function
	leftpartexp := H`ZeroSequenceOfLengthDim;
	rightpartexp := mu;

	//now, compute the commutator
	for t:=1 to n do
		if assigned W`CherednikCoefficients then
			checoeff := W`CherednikCoefficients[i][t][s`ReflectionNumber]; //may not be defined
		else
			checoeff := CherednikCoefficient(i,t,s);
		end if;
		if IsZero(checoeff) then
			continue;
		end if;

		leftpart := Monomial(A,leftpartexp); //x_1^{\mu_1}*...*x_{t-1}^{\mu_{t-1}}

		rightpartexp[t] := 0;
		rightpart := Monomial(A,rightpartexp); //x_{t+1}^{\mu_{t+1}}*...*x_{n}^{\mu_n}

		//set leftpartexp for next round
		leftpartexp[t] := mu[t];

		rightparts := Action(rightpart, s`DualElement : Dual:=false); //will be non-zero since obtained from group action on non-zero element

		midpart := A`Zero;	//the middle part
		u := mu[t]-1; //will be the exponent in the following
		for l:=0 to u do
        	midpart +:= A`Generators[t]^l*Action(A`Generators[t]^u, s`DualElement : Dual:=false);
        	u -:= 1;
    	end for;
    	if IsZero(midpart) then
    		continue;
    	end if;
		comm +:= checoeff*leftpart*midpart*rightparts;
	end for;

	comm := H`cParameter(s`ReflectionClass)*H`CoinvariantAlgebraToxAlgebra(comm);

	if H`UseCommutatorsTable then
		H`CommutatorsTable[i][r][mu] := comm;
	end if;

	return comm;

end function;

//=========================================================================
intrinsic Commutator(H::AlgCheRes, i::RngIntElt, mu::SeqEnum, k::RngIntElt) -> RngMPolElt
{The coefficient [y_i, x^mu]_(s_k) of the reflection s_k in the commutator [y_i, x^mu].}

	return CommutatorOfYVariableWithMonomialInX(H,i,mu,k);

end intrinsic;

//=========================================================================
intrinsic CommutatorGradedVector(H::AlgCheRes, i::RngIntElt, mu::SeqEnum, k::RngIntElt) -> Tup
{The coefficient [y_i, x^mu]_(s_k) of the reflection s_k in the commutator [y_i, x^mu] as a graded vector of R[V]_G.}

	v:= H`xAlgebra`GradedVectorSpaceMap(CommutatorOfYVariableWithMonomialInX(H,i,mu,k));
	if #v eq 0 then
		return 0;
	end if;

	return v[1][2]; //v is homogeneous of degree 1 less than mu, so we just return this component


end intrinsic;

//=========================================================================
intrinsic Commutator(H::AlgCheRes, i::RngIntElt, p::RngMPolElt, k::RngIntElt) -> RngMPolElt
{The coefficient [y_i, p]_(s_k) of the reflection s_k in the commutator [y_i, p] for arbitrary p in R[V^*+V].}

	comm := Zero(H`yxAlgebra);
	n := H`GroupDimension;
	s := H`Group`ReflectionLibraryFlat[k];

	for mon in Monomials(p) do
		coeff := MonomialCoefficient(p, mon);
		exp := Exponents(mon);
    	mu := exp[H`xSeq];

    	//does not seem to be quicker than above
    	lambda := exp[H`ySeq];
    	comm +:= coeff*H`yxAlgebra`yEmbedding(Action(Monomial(H`yAlgebra, lambda), s`Element : Dual:=false))*H`yxAlgebra`xEmbedding(CommutatorOfYVariableWithMonomialInX(H,i,mu,k));

    end for;

    return comm;

end intrinsic;

//=========================================================================
function MultiplyMonomialInYXWithYVariable(H, p, i)
/*
	Multiply an element p of R[V + V^*] with y_i as in Lemma 1.5 in [Thi15].
*/

	W := H`Group;
	n := W`Dimension;

	exp := Exponents(p);
	pyexp := exp[H`ySeq];
	py := Monomial(H`yAlgebra, pyexp);
    pxexp := exp[H`xSeq];
    px := Monomial(H`xAlgebra, pxexp);

    if H`UseProductTable and IsDefined(H`ProductTable[i], pxexp) then
        xprod := H`ProductTable[i][pxexp];
    else

        xprod := New(AlgCheResElt);
        xprod`Parent := H;
        xprod`Element := H`GroupAlgebra!(H`yxAlgebra`Generators[i]*H`yxAlgebra`xEmbedding(px));

        n := H`GroupDimension;

        for k:=1 to #H`Group`ReflectionLibraryFlat do
        	s := H`Group`ReflectionLibraryFlat[k];
        	xprod`Element +:= H`GroupAlgebra!(H`yxAlgebra`xEmbedding(CommutatorOfYVariableWithMonomialInX(H,i,pxexp,k)))*s`Element;
        end for;

        if H`UseProductTable then
            H`ProductTable[i][pxexp] := xprod;
        end if;
    end if;

    //now, multiply with y and push group elements to the left
    if py eq 1 then
        return xprod;
    else
        prod := Zero(H);
        for g in Support(xprod`Element) do
            pyg := H`yxAlgebra`yEmbedding(Action(py, g : Dual:=false));
            prod`Element +:= pyg*Coefficient(xprod`Element,g)*H`GroupAlgebra!g;
        end for;
    end if;

    return prod;

end function;

//=========================================================================
function MultiplyPolynomialInYXWithYVariable(H, p, i)

    prod := Zero(H);

    //here, we might be able to optimize using Gordon operator!

    for mon in Monomials(p) do
        coeff := MonomialCoefficient(p,mon);
        prod`Element +:= (H`GroupAlgebra!(H`yxAlgebra!coeff))*MultiplyMonomialInYXWithYVariable(H, mon, i)`Element;
    end for;

    return prod;

end function;

//=========================================================================
function MultiplyWithYVariable(H, h, i)

    prod := Zero(H);

    for g in Support(h`Element) do
        c := Coefficient(h`Element, g);
        prod`Element +:= g*MultiplyPolynomialInYXWithYVariable(H,c,i)`Element;
    end for;

    return prod;

end function;


//=========================================================================
function MultiplyWithMonomialInY(H,h,mu)

    prod := Zero(H);	//need new element, see below
    prod`Element := h`Element;  //avoids overwriting h as done by prod=h

    for i:=1 to H`GroupDimension do
        for j:=1 to mu[i] do
            prod`Element := MultiplyWithYVariable(H, prod, i)`Element;
        end for;
    end for;

    return prod;

end function;


//=========================================================================
function MultiplyWithPolynomialInYX(H, h, p)

    prod := Zero(H);
	n := H`GroupDimension;

    for mon in Monomials(p) do
        coeff := MonomialCoefficient(p,mon);
        monexp := Exponents(mon);
        monyexp := monexp[H`ySeq];
        monxexp := monexp[H`xSeq];
        monx := H`yxAlgebra`xEmbedding(Monomial(H`xAlgebra,monxexp));

        //there is a bug in Magma with scalar multiplication in group algebras over quotient rings (up to 2.21-10 at least!)
        //must take H`GroupAlgebra!monx here
        //prod`Element +:= coeff*MultiplyWithMonomialInY(H,h,monyexp)`Element*monx;
        prod`Element +:= coeff*MultiplyWithMonomialInY(H,h,monyexp)`Element*H`GroupAlgebra!monx;
    end for;

    return prod;

end function;

//=========================================================================
function Multiply(H, h1,h2)

    prod := Zero(H);

    for g in Support(h2`Element) do
        c := Coefficient(h2`Element,g);
        prod`Element +:= MultiplyWithPolynomialInYX(H,MultiplyWithG(H,h1,g), c)`Element;
    end for;

    return prod;

end function;

//=========================================================================
intrinsic '*'(h1::AlgCheResElt, h2::AlgCheResElt) -> AlgCheResElt
{}

    H := h1`Parent;
    prod := Zero(H);

    return Multiply(H, h1,h2);

end intrinsic;
