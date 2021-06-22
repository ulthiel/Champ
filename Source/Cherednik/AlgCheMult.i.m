/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
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
function CommutatorOfYVariableWithMonomialInX(H, i, mu, k)
/*
	The coefficient [y_i, x^\mu]_{s_k} of the reflection s_k in the commutator [y_i, x^\mu] as in [Thi15], Lemma 1.15. This is an element of R[V] but we coerce it already into R[V^*+V].
*/

	W := H`Group;
	s := W`ReflectionLibraryFlat[k];

	comm := H`xAlgebra`Zero;

	if IsEmpty(mu) then
		return comm;
	end if;

	//if UseCommutatorsTable, try to find commutator in the table
	if H`UseCommutatorsTable then
		if IsDefined(H`CommutatorsTable[i][k], mu) then
			return H`CommutatorsTable[i][k][mu];
		end if;
	end if;

	//otherwise, compute it

	n := H`GroupDimension;

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

		leftpart := Monomial(H`xAlgebra,leftpartexp); //x_1^{\mu_1}*...*x_{t-1}^{\mu_{t-1}}

		rightpartexp[t] := 0;
		rightpart := Monomial(H`xAlgebra,rightpartexp); //x_{t+1}^{\mu_{t+1}}*...*x_{n}^{\mu_n}

		//set leftpartexp for next round
		leftpartexp[t] := mu[t];

		rightparts := Action(rightpart, s`DualElement : Dual:=false); //will be non-zero since obtained from group action on non-zero element

		midpart := H`xAlgebra`Zero;	//the middle part
		u := mu[t]-1; //will be the exponent in the following
		for l:=0 to u do
        	midpart +:= H`xAlgebra`Generators[t]^l*Action(H`xAlgebra`Generators[t]^u, s`DualElement : Dual:=false);
        	u -:= 1;
    	end for;
    	if IsZero(midpart) then
    		continue;
    	end if;
		comm +:= checoeff*leftpart*midpart*rightparts;
	end for;

	comm *:= H`cParameter(s`ReflectionClass);

	comm := H`yxAlgebra`xEmbedding(comm); //coerce into R[V^*+V]

	if H`UseCommutatorsTable then
		H`CommutatorsTable[i][k][mu] := comm;
	end if;

	return comm;

end function;

//=========================================================================
intrinsic Commutator(H::AlgChe, i::RngIntElt, p::RngMPolElt, k::RngIntElt) -> RngMPolElt
{The coefficient [y_i, p]_(s_k) of the reflection s_k in the commutator [y_i, p].}

	comm := Zero(H`yxAlgebra);
	n := H`GroupDimension;
	s := H`Group`ReflectionLibraryFlat[k];

	for mon in Monomials(p) do
		coeff := MonomialCoefficient(p, mon);
		exp := Exponents(mon);
    	mu := exp[H`xSeq];

    	//does not seem to be quicker than above
    	lambda := exp[H`ySeq];
    	comm +:= coeff*H`yxAlgebra`yEmbedding(Action(Monomial(H`yAlgebra, lambda), s`Element : Dual:=false))*CommutatorOfYVariableWithMonomialInX(H,i,mu,k);

    end for;

    return comm;

end intrinsic

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
        xprod := New(AlgCheElt);
        xprod`Parent := H;
        xprod`Element := H`GroupAlgebra!(H`yxAlgebra`Generators[i]*H`yxAlgebra`xEmbedding(px));

        n := H`GroupDimension;

        for k:=1 to #H`Group`ReflectionLibraryFlat do
        	s := H`Group`ReflectionLibraryFlat[k];
        	xprod`Element +:= H`GroupAlgebra!(CommutatorOfYVariableWithMonomialInX(H,i,pxexp,k))*s`Element;
        end for;

        if H`tParameter ne 0 then
            tpart := Zero(H);
            for j:=1 to n do
                if j ne i then
                    continue;
                end if;
                if pxexp[j] lt 1 then
                    continue;
                end if;
                pxexpj := pxexp;
                pxexpj[j] -:= 1;
                xpart := H`yxAlgebra`xEmbedding(Monomial(H`xAlgebra, pxexpj));
                tpart`Element +:= H`tParameter*pxexp[j]*xpart;
            end for;
            xprod`Element +:= tpart`Element;
        end if;

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
        prod`Element +:= (H`yxAlgebra!coeff)*MultiplyMonomialInYXWithYVariable(H, mon, i)`Element;
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

	//print "In ";
	//print h,p;

    for mon in Monomials(p) do
        coeff := MonomialCoefficient(p,mon);
        monexp := Exponents(mon);
        monyexp := monexp[H`ySeq];
        monxexp := monexp[H`xSeq];
        monx := H`yxAlgebra`xEmbedding(Monomial(H`xAlgebra,monxexp));
        prod`Element +:= coeff*MultiplyWithMonomialInY(H,h,monyexp)`Element*monx;
    end for;

	//print "Out";
	//print prod;

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
intrinsic '*'(h1::AlgCheElt, h2::AlgCheElt) -> AlgCheElt
{}

    H := h1`Parent;
    prod := Zero(H);

    return Multiply(H, h1,h2);

end intrinsic;
