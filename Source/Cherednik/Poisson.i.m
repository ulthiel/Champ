/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


//=========================================================================
intrinsic LiftElementToPoissonAlgebra(h::AlgCheElt) -> AlgCheElt
{}

    H := Parent(h);
    hP := Zero(H`PoissonAlgebra);
    for g in Support(h`Element) do
        coeff := Coefficient(h`Element,g);
        coeffPoisson := Zero(H`PoissonAlgebra`yxAlgebra);
        for t in Terms(coeff) do
            tcoeff := Coefficients(t)[1];
            tmon := Monomials(t)[1];
            tcoeffPoisson := H`PoissonBaseRingEmbedding(tcoeff);
            coeffPoisson +:= tcoeffPoisson*Monomial(H`PoissonAlgebra`yxAlgebra, Exponents(tmon));
        end for;
        hP`Element +:= coeffPoisson*H`PoissonAlgebra`GroupAlgebra!g;
    end for;

    return hP;

end intrinsic;;

//=========================================================================
intrinsic ProjectElementFromPoissonAlgebra(H::AlgChe, hP::AlgCheElt) -> AlgCheElt
{}
	//cut out the t-Part of hP

    HP := Parent(hP);
    h := Zero(H);

	for g in Support(hP`Element) do
		coeffP := Coefficient(hP`Element,g);
		coeff := Zero(H`yxAlgebra);
		//extract t-Part from coeffP. This lives in R[t/t^2][y1,...,yn,x1,...,xn].
		//have to go through all terms of this and extract t-part.
		for term in Terms(coeffP) do
			termcoeff := Coefficients(term)[1];	//coefficient of this term (is in the t-algebra)
			termmon := Monomials(term)[1];	//monomial of this term (is in the x's and y's)
			coeff +:= H`PoissonBaseRingProjection(termcoeff)*Monomial(H`yxAlgebra, Exponents(termmon));
		end for;
		h`Element +:= coeff*H`GroupAlgebra!g;
	end for;

    return h;

end intrinsic;

//=========================================================================
intrinsic PoissonBracket(h1::AlgCheElt, h2::AlgCheElt) -> AlgCheElt
{}
    H := h1`Parent;

    if not H`Poisson then
        error "Please create the Cherednik algebra with the option \"Poisson\" enabled to compute Poisson brackets.";
    end if;

    //lift to PoissonAlgebra
    h1P := LiftElementToPoissonAlgebra(h1);
    h2P := LiftElementToPoissonAlgebra(h2);

    commP := h1P*h2P - h2P*h1P;

    return ProjectElementFromPoissonAlgebra(H, commP);

end intrinsic;
