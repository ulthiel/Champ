/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

/*
	Representation theory of restricted rational Cherednik algebras
*/


declare attributes AlgCheRes:
	CommutatorVectorsTable, //for computing standard modules
	StandardModules, //the standard modules
	SimpleModules, //their simple heads
	StandardsInSimples,
	StandardsInSimplesQuantum,
	StandardsInGroupSimples,
	StandardsInGroupSimplesQuantum,
	SimplesInGroupSimples,
	SimplesInGroupSimplesQuantum,
	ProjectivesInStandards,
	ProjectivesInStandardsQuantum,
	ProjectivesInSimples,
	ProjectivesInSimplesQuantum,
	ProjectivesInGroupSimples,
	ProjectivesInGroupSimplesQuantum,
	SimplesGradedDimension,
	StandardsGradedDimension,
	ProjectivesGradedDimension,
	StandardsAtBottomOfProjectives,
	qCharacterField;

//============================================================================
intrinsic StandardModules(H::AlgCheRes, rho::HomGrp : Rep:="Sparse",  Verbose:=true) -> ModGrOld
{
	The Standard module for rho for the restricted rational Cherednik algebra H.
}

	//We compute here Delta(rho) = K[V]_G \otimes rho as in [BR] or in [Thi15].

  if Rep ne "Sparse" and Rep ne "Dense" then
  	error "Rep has to be one of \"Sparse\" (default) or \"Dense\".";
  end if;

	Initialize(~H);

	W := H`Group;

	//we need a monomial basis of the coinvariant algebra which is sorted by degree, otherwise things will mess up. we check this here even though Basis will do this automatically.
	if exists{f : f in W`CoinvariantAlgebra`Basis | IsMonomial(f) eq false} then
		error "Basis of the coinvariant algebra has to be a monomial basis.";
	end if;

	//check if bases are compatible
	if not H`xAlgebra`Basis eq {@ H`CoinvariantAlgebraToxAlgebra(b) : b in W`CoinvariantAlgebra`Basis @} then
		error "Bases of coinvariant algebra and xAlgebra are not equal.";
	end if;

	//G-action on coinvariant algebra
	CoinvariantAlgebraGradedGModule(~W : Rep:=Rep, Verbose:=Verbose);

	//coinvariant algebra as a module over itself
	CoinvariantAlgebraGradedModule(~W : Rep:=Rep, Verbose:=Verbose);

	//some variables
	A := W`CoinvariantAlgebra;
	AW := W`CoinvariantAlgebraGradedGModule;
	AA := W`CoinvariantAlgebraGradedModule;
	Adim := Dimension(A);
	rhodim := Ncols(rho(W.1));
	Wdim := Dimension(W);
	Worder := Order(W);
	R := Codomain(H`cParameter);
	V := VectorSpace(W);

	rhospace := VectorSpace(R, rhodim); //space on which rho acts (extended to R)

	//generator gradings (w,y,x)
	alggrading := [ 0 : i in [1..Ngens(W)] ] cat [ -1 : i in [1..Wdim] ] cat [ 1 : i in [1..Wdim] ];

	//initialize matrices
	if Rep eq "Dense" then
		opmats := [* <i,d,ZeroMatrix(R,Dimension(AA`ModuleComponents[d])*rhodim,Dimension(AA`ModuleComponents[d+alggrading[i]])*rhodim)> : d in AA`Support, i in [1..#alggrading] | d+alggrading[i] in AA`Support *];
	else
		opmats := [* <i,d,SparseMatrix(R,Dimension(AA`ModuleComponents[d])*rhodim,Dimension(AA`ModuleComponents[d+alggrading[i]])*rhodim)> : d in AA`Support, i in [1..#alggrading] | d+alggrading[i] in AA`Support *];
	end if;

	//component dimensions of Standard module
	componentdimensions := [ <d,Dimension(AA`ModuleComponents[d])*rhodim> : d in AA`Support ];

	//operation matrices of rho in same type as Rep
	if Rep eq "Dense" then
		rhomats := [rho(W.i) : i in [1..Ngens(W)]];
	else
		rhomats := [SparseMatrix(rho(W.i)) : i in [1..Ngens(W)]];
	end if;

	//index blocks
	wblock := [1..Ngens(W)];
	yblock := [Ngens(W)+1..Ngens(W)+Wdim];
	xblock := [Ngens(W)+Wdim+1..Ngens(W)+2*Wdim];

	//compute action matrices of Standard module
	if Verbose then
		print "Computing Standard module.";
	end if;

	for m:=1 to #opmats do

		if Verbose then
			PrintPercentage(m, #opmats);
		end if;

		i:=opmats[m][1];	//algebra generator
		d:=opmats[m][2];	//source degree

		//g-action
		if i in wblock then
			opmats[m][3] := TensorProductMixedType( ChangeRing(AW`Matrices[i][d],R),ChangeRing(rhomats[i],R) : Rep:=Rep); //g-action on coinvariant algebra tensor g-action; it's that simple

		//x-action
		elif i in xblock then
			i:=i-Wdim-Ngens(W);	//x_i
			opmats[m][3] := TensorProductMixedType(ChangeRing(AA`Matrices[i][d],R),IdentityMatrix(R,rhodim) : Rep:=Rep); //x-action on coinvariant algebra tensor identity; it's that simple

		//y-action
		elif i in yblock then
			i:=i-Ngens(W);	//y_i
			for j:=1 to Dimension(AA`ModuleComponents[d]) do
				//action of y_i on j-th basis vector of degree d component of coinvariant algebra
				xmu := A`Basis[A`ComponentBasis[d][j]]; ////j-th basis vector of degree d component of coinvariant algebra is x^mu
				mu := Exponents(xmu);
				for r:=1 to #W`ReflectionLibraryFlat do
					v := CommutatorGradedVector(H,i,mu,r); //graded vector of [y_i,x^mu]_{s_r}; this is homogeneous of degree d-1
					if IsZero(v) then //zero vector
						continue;
					end if;
					rhos := rho(W`ReflectionLibraryFlat[r]`Element);
					for k:=1 to rhodim do
						u := TensorProduct(v, ChangeRing(rhos[k],R));
						for l in Support(u) do
							opmats[m][3][(j-1)*rhodim+k][l] +:= u[l];
						end for;
					end for;
				end for;
			end for;
		end if;

	end for;

	return GradedModuleOld(GradedModule(R, alggrading, componentdimensions, opmats : Rep:=Rep));

end intrinsic;

//============================================================================
intrinsic IsModule(H::AlgCheRes, V::ModGrOld) -> BoolElt
/*
    Intrinsic: IsModuleForRRCA

    History:
        Thursday, September 19, 2013 13:39:59: Initial.
*/
{Checks if V satisfies all relations to make it to a module for the restricted rational Cherednik algebra for G and parameter c.}

		G := H`Group;
		c := H`cParameter;
		V := RModule(V);
		matrices := ActionGenerators(V); //g, y, x
		rep := "Dense";

    FPGroup(~G);
    Gdim := Dimension(G);
    Gspace := VectorSpace(G);
    NGrels := #Relations(G`FPGroup);
    CoinvariantAlgebra(~G);
    CoinvariantAlgebra(~G`DualGroup);
    NCoinvrels := #Basis(G`HilbertIdeal);
    NDualCoinvrels := #Basis(G`DualGroup`HilbertIdeal);
    nops := 2*Gdim^2 + NGrels + 2*Gdim*Ngens(G) + Gdim*Gdim + NCoinvrels + NDualCoinvrels;

		//index blocks
		gblock := [1..Ngens(G)];
		yblock := [Ngens(G)+1..Ngens(G)+Gdim];
		xblock := [Ngens(G)+Gdim+1..Ngens(G)+2*Gdim];

    //y-y operation
    //check if y_i*y_j action is equal to y_j*y_i action
    for i in yblock do
        for j in yblock do
            if matrices[i]*matrices[j] ne matrices[j]*matrices[i] then
                print "(IsModuleForRRCA). x"*Sprint(i)*"-x"*Sprint(j)*" operation is wrong.";
                return false;
            end if;
            PrintPercentage( (i-1)*j, nops);
        end for;
    end for;

    //x-x operation
    //check if x_i*x_j action is equal to x_j*x_i action
    for i in xblock do
        for j in xblock do
            if matrices[i]*matrices[j] ne matrices[j]*matrices[i] then
                print "(IsModuleForRRCA). y"*Sprint(i)*"-y"*Sprint(j)*" operation is wrong.";
                return false;
            end if;
            PrintPercentage(Gdim^2 + (i-1)*j, nops);
        end for;
    end for;

    //g-g operations
    //check group relations (that's a bit difficult)
    relcount := 0;
    for rel in Relations(G`FPGroup) do
        relcount +:= 1;
        lhs := Eltseq(LHS(rel));
        rhs := Eltseq(RHS(rel));
        if rep eq "Sparse" then
            Standardlhs := IdentitySparseMatrix(BaseRing(V), Dimension(V));
        else
            Standardlhs := IdentityMatrix(BaseRing(V), Dimension(V));
        end if;
        for i:=1 to #lhs do
            if lhs[i] gt 0 then
                j := lhs[i];
                Standardlhs *:= matrices[j];
            else
                j := Abs(lhs[i]);
                Standardlhs *:= matrices[j]^(Order(G.j)-1);
            end if;
        end for;

        if rep eq "Sparse" then
            Standardrhs := IdentitySparseMatrix(BaseRing(V), Dimension(V));
        else
            Standardrhs := IdentityMatrix(BaseRing(V), Dimension(V));
        end if;
        for i:=1 to #rhs do
            if rhs[i] gt 0 then
                j := rhs[i];
                Standardrhs *:= matrices[j];
            else
                j := Abs(rhs[i]);
                Standardrhs *:= matrices[j]^(Order(G.j)-1);
            end if;
        end for;

        if Standardlhs ne Standardrhs then
            print "(IsModuleForRRCA). Group relation "*Sprint(rel)*" not respected.";
            return false;
        end if;

        PrintPercentage(2*Gdim^2 + relcount, nops);
    end for;


    //now we come to the interactions

    //x-g operation
    //check if x_i*g_j action is equal to g_j*x_i^g_j action (we're acting from the right in Magma)
    for i in xblock do
        for j in gblock do
            A := matrices[i]*matrices[j]; //x_i*g_j

            //B will be action of g_j*x_i^g_j
            if rep eq "Dense" then
                B := ZeroMatrix(BaseRing(V), Dimension(V), Dimension(V));
            else
                B := SparseMatrix(BaseRing(V), Dimension(V), Dimension(V));
            end if;
            v := (Transpose(G.j)^-1)[i-Ngens(G)-Gdim]; // action of g_j on x_i
            for k in Support(v) do
                B +:= v[k]*matrices[xblock[k]];
            end for;
            B := matrices[j]*B;
            if A ne B then
                print "(IsModuleForRRCA). x"*Sprint(i)*"-g"*Sprint(j)*" operation is wrong.";
                return false;
            end if;

            PrintPercentage(2*Gdim^2 + NGrels + (i-1)*Ngens(G), nops);
        end for;
    end for;

    //y-g operation
    //check if y_i*g_j action is equal to g_j*y_i^g_j action (we're acting from the right in Magma)
    for i in yblock do
        for j in gblock do
            A := matrices[i]*matrices[j]; //y_i*g_j

            //B will be action of g_j*y_i^g_j
            if rep eq "Dense" then
                B := ZeroMatrix(BaseRing(V), Dimension(V), Dimension(V));
            else
                B := SparseMatrix(BaseRing(V), Dimension(V), Dimension(V));
            end if;
            v := (G.j)[i-Ngens(G)]; // action of g_j on y_i
            for k in Support(v) do
                B +:= v[k]*matrices[yblock[k]];
            end for;
            B := matrices[j]*B;
            if A ne B then
                print "(IsModuleForRRCA). y"*Sprint(i)*"-g"*Sprint(j)*" operation is wrong.";
                return false;
            end if;

            PrintPercentage(2*Gdim^2 + NGrels + Gdim*Ngens(G) + (i-1)*Ngens(G), nops);
        end for;
    end for;

    //x-y operation
    //check if x_i*y_j action is equal to the [x_i,y_j]+(y_j*x_i) action (remember: we are in the opposite algebra!)
    for i in xblock do
        for j in yblock do
            B := matrices[i]*matrices[j]; //x_i*y_j
            A := matrices[j]*matrices[i]; //y_j*x_i

            //now we have to add the action of the commutator to B
            for s:=1 to #G`ReflectionLibraryFlat do
                sword := ElementToWord(G`ReflectionLibraryFlat[s]`Element); //reflection s as a word in the generators

                //matrix C will be the action of the reflection s on the module
                if rep eq "Dense" then
                    C := IdentityMatrix(BaseRing(V), Dimension(V));
                else
                    C := IdentitySparseMatrix(BaseRing(V), Dimension(V));
                end if;

                for k:=1 to #sword do
                    if sword[k] gt 0 then
                        p := sword[k];
                        C *:= matrices[p];
                    else
                        p := Abs(sword[k]);
                        C *:= matrices[p]^(Order(G.p)-1);
                    end if;
                end for;

                B -:= c(G`ReflectionLibraryFlat[s]`ReflectionClass)*CherednikCoefficient(Gspace.(j-Ngens(G)),Gspace.(i-Ngens(G)-Gdim),G`ReflectionLibraryFlat[s])*C;
            end for;

            if A ne B then
                print "(IsModuleForRRCA). x"*Sprint(i)*"-y"*Sprint(j)*" operation is wrong.";
                return false;
            end if;

            PrintPercentage(2*Gdim^2 + NGrels + 2*Gdim*Ngens(G) + (i-1)*Gdim, nops);

        end for;
    end for;

    //coinvariant algebra (x)
    relcount := 0;
    for rel in Basis(G`HilbertIdeal) do
        relcount +:= 1;
        mons := Monomials(rel);
        coeffs := Coefficients(rel);
        if rep eq "Dense" then
            A := ZeroMatrix(BaseRing(V), Dimension(V), Dimension(V));
        else
            A := SparseMatrix(BaseRing(V), Dimension(V), Dimension(V));
        end if;
        for i:=1 to #mons do
            exp := Exponents(mons[i]);
            if rep eq "Dense" then
                C := IdentityMatrix(BaseRing(V), Dimension(V));
            else
                C := IdentitySparseMatrix(BaseRing(V), Dimension(V));
            end if;
            for j:=1 to #exp do
                C *:= matrices[xblock[j]]^exp[j];
            end for;

            A +:= coeffs[i]*C;
        end for;

        if not IsZero(A) then
            print "(IsModuleForRRCA). Coinvariant algebra relations not respected.";
            return false;
        end if;

        PrintPercentage(2*Gdim^2 + NGrels + 2*Gdim*Ngens(G) + Gdim^2 + relcount, nops);

    end for;

    //dual coinvariant algebra (y)
    relcount := 0;
    for rel in Basis(G`DualGroup`HilbertIdeal) do
        relcount +:= 1;
        mons := Monomials(rel);
        coeffs := Coefficients(rel);
        if rep eq "Dense" then
            A := ZeroMatrix(BaseRing(V), Dimension(V), Dimension(V));
        else
            A := SparseMatrix(BaseRing(V), Dimension(V), Dimension(V));
        end if;
        for i:=1 to #mons do
            exp := Exponents(mons[i]);
            if rep eq "Dense" then
                C := IdentityMatrix(BaseRing(V), Dimension(V));
            else
                C := IdentitySparseMatrix(BaseRing(V), Dimension(V));
            end if;
            for j:=1 to #exp do
                C *:= matrices[yblock[j]]^exp[j];
            end for;

            A +:= coeffs[i]*C;
        end for;

        if not IsZero(A) then
            print "(IsModuleForRRCA). Dual coinvariant algebra relations not respected.";
            return false;
        end if;

        PrintPercentage(2*Gdim^2 + NGrels + 2*Gdim*Ngens(G) + Gdim^2 + NCoinvrels + relcount, nops);

    end for;

    print "Everything OK.";

    return true;

end intrinsic;

//============================================================================
intrinsic StandardModules(~H::AlgCheRes, i::RngIntElt)
{}

	W := H`Group;
	Representations(~W);
	if not assigned H`StandardModules then
		H`StandardModules := AssociativeArray({1..#W`Representations[0]});
	end if;
	if not IsDefined(H`StandardModules, i) then
		H`StandardModules[i] := StandardModules(H, W`Representations[0][i]);
	end if;

end intrinsic;

intrinsic StandardModules(H::AlgCheRes, i::RngIntElt) -> ModGrOld
{}

	StandardModules(~H,i);
	return H`StandardModules[i];

end intrinsic;

intrinsic StandardModules(~H::AlgCheRes)
{}

	W := H`Group;
	Representations(~W);
	for i:=1 to #W`Representations[0] do
		StandardModules(~H, i);
	end for;

end intrinsic;

//============================================================================
intrinsic SimpleModules(~H::AlgCheRes, i::RngIntElt)
{}

	W := H`Group;
	Representations(~W);
	if not assigned H`SimpleModules then
		H`SimpleModules := AssociativeArray({1..#W`Representations[0]});
	end if;
	if not IsDefined(H`SimpleModules, i) then
		StandardModules(~H, i);
		res,L,J,P := HeadOfLocalModule(H`StandardModules[i] : Rounds:=10);
		if res then
			H`SimpleModules[i] := L;
		end if;
	end if;

end intrinsic;

intrinsic SimpleModules(~H::AlgCheRes)
{}

	CharacterTable(~H`Group);
	for i:=1 to #H`Group`CharacterTable do
		SimpleModules(~H,i);
	end for;

end intrinsic;


//============================================================================
intrinsic StandardsInGroupSimples(~H::AlgCheRes, i::RngIntElt)
//
//	This uses the formula in my proceedings paper.
{}
	W := H`Group;
	CharacterTable(~W);
	FakeDegrees(~W);
	N := #W`CharacterTable;
	if not assigned H`StandardsInGroupSimples then
		H`StandardsInGroupSimples := AssociativeArray({1..N});
	end if;
	if not IsDefined(H`StandardsInGroupSimples, i) then
		V := KSpace(H`qField, N);
		lambda := W`CharacterTable[i];
		D := Zero(V);
		for j:=1 to N do
			mu := W`CharacterTable[j];
			fmu := W`FakeDegrees[j];
			muduallambda := V!Decomposition(W`CharacterTable, ComplexConjugate(mu)*lambda); // decomposition of mu* \otimes lambda
			D +:= fmu*muduallambda;
		end for;

		H`StandardsInGroupSimples[i] := D;
	end if;

end intrinsic;

intrinsic StandardsInGroupSimples(H::AlgCheRes, i::RngIntElt) -> ModTupFldElt
{}

	StandardsInGroupSimples(~H,i);
	return H`StandardsInGroupSimples[i];

end intrinsic;

intrinsic StandardsInGroupSimples(~H::AlgCheRes)
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		StandardsInGroupSimples(~H,i);
	end for;

end intrinsic;

intrinsic StandardsInGroupSimples(H::AlgCheRes) -> AlgMatElt
{}

	StandardsInGroupSimples(~H);
	N := #H`Group`CharacterTable;
	CDelta := Matrix(H`qField, N, N, [StandardsInGroupSimples(H,i) : i in [1..N]]);
	return CDelta;

end intrinsic;


//============================================================================
intrinsic InGroupSimples(H::AlgCheRes, M::ModGrOld) -> ModRngElt
{}

	W := H`Group;
	CharacterTable(~W);
	V := KSpace(H`qField, #W`CharacterTable);
	return V!DecompositionInGradedGrothendieckGroup(M, H`Group, [1..Ngens(H`Group)]);

end intrinsic;

intrinsic SimplesInGroupSimples(~H::AlgCheRes, i::RngIntElt)
{}

	W := H`Group;
	if not assigned H`SimplesInGroupSimples then
		CharacterTable(~W);
		H`SimplesInGroupSimples := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`SimplesInGroupSimples, i) then
		SimpleModules(~H,i);
		H`SimplesInGroupSimples[i] := InGroupSimples(H, H`SimpleModules[i]);
	end if;

end intrinsic;

intrinsic SimplesInGroupSimples(H::AlgCheRes, i::RngIntElt) -> ModTupFldElt
{}

	SimplesInGroupSimples(~H,i);
	return H`SimplesInGroupSimples[i];

end intrinsic;

intrinsic SimplesInGroupSimples(~H::AlgCheRes)
{}
	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		SimplesInGroupSimples(~H, i);
	end for;
end intrinsic;

intrinsic SimplesInGroupSimples(H::AlgCheRes) -> AlgMat
{}

	SimplesInGroupSimples(~H);
	W := H`Group;
	N := #H`Group`CharacterTable;
	V := KSpace(H`qField, N);
	CL := Matrix(H`qField, N, N, [V!SimplesInGroupSimples(H,i) : i in [1..N]]);
	return CL;

end intrinsic;

//============================================================================
intrinsic SimplesGradedDimension(~H::AlgCheRes, i::RngIntElt)
{}

	if not assigned H`SimplesGradedDimension then
		W := H`Group;
		CharacterTable(~W);
		H`SimplesGradedDimension := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`SimplesGradedDimension, i) then
		SimpleModules(~H,i);
		H`SimplesGradedDimension[i] := H`qField!(PoincareSeries(H`SimpleModules[i]));
	end if;

end intrinsic;

intrinsic SimplesGradedDimension(H::AlgCheRes, i::RngIntElt) -> FldElt
{}

	SimplesGradedDimension(~H,i);
	return H`SimplesGradedDimension[i];

end intrinsic;

intrinsic SimplesGradedDimension(~H::AlgCheRes)
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		SimplesGradedDimension(~H,i);
	end for;

end intrinsic;

intrinsic SimplesGradedDimension(H::AlgCheRes) -> SeqEnum
{}

	SimplesGradedDimension(~H);
	W := H`Group;
	CharacterTable(~W);
	return [ H`SimplesGradedDimension[i] : i in [1..#W`CharacterTable] ];

end intrinsic;

//============================================================================
intrinsic ProjectivesGradedDimension(~H::AlgCheRes, i::RngIntElt)
{}

	if not assigned H`ProjectivesGradedDimension then
		W := H`Group;
		CharacterTable(~W);
		H`ProjectivesGradedDimension := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`ProjectivesGradedDimension, i) then
		ProjectivesInGroupSimples(~H);
		f := Zero(H`qField);
		dec := H`ProjectivesInGroupSimples[i];
		for i:=1 to #Eltseq(dec) do
			f +:= dec[i]*Degree(H`Group`CharacterTable[i]);
		end for;
		H`ProjectivesGradedDimension[i] := f;
	end if;

end intrinsic;

intrinsic ProjectivesGradedDimension(H::AlgCheRes, i::RngIntElt) -> FldElt
{}

	ProjectivesGradedDimension(~H,i);
	return H`ProjectivesGradedDimension[i];

end intrinsic;

intrinsic ProjectivesGradedDimension(~H::AlgCheRes)
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		ProjectivesGradedDimension(~H,i);
	end for;

end intrinsic;

intrinsic ProjectivesGradedDimension(H::AlgCheRes) -> SeqEnum
{}

	ProjectivesGradedDimension(~H);
	W := H`Group;
	CharacterTable(~W);
	return [ H`ProjectivesGradedDimension[i] : i in [1..#W`CharacterTable] ];

end intrinsic;



//============================================================================
intrinsic IdentifyModule(H::AlgCheRes, M::ModGrOld : UseCharacters:=true) -> RngIntElt
{}

	d := Minimum(SequenceToSet(M`RowDegrees));
	Md := GModule(H`Group, [ Matrix(M`HomogeneousComponentMatrices[<i,d>]) : i in [1..Ngens(H`Group)]]);
	ddec := DecompositionInGrothendieckGroup(Md : UseCharacters:=UseCharacters);
	supp := Support(ddec);
	if #supp gt 1 then
		error "No unique lowest component";
	end if;
	supp := SetToSequence(supp)[1];
	if ddec[supp] ne 1 then
		error "No unique lowest component";
	end if;
	return supp;

end intrinsic;


//============================================================================
intrinsic StandardsInSimples(~H::AlgCheRes, i::RngIntElt)
{}

	if not assigned H`StandardsInSimples then
		W := H`Group;
		CharacterTable(~W);
		H`StandardsInSimples := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`StandardsInSimples, i) then
		SimplesInGroupSimples(~H);
		N := #H`Group`CharacterTable;
		V := KSpace(H`qField, N);
		CLspace := KSpaceWithBasis([ V!(H`SimplesInGroupSimples[i]) : i in [1..N] ]);
		H`StandardsInSimples[i] := V!Coordinates(CLspace, StandardsInGroupSimples(H,i));
	end if;

end intrinsic;

intrinsic StandardsInSimples(H::AlgCheRes, i::RngIntElt) -> ModTupFldElt
{}

	StandardsInSimples(~H,i);
	return H`StandardsInSimples[i];

end intrinsic;

intrinsic StandardsInSimples(~H::AlgCheRes)
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		StandardsInSimples(~H,i);
	end for;

end intrinsic;

intrinsic StandardsInSimples(H::AlgCheRes) -> AlgMat
{}

	StandardsInSimples(~H);
	W := H`Group;
	CharacterTable(~W);
	N := #W`CharacterTable;
	V := KSpace(H`qField, N);
	D := Matrix(H`qField, N, N, [V!StandardsInSimples(H,i) : i in [1..N]]);
	return D;

end intrinsic;


//=============================================================================
intrinsic ProjectivesInSimples(~H::AlgCheRes)
{}

	if not assigned H`ProjectivesInSimples then
		D := StandardsInSimples(H);
		H`ProjectivesInSimples := Transpose(D)*D;
	end if;

end intrinsic;

intrinsic ProjectivesInSimples(H::AlgCheRes) -> AlgMat
{}

	ProjectivesInSimples(~H);
	return H`ProjectivesInSimples;

end intrinsic;

intrinsic ProjectivesInSimples(H::AlgCheRes,i::RngIntElt) -> ModTupFldElt
{}

	ProjectivesInSimples(~H);
	return H`ProjectivesInSimples[i];

end intrinsic;

//=============================================================================
intrinsic ProjectivesInStandards(~H::AlgCheRes)
{}

	if not assigned H`ProjectivesInStandards then
		D := StandardsInSimples(H);
		H`ProjectivesInStandards := Transpose(D);
	end if;

end intrinsic;

intrinsic ProjectivesInStandards(H::AlgCheRes) -> AlgMat
{}

	ProjectivesInStandards(~H);
	return H`ProjectivesInStandards;

end intrinsic;

intrinsic ProjectivesInStandards(H::AlgCheRes,i::RngIntElt) -> ModTupFldElt
{}

	ProjectivesInStandards(~H);
	return H`ProjectivesInStandards[i];

end intrinsic;

//=============================================================================
procedure qCharacterField(~H)

	if not assigned H`qCharacterField then
		CharacterNames(~H`Group);
		R:=PolynomialRing(Integers(), #H`Group`CharacterTable);
		AssignNames(~R,H`Group`CharacterNames);
		K:=RationalFunctionField(R,1);
		AssignNames(~K,["q"]);
		H`qCharacterField := K;
	end if;

end procedure;

//=============================================================================
intrinsic ProjectivesInStandardsQuantum(~H::AlgCheRes,i::RngIntElt)
{}

	if not assigned H`ProjectivesInStandardsQuantum then
		W := H`Group;
		CharacterTable(~W);
		H`ProjectivesInStandardsQuantum := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`ProjectivesInStandardsQuantum, i) then
		ProjectivesInStandards(~H);
		qCharacterField(~H);
		phi := hom<H`qField -> H`qCharacterField | [H`qCharacterField.1]>;
		f:=Zero(H`qCharacterField);
		dec :=  ProjectivesInStandards(H,i);
		R := BaseRing(H`qCharacterField);
		for j in Support(dec) do
			f +:= phi(dec[j])*R.j;
		end for;
		H`ProjectivesInStandardsQuantum[i] := f;
	end if;

end intrinsic;

intrinsic ProjectivesInStandardsQuantum(H::AlgCheRes,i::RngIntElt) -> FldFunratElt
{}

		ProjectivesInStandardsQuantum(~H,i);
		return H`ProjectivesInStandardsQuantum[i];

end intrinsic;

intrinsic ProjectivesInStandardsQuantum(H::AlgCheRes) -> SeqEnum
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		ProjectivesInStandardsQuantum(~H,i);
	end for;
	return [H`ProjectivesInStandardsQuantum[i] : i in [1..#W`CharacterTable]];

end intrinsic;

//=============================================================================
intrinsic ProjectivesInGroupSimples(~H::AlgCheRes)
{}

	if not assigned H`ProjectivesInGroupSimples then
		W := H`Group;
		CharacterTable(~W);
		H`ProjectivesInGroupSimples := AssociativeArray({1..#W`CharacterTable});
	end if;
	C := ProjectivesInSimples(H);
	CL := SimplesInGroupSimples(H);
	N := #H`Group`CharacterTable;
	V := KSpace(H`qField, N);
	H`ProjectivesInGroupSimples := C*CL;

end intrinsic;


intrinsic ProjectivesInGroupSimples(H::AlgCheRes,i::RngIntElt) -> AlgMat
{}

	ProjectivesInGroupSimples(~H);
	return H`ProjectivesInGroupSimples[i];

end intrinsic;

intrinsic ProjectivesInGroupSimples(H::AlgCheRes) -> AlgMat
{}

	ProjectivesInGroupSimples(~H);
	return H`ProjectivesInGroupSimples;

end intrinsic;

//=============================================================================
intrinsic StandardsInSimplesQuantum(~H::AlgCheRes,i::RngIntElt)
{}

	if not assigned H`StandardsInSimplesQuantum then
		W := H`Group;
		CharacterTable(~W);
		H`StandardsInSimplesQuantum := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`StandardsInSimplesQuantum, i) then
		StandardsInSimples(~H);
		qCharacterField(~H);
		phi := hom<H`qField -> H`qCharacterField | [H`qCharacterField.1]>;
		f:=Zero(H`qCharacterField);
		dec :=  StandardsInSimples(H,i);
		R := BaseRing(H`qCharacterField);
		for j in Support(dec) do
			f +:= phi(dec[j])*R.j;
		end for;
		H`StandardsInSimplesQuantum[i] := f;
	end if;

end intrinsic;


intrinsic StandardsInSimplesQuantum(H::AlgCheRes,i::RngIntElt) -> FldFunratElt
{}

		StandardsInSimplesQuantum(~H,i);
		return H`StandardsInSimplesQuantum[i];

end intrinsic;

intrinsic StandardsInSimplesQuantum(H::AlgCheRes) -> SeqEnum
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		StandardsInSimplesQuantum(~H,i);
	end for;
	return [H`StandardsInSimplesQuantum[i] : i in [1..#W`CharacterTable]];

end intrinsic;


//=============================================================================
intrinsic ProjectivesInGroupSimplesQuantum(~H::AlgCheRes,i::RngIntElt)
{}

	if not assigned H`ProjectivesInGroupSimplesQuantum then
		W := H`Group;
		CharacterTable(~W);
		H`ProjectivesInGroupSimplesQuantum := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`ProjectivesInGroupSimplesQuantum, i) then
		ProjectivesInGroupSimples(~H);
		qCharacterField(~H);
		phi := hom<H`qField -> H`qCharacterField | [H`qCharacterField.1]>;
		f:=Zero(H`qCharacterField);
		dec :=  ProjectivesInGroupSimples(H,i);
		R := BaseRing(H`qCharacterField);
		for j in Support(dec) do
			f +:= phi(dec[j])*R.j;
		end for;
		H`ProjectivesInGroupSimplesQuantum[i] := f;
	end if;

end intrinsic;

intrinsic ProjectivesInGroupSimplesQuantum(H::AlgCheRes,i::RngIntElt)
{}

	ProjectivesInGroupSimplesQuantum(~H,i);

end intrinsic;

intrinsic ProjectivesInGroupSimplesQuantum(~H::AlgCheRes)
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		ProjectivesInGroupSimplesQuantum(~H,i);
	end for;
	H`ProjectivesInGroupSimplesQuantum := [H`ProjectivesInGroupSimplesQuantum[i] : i in [1..#W`CharacterTable]];

end intrinsic;

intrinsic ProjectivesInGroupSimplesQuantum(H::AlgCheRes) -> SeqEnum
{}

	ProjectivesInGroupSimplesQuantum(~H);
	return H`ProjectivesInGroupSimplesQuantum;

end intrinsic;

//=============================================================================
intrinsic SimplesInGroupSimplesQuantum(~H::AlgCheRes,i::RngIntElt)
{}

	if not assigned H`SimplesInGroupSimplesQuantum then
		W := H`Group;
		CharacterTable(~W);
		H`SimplesInGroupSimplesQuantum := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`SimplesInGroupSimplesQuantum, i) then
		SimplesInGroupSimples(~H);
		qCharacterField(~H);
		phi := hom<H`qField -> H`qCharacterField | [H`qCharacterField.1]>;
		f:=Zero(H`qCharacterField);
		dec :=  SimplesInGroupSimples(H,i);
		R := BaseRing(H`qCharacterField);
		for j in Support(dec) do
			f +:= phi(dec[j])*R.j;
		end for;
		H`SimplesInGroupSimplesQuantum[i] := f;
	end if;

end intrinsic;

intrinsic SimplesInGroupSimplesQuantum(H::AlgCheRes,i::RngIntElt)
{}

	SimplesInGroupSimplesQuantum(~H,i);

end intrinsic;

intrinsic SimplesInGroupSimplesQuantum(~H::AlgCheRes)
{}

	W := H`Group;
	CharacterTable(~W);
	for i:=1 to #W`CharacterTable do
		SimplesInGroupSimplesQuantum(~H,i);
	end for;
	H`SimplesInGroupSimplesQuantum := [H`SimplesInGroupSimplesQuantum[i] : i in [1..#W`CharacterTable]];

end intrinsic;

intrinsic SimplesInGroupSimplesQuantum(H::AlgCheRes) -> SeqEnum
{}

	SimplesInGroupSimplesQuantum(~H);
	return H`SimplesInGroupSimplesQuantum;

end intrinsic;


//=============================================================================
intrinsic StandardsAtBottomOfProjectives(~H::AlgCheRes, i::RngIntElt)
{}

	if not assigned H`StandardsAtBottomOfProjectives then
		W := H`Group;
		CharacterTable(~W);
		H`StandardsAtBottomOfProjectives := AssociativeArray({1..#W`CharacterTable});
	end if;
	if not IsDefined(H`StandardsAtBottomOfProjectives, i) then
		f := ProjectivesInStandardsQuantum(H,i);
		if Denominator(f) ne 1 then
			error "Some problem with LaurentSeriesRing/RationalFunctionField crap";
		end if;
		f:=Numerator(f);
		lc := LeadingCoefficient(f);
		deg := Degree(f);
		//make sure there's really just one at the bottom! Should be true if didn't
		//mess up the proof!
		if #Terms(lc) ne 1 or Coefficients(lc) ne [1] then
			error "Oh no: more than one standard at bottom! You proved nonsense!";
		end if;
		num := Eltseq(lc)[1];
		H`StandardsAtBottomOfProjectives[i] := <num,deg>;
	end if;

end intrinsic;

intrinsic StandardsAtBottomOfProjectives(H::AlgCheRes, i::RngIntElt) -> Tup
{}

	StandardsAtBottomOfProjectives(~H,i);
	return H`StandardsAtBottomOfProjectives[i];

end intrinsic;

intrinsic StandardsAtBottomOfProjectives(~H::AlgCheRes)
{}

	if not assigned H`StandardsAtBottomOfProjectives then
		W := H`Group;
		CharacterTable(~W);
		for i:=1 to #W`CharacterTable do
			StandardsAtBottomOfProjectives(~H,i);
		end for;
	end if;

end intrinsic;

intrinsic StandardsAtBottomOfProjectives(H::AlgCheRes) -> SeqEnum
{}

	StandardsAtBottomOfProjectives(~H);
	W := H`Group;
	CharacterTable(~W);
	return [H`StandardsAtBottomOfProjectives[i] : i in [1..#W`CharacterTable]];

end intrinsic;

intrinsic MediaWiki(H::AlgCheRes)
{}

	StandardsInSimples(~H);
	fams,sigma := Families(StandardsInSimples(H));
	charnames := [H`Group`CharacterNames[i] : i in [1..#H`Group`CharacterTable]];

	cparam := "[";
	counter := 0;
	for i in Domain(H`cParameter) do
		cparam *:= Sprint(H`cParameter(i));
		counter +:= 1;
		if counter lt #Domain(H`cParameter) then
			cparam *:= ",";
		end if;
	end for;
	cparam *:= "]";
	print "";

	printf "= Parameter $c=%o$ =\n\n", cparam;

	printf "== Representation theory ==\n\n";

	printf "=== Block structure ===\n";

	printf "==== Families ====\n";
	famtable := [];
	for i:=1 to #fams do
		famstr := "";
		famstr *:="{";
		for j:=1 to #fams[i] do
			famstr*:=H`Group`CharacterNames[fams[i][j]];
			if j lt #fams[i] then
				famstr*:=",";
			end if;
		end for;
		famstr*:="}";
		famtable cat:=[ [Sprint(#fams[i]), famstr] ];
	end for;
	MediaWiki(famtable, "Default" : ColHeader:=["Size", "Characters"], RowHeader:=[Sprint(i) : i in [1..#fams]], Legend:="#", ScrollingTable:=true);
	printf "\n";

	printf "==== Tilting permutation ====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	bottom := StandardsAtBottomOfProjectives(H);
	for F in fams do
		bottomF := [ [Sprint(H`Group`CharacterNames[x[1]])*"["*Sprint(x[2])*"]"] : x in [bottom[i] : i in F] ];
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(bottomF, "Default" : RowHeader:=charnamesF, ColHeader:=["&lambda;<sup>h</sup>"], Legend:="&lambda;",ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";


	printf "=== Projectives ===\n";

	printf "==== in standards ====\n";

	printf "===== by character =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	D := ProjectivesInStandards(H);
	for F in fams do
		DF := Submatrix(D, IndexedSetToSequence(F), IndexedSetToSequence(F));
		charnamesF := [H`Group`CharacterNames[i] : i in F];
		MediaWiki(DF, "Latex" : ColHeader:=charnamesF, RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "===== by degree =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	D := ProjectivesInStandardsQuantum(H);
	R:=PolynomialRing(BaseRing(H`qCharacterField));
	for F in fams do
		for i in F do
			if Denominator(D[i]) ne 1 then
				error "Something wrong with LaurentSeries/RationalFunction crap";
			end if;
		end for;
		DF := [R!(Numerator(D[i])) : i in F];
		degrees := {};
		for f in DF do
			degrees join:=SequenceToSet(Support(f));
		end for;
		degrees := SetToSequence(degrees);
		A := ZeroMatrix(BaseRing(R),#F,#degrees);
		for i:=1 to #F do
			for j:=1 to #degrees do
				A[i,j] := Coefficient(DF[i],degrees[j]);
			end for;
		end for;
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(A, "Default" : ColHeader:=[Sprint(degrees[i]) : i in [1..#degrees]], RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "==== in group simples ====\n";

	printf "===== by character =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	for F in fams do
		DF := [ Eltseq(ProjectivesInGroupSimples(H,i)) : i in F ];
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(DF, "Latex" : ColHeader:=charnames, RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";printf "</div>\n\n";

	printf "===== by degree =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	D := ProjectivesInGroupSimplesQuantum(H);
	R:=PolynomialRing(BaseRing(H`qCharacterField));
	for F in fams do
		for i in F do
			if Denominator(D[i]) ne 1 then
				error "Something wrong with LaurentSeries/RationalFunction crap";
			end if;
		end for;
		DF := [R!(Numerator(D[i])) : i in F];
		degrees := {};
		for f in DF do
			degrees join:=SequenceToSet(Support(f));
		end for;
		degrees := SetToSequence(degrees);
		A := ZeroMatrix(BaseRing(R),#F,#degrees);
		for i:=1 to #F do
			for j:=1 to #degrees do
				A[i,j] := Coefficient(DF[i],degrees[j]);
			end for;
		end for;
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(A, "Default" : ColHeader:=[Sprint(degrees[i]) : i in [1..#degrees]], RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "==== Graded dimension ====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	dims:=ProjectivesGradedDimension(H);
	for F in fams do
		dimsF := [ Numerator(dims[i]) : i in F ];
		maxdeg := 0;
		for f in dimsF do
			if Degree(f) gt maxdeg then
				maxdeg := Degree(f);
			end if;
		end for;
		A := ZeroMatrix(Integers(),#F, maxdeg+1);
		for i:=1 to #F do
			for j:=0 to maxdeg do
				A[i][j+1] := Coefficient(dimsF[i],j);
			end for;
		end for;
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(A, "Default" : ColHeader:=[Sprint(i) : i in [0..maxdeg]], RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "=== Standards ===\n";

	printf "==== in simples ====\n";

	printf "===== by character =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	D := StandardsInSimples(H);
	for F in fams do
		DF := Submatrix(D, IndexedSetToSequence(F), IndexedSetToSequence(F));
		charnamesF := [H`Group`CharacterNames[i] : i in F];
		MediaWiki(DF, "Latex" : ColHeader:=charnamesF, RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "===== by degree =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	D := StandardsInSimplesQuantum(H);
	R:=PolynomialRing(BaseRing(H`qCharacterField));
	for F in fams do
		for i in F do
			if Denominator(D[i]) ne 1 then
				error "Something wrong with LaurentSeries/RationalFunction crap";
			end if;
		end for;
		DF := [R!(Numerator(D[i])) : i in F];
		degrees := {};
		for f in DF do
			degrees join:=SequenceToSet(Support(f));
		end for;
		degrees := SetToSequence(degrees);
		A := ZeroMatrix(BaseRing(R),#F,#degrees);
		for i:=1 to #F do
			for j:=1 to #degrees do
				A[i,j] := Coefficient(DF[i],degrees[j]);
			end for;
		end for;
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(A, "Default" : ColHeader:=[Sprint(degrees[i]) : i in [1..#degrees]], RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "==== in group simples ====\n";
	printf "===== by character =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	for F in fams do
		DF := [ Eltseq(SimplesInGroupSimples(H,i)) : i in F ];
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(DF, "Latex" : ColHeader:=charnames, RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "===== by degree =====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	D := SimplesInGroupSimplesQuantum(H);
	R:=PolynomialRing(BaseRing(H`qCharacterField));
	for F in fams do
		for i in F do
			if Denominator(D[i]) ne 1 then
				error "Something wrong with LaurentSeries/RationalFunction crap";
			end if;
		end for;
		DF := [R!(Numerator(D[i])) : i in F];
		degrees := {};
		for f in DF do
			degrees join:=SequenceToSet(Support(f));
		end for;
		degrees := SetToSequence(degrees);
		A := ZeroMatrix(BaseRing(R),#F,#degrees);
		for i:=1 to #F do
			for j:=1 to #degrees do
				A[i,j] := Coefficient(DF[i],degrees[j]);
			end for;
		end for;
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(A, "Default" : ColHeader:=[Sprint(degrees[i]) : i in [1..#degrees]], RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

	printf "==== Graded dimension ====\n";
	//printf "<div class=\"mw-collapsible mw-collapsed\">\n";
	dims:=SimplesGradedDimension(H);
	for F in fams do
		dimsF := [ Numerator(dims[i]) : i in F ];
		maxdeg := 0;
		for f in dimsF do
			if Degree(f) gt maxdeg then
				maxdeg := Degree(f);
			end if;
		end for;
		A := ZeroMatrix(Integers(),#F, maxdeg+1);
		for i:=1 to #F do
			for j:=0 to maxdeg do
				A[i][j+1] := Coefficient(dimsF[i],j);
			end for;
		end for;
		charnamesF := [ H`Group`CharacterNames[i] : i in F ];
		MediaWiki(A, "Default" : ColHeader:=[Sprint(i) : i in [0..maxdeg]], RowHeader:=charnamesF, ScrollingTable:=true);
	end for;
	//printf "</div>\n";
	printf "\n";

end intrinsic;
