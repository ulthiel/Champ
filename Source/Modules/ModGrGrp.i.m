/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	Graded G-modules (with group algebra sitting in degree zero). This is a derivative of ModGr.
*/

declare type ModGrGrp[ModGrGrpElt];

declare attributes ModGrGrp:
    BaseRing,
    ModuleComponents,
    Support,
    AlgebraDegrees,
    Dimension,
    ComponentDimension,
    Matrices,
    Group,
    GModules, //the GModules (of type ModGrp) for each component
    RModules, //the GModules (of type ModRng) for each component
    GradedModule, //this will be the corresponding ModGr.
    DecompositionInGradedGrothendieckGroup; //G-characters of the components

declare attributes ModGrGrpElt:
	Parent,
	Components,
	Support,
	GradedModuleElement; //this will be the corresponding ModGrElt

//============================================================================
intrinsic GradedGModule(R::Rng, G::Grp, Components::SeqEnum[Tup], Q::List : Rep:="Dense") -> ModGrGrp
{
	The graded G-module defined by the given data (with the group algebra sitting in degree zero). Internally, this is of type ModGr.
}

	M := New(ModGrGrp);
	M`GradedModule := GradedModule(R, [0 : i in [1..Ngens(G)]], Components, Q : Rep:=Rep);
	M`BaseRing := M`GradedModule`BaseRing;
	M`ModuleComponents := M`GradedModule`ModuleComponents;
	M`AlgebraDegrees := M`GradedModule`AlgebraDegrees;
	M`Dimension := M`GradedModule`Dimension;
	M`ComponentDimension := M`GradedModule`ComponentDimension;
	M`Matrices := M`GradedModule`Matrices;
	M`Support := M`GradedModule`Support;
	M`Group := G;

	return M;

end intrinsic;

//============================================================================
procedure AttachGradedModuleElement(~m)
	/*
		Here, we simply assign m`GradedModuleElement from the data in m.
	*/

	m`GradedModuleElement := New(ModGrElt);
	m`GradedModuleElement`Parent := m`Parent`GradedModule;
	m`GradedModuleElement`Support := m`Support;
	m`GradedModuleElement`Components := m`Components;

end procedure;

//============================================================================
procedure AttachGradedGModuleElement(~m)
	/*
		Here, we simply assign the data of m from the data in m`GradedModuleElement.
	*/

	m`Support := m`GradedModuleElement`Support;
	m`Components := m`GradedModuleElement`Components;

end procedure;

//============================================================================
intrinsic AttachGModules(~M::ModGrGrp : Verbose:=false)
{
	Attaches the GModules.
}

	if assigned M`GModules then
		return;
	end if;

	M`GModules := AssociativeArray(M`Support);

	count := 0;
	total := #M`Support;

	for d in M`Support do
		M`GModules[d] := GModule(M`Group, [ Matrix(M`Matrices[i][d]) : i in [1..Ngens(M`Group)] ]);
		if Verbose then
			count +:= 1;
			PrintPercentage(count, total);
		end if;
	end for;

	if Verbose then
		print "";
	end if;

end intrinsic;


//============================================================================
intrinsic AttachRModules(~M::ModGrGrp : Verbose:=false)
{
	Attaches the GModules.
}

	if assigned M`RModules then
		return;
	end if;

	M`RModules := AssociativeArray(M`Support);

	count := 0;
	total := #M`Support;

	for d in M`Support do
		M`RModules[d] := RModule([ Matrix(M`Matrices[i][d]) : i in [1..Ngens(M`Group)] ]);
		if Verbose then
			count +:= 1;
			PrintPercentage(count, total);
		end if;
	end for;

	if Verbose then
		print "";
	end if;

end intrinsic;


//============================================================================
intrinsic Print(M::ModGrGrp)
{}

   	printf "Graded G-module of dimension "*Sprint(&+[M`ComponentDimension[i] : i in [1..#M`Support]])*" and components "*Sprint( [<M`Support[i],M`ComponentDimension[i] > : i in [1..#M`Support]] )*".";

end intrinsic;

//============================================================================
intrinsic Parent(m::ModGrGrpElt) -> ModGrGrp
{}

	return m`Parent;

end intrinsic;

//============================================================================
intrinsic Print(m::ModGrGrpElt)
{}

   Print(m`GradedModuleElement);

end intrinsic;

//============================================================================
intrinsic Zero(M::ModGrGrp) -> ModGrGrpElt
{Returns the zero element of a graded module.}

	z := New(ModGrGrpElt);
	z`GradedModuleElement := Zero(M`GradedModule);
	z`Parent := M;
	AttachGradedGModuleElement(~z);
	return z;

end intrinsic;

//============================================================================
intrinsic '.' (M::ModGrGrp, j::RngIntElt) -> ModGrGrpElt
{}

	m := New(ModGrGrpElt);
	m`Parent := M;
	m`GradedModuleElement := (M`GradedModule).j;
	AttachGradedGModuleElement(~m);
	return m;

end intrinsic;

//============================================================================
intrinsic '.' (M::ModGrGrp, q::Tup) -> ModGrGrpElt
{For a pair q := <d,j> this intrinsic returns the j-th basis vector of the component of degree d.}

	m := New(ModGrGrpElt);
	m`Parent := M;
	m`GradedModuleElement := (M`GradedModule).q;
	AttachGradedGModuleElement(~m);
	return m;

end intrinsic;

//============================================================================
intrinsic '^' (m::ModGrGrpElt, i::RngIntElt) -> ModGrElt
{Operation of the i-th Algebra generator on the graded module element m.}

	x := New(ModGrGrpElt);
	x`Parent := m`Parent;
	x`GradedModuleElement := (m`GradedModuleElement)^i;
	AttachGradedGModuleElement(~x);
	return x;

end intrinsic;

//============================================================================
intrinsic Spin(M::ModGrGrp, U::List : UseModules:=true, Verbose:=false) -> List
{
	Returns a basis of the graded submodule of a graded module M generated by homogenous elements listet in U.
}

	if not UseModules then
		B := Spin(M`GradedModule, [* u`GradedModuleElement : u in U *] : Verbose:=Verbose);
		Bnew := [* *];
		for b in B do
			bnew := New(ModGrGrpElt);
			bnew`Parent := M;
			bnew`GradedModuleElement := b;
			AttachGradedGModuleElement(~bnew);
			Append(~Bnew, bnew);
		end for;
		if Verbose then
			print "";
		end if;
		return B;
	end if;

	//the construction of the GModules takes some time, so we prefer attaching the RModules if the GModules aren't assigned already
	if assigned M`RModules then
		use := "R";
	elif assigned M`GModules then
		use := "G";
	else
		AttachRModules(~M);
		use := "R";
	end if;

	for u in U do
		if not #u`Support eq 1 then
			error "Given list does contain non-homogeneous vectors.";
		end if;
	end for;

	degrees := {u`Support[1] : u in U};
	span := AssociativeArray(degrees);

	total := #degrees;
	count := 0;

	for d in degrees do
		Ud := {u`Components[d] : u in U | u`Support[1] eq d};
		if use eq "R" then
			span[d] := sub<M`RModules[d]|Ud>;
		elif use eq "G" then
			span[d] := sub<M`GModules[d]|Ud>;
		else
			error "Error."; //shouldn't come here
		end if;
		if IsZero(span[d]) then
			Remove(~span, d);
		end if;
		if Verbose then
			count +:= 1;
			PrintPercentage(count, total);
		end if;
	end for;

	B := [* *];
	for d in degrees do
		for b in Basis(span[d]) do
			x := New(ModGrGrpElt);
			x`Parent := M;
			x`Support := {d};
			x`Components := AssociativeArray(M`Support);
			x`Components[d] := b;
			AttachGradedModuleElement(~x);
			Append(~B, x);
		end for;
	end for;

	if Verbose then
		print "";
	end if;

	return B;

end intrinsic;


//============================================================================
intrinsic DecompositionInGradedGrothendieckGroup(~M::ModGrGrp : Verbose:=false, UseCharacters:=true)
{}

	if assigned M`DecompositionInGradedGrothendieckGroup then
		return;
	end if;

	AttachGModules(~M);

	P<q> := PolynomialRing(Integers());
    p := Characteristic(M`BaseRing);
    G := M`Group;
    Modules(~G,p);
    V := RModule(P, #G`Modules[p]);

    dec := Zero(V);
    for d in M`Support do
        Mddec := DecompositionInGrothendieckGroup(M`GModules[d] : UseCharacters:=UseCharacters);
        for j:=1 to #M`Group`Modules[p] do
            dec[j] +:= Mddec[j]*q^d;
        end for;
    end for;

    M`DecompositionInGradedGrothendieckGroup := dec;

end intrinsic;

//============================================================================
intrinsic DecompositionInGradedGrothendieckGroup(M::ModGrGrp : Verbose:=false, UseCharacters:=true) -> ModRngElt
{}

	DecompositionInGradedGrothendieckGroup(~M : Verbose:=Verbose, UseCharacters:=UseCharacters);

	return M`DecompositionInGradedGrothendieckGroup;

end intrinsic;
