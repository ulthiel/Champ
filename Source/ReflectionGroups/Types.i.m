//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Constructors for (irreducible) complex reflection groups.
// The functions will simply retrieve the particular model from the database
// (default) or use Magma's function.
//
//==============================================================================

//==============================================================================
// For complex reflection groups that are Coxeter groups we prefer the Coxeter
// name. Also for cyclic groups we prefer "Cyc".
//==============================================================================
intrinsic ComplexReflectionGroupName(n::RngIntElt) -> MonStgElt
{}

	if n eq 23 then
		name := "H3";
	elif n eq 28 then
		name := "F4";
	elif n eq 30 then
		name := "H4";
	elif n eq 35 then
		name := "E6";
	elif n eq 36 then
		name := "E7";
	elif n eq 37 then
		name := "E8";
	else
		name := "G"*Sprint(n);
	end if;

	return name;

end intrinsic;

intrinsic ComplexReflectionGroupName(m::RngIntElt, p::RngIntElt, n::RngIntElt) -> MonStgElt
{}

	if m eq 1 and p eq 1 then
		name := "A"*Sprint(n-1);
	elif m eq 2 and p eq 1 then
		name := "B"*Sprint(n);
	elif m eq 2 and p eq 2 then
		name := "D"*Sprint(n);
	elif p eq 1 and n eq 1 then
		name := "Cyc"*Sprint(m);
	elif m eq p and n eq 2 then
		name := "Dih"*Sprint(m);
	else
		name := Sprintf("G%o_%o_%o", m, p, n);
	end if;

	return name;

end intrinsic;

//==============================================================================
// By default the following function ComplexReflectionGroup will return the
// exact same model as the smae command in GAP3. We load the model from the
// database.
//==============================================================================
intrinsic ComplexReflectionGroup(n:: RngIntElt : Model:="CHEVIE") -> GrpMat
{Returns an explicit model of the Shephard–Todd exceptional complex reflection group G_n. The default model is "CHEVIE". The Model "LT" (Lehrer-Taylor) is the one used in Magma.}

	name := ComplexReflectionGroupName(n)*"_"*Model;

	if Model eq "CHEVIE" then
		G := CHAMP_GetFromDB("ReflectionGroups", name);
		G`GAP3_Code := Sprintf("ComplexReflectionGroup(%o);", n);
	elif Model eq "LT" then
		G := ShephardTodd(n);
		G`DBName := name;
	else
		G := CHAMP_GetFromDB("ReflectionGroups", name);
	end if;

	G`IsReflectionGroup := true;

	return G;

end intrinsic;

intrinsic ComplexReflectionGroup(m::RngIntElt, p::RngIntElt, n::RngIntElt : Model:="CHEVIE") -> GrpMat
{Returns an explicit model of an irreducible reflection representation of the Shephard-Todd complex reflection group G(m,p,n). The default model is CHEVIE. Note that for m=p=1, this group is the symmetric group and the irreducible reflection representation is of dimension n-1; otherwise it is of dimension n. This differs from Magma's ShephardTodd(1,1,n).}


	name := ComplexReflectionGroupName(m,p,n)*"_"*Model;

	if Model eq "CHEVIE" then
		G := CHAMP_GetFromDB("ReflectionGroups", name);
		G`GAP3_Code := Sprintf("ComplexReflectionGroup(%o,%o,%o);", m,p,n);
	elif Model eq "LT" then
		if m eq 1 and p eq 1 then
			//I could take the quotient by the line but I don't want to implement
			//this now.
			error "LT has no irreducible model for the symmetric group";
		end if;
		G := ShephardTodd(n);
		G`DBName := name;
	else
		G := CHAMP_GetFromDB("ReflectionGroups", name);
	end if;

	G`IsReflectionGroup := true;

	return G;

end intrinsic;

//==============================================================================
// Shortcut functions
intrinsic TypeAReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat
{}
	return ComplexReflectionGroup(1,1,n+1:Model:=Model);
end intrinsic;

intrinsic SymmetricReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat
{}
	return ComplexReflectionGroup(1,1,n:Model:=Model);
end intrinsic;

intrinsic TypeBReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat
{}
	return ComplexReflectionGroup(2,1,n:Model:=Model);
end intrinsic;

intrinsic TypeDReflectionGroup(n::RngIntElt : Model:="CHEVIE") -> GrpMat
{}
	return ComplexReflectionGroup(2,2,n:Model:=Model);
end intrinsic;

intrinsic DihedralReflectionGroup(m::RngIntElt : Model:="CHEVIE") -> GrpMat
{}
	return ComplexReflectionGroup(m,m,2:Model:=Model);
end intrinsic;

intrinsic CyclicReflectionGroup(m::RngIntElt : Model:="CHEVIE") -> GrpMat
{}
	return ComplexReflectionGroup(m,1,1:Model:=Model);
end intrinsic;
