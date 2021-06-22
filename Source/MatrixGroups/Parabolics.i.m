/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

/*
	A parabolic subgroup of a matrix group (G,V) is the stabilizer G_W of a subspace W of V.
*/

//============================================================================
declare attributes GrpMat:
	Subgroups, //the list of conjugacy classes of subgroups
	ParabolicSubgroups; //the list of conjugacy classes of parabolic subgroups

//============================================================================
intrinsic Subgroups(~G::GrpMat)
{
	Attaches to G`Subgroups the list of conjugacy classes of subgroups of G.
}

	if assigned G`Subgroups then
		return;
	end if;

	G`Subgroups := Subgroups(G);

end intrinsic;

//============================================================================
intrinsic FixedSpace(G::GrpMat) -> ModTupFld
{
	If G is a subgroup of GL(V), return the subspace W of V fixed by G.
}

	V := VectorSpace(G);
	for g in Generators(G) do
		V meet:=FixedSpace(g);
	end for;

	return V;

end intrinsic;

//============================================================================
intrinsic PointwiseStabilizer(G::GrpMat, U::ModTupFld) -> GrpMat
{
	The subgroup of G fixing U pointwise.
}

	P := G;
	for u in Basis(U) do
		P meet:=Stabilizer(G, u);
	end for;
	return P;

end intrinsic;

//============================================================================
intrinsic IsParabolic(G::GrpMat, H::GrpMat) -> BoolElt
{
	Returns true iff H is a parabolic subgroup of G, i.e., iff H is the pointwise stabilizer of its fixed space.
}

	if H ne PointwiseStabilizer(G, FixedSpace(H)) then
		return false;
	end if;

	return true;

end intrinsic;

//============================================================================
intrinsic ParabolicSubgroups(~G::GrpMat)
{
	Attaches the list of conjugacy classes of parabolic subgroups.
}

	if assigned G`ParabolicSubgroups then
		return;
	end if;

	Subgroups(~G);
	G`ParabolicSubgroups := {};

	for i:=1 to #G`Subgroups do
		if IsParabolic(G, G`Subgroups[i]`subgroup) then
			G`ParabolicSubgroups join:={i};
		end if;
	end for;

end intrinsic;
