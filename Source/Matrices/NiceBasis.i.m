/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

//============================================================================
intrinsic FindNiceBasis(Mats::SeqEnum[Mtrx]) -> Mtrx
{Tries to find a basis in which the matrices in X have many zeros.}


	K := BaseRing(Mats[1]);
	n := Ncols(Mats[1]);

	entries := {};
	for M in Mats do
		entries join:=SequenceToSet(Eltseq(M));
	end for;

	zeros := &+[ NumberOfNonZeroEntries(M) : M in Mats ];

	totalentries := #Mats*n^2;

	print "Matrices have "*Sprint(zeros)*" zeros ("*Sprint(ChangePrecision(RealField()!(zeros/totalentries*100), 4))*"%)";

	print "";

	goodT := IdentityMatrix(K, n);

	for i:=1 to #Mats do
		//take random invertible matrix for conjugating
		J,T := JordanForm(Mats[i]);
		MatsT := [ T*M*T^-1 : M in Mats ];
		newzeros := &+[ NumberOfNonZeroEntries(M) : M in MatsT ];
		if newzeros gt zeros then
			print "Found basis with "*Sprint(newzeros)*" zeros ("*Sprint(ChangePrecision(RealField()!(newzeros/totalentries*100), 4))*"%)";
			goodT := T;
			zeros := newzeros;
			IndentPush();
			print T;
			IndentPop();
			print "";
		end if;
	end for;

	return goodT;

end intrinsic;

//============================================================================
intrinsic FindNiceBasis(G::GrpMat) -> Mtrx
{}

	T:= FindNiceBasis([Matrix(g) : g in G]);

	return T, MatrixGroup<Dimension(G),BaseRing(G)| [ T*G.i*T^-1 : i in [1..Ngens(G)]]>;
end intrinsic;
