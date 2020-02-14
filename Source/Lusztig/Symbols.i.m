/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/



//============================================================================
intrinsic LusztigSymbolB(lambda::SeqEnum[SeqEnum[RngIntElt]], weight::SeqEnum, N::RngIntElt) -> SeqEnum
{
	Lusztig's symbol (as defined in [Lu03, 22]) for the bipartition lambda of n at parameters (b,a) with b=a-a-...-a for type Bn with respect to N.
}
	n := &+lambda[1] + &+lambda[2];

	b := weight[1];
	a := weight[2];

	r := b mod a;
	m := b div a;

	if not (N+m ge #lambda[1] and N ge #lambda[2]) then error "N not large enough.";
	end if;

	S1 := [];
	S2 := [];

	//do padding
	while #lambda[1] lt N+m do
		Append(~lambda[1], [0]);
	end while;

	while #lambda[2] lt N do
		Append(~lambda[2], [0]);
	end while;

	for i:=1 to N+m
		do Append(~S1, a*(lambda[1][N+m-i+1] + i - 1) + r);
	end for;

	for j:=1 to N do
		Append(~S2, a*(lambda[2][N-j+1] + j - 1));
	end for;

	return [S1,S2];

end intrinsic;

//============================================================================
intrinsic LusztigSymbolB(lambda::SeqEnum[SeqEnum[RngIntElt]], weight::SeqEnum[RngIntElt]) ->
SeqEnum
{
	Minimal representative of LusztigSymbolB.
}

	b := weight[1];
	a := weight[2];

	r := b mod a;
	m := b div a;

	N := Maximum({#lambda[1]+m,#lambda[2]}) + 1;
	return LusztigSymbolBReduce(LusztigSymbolB(lambda, weight, N), weight);

end intrinsic;

//============================================================================
intrinsic LusztigSymbolBShift(S::SeqEnum[SeqEnum[RngIntElt]], weight::SeqEnum[RngIntElt], s::RngIntElt) -> SeqEnum
{
    Shift by d of a Lusztig symbol for type B and given weight.

    Weight is (b,a) with b=a-a-...-a.
}


	b := weight[1];
	a := weight[2];

	r := b mod a;
	m := b div a;

	if s ge 0 then
		for z:=1 to s do
			S[1] := [r] cat [ S[1][i]+a : i in [1..#S[1]] ];
			S[2] := [0] cat [ S[2][i]+a : i in [1..#S[2]] ];
		end for;
	else
		for z:=1 to -s do
			if S[1][1] ne r then
				error "Cannot shift further down.";
			end if;
			S[1] := [ S[1][i]-a : i in [2..#S[1]] ];
			S[2] := [ S[2][i]-a : i in [2..#S[2]] ];
		end for;
	end if;

	return S;

end intrinsic;

//============================================================================
intrinsic LusztigSymbolBReduce(S::SeqEnum, weight::SeqEnum[RngIntElt]) -> SeqEnum
{
    Shift symbol as far down as possible. This is then a minima representative in the limit.
}

	b := weight[1];
	a := weight[2];

	r := b mod a;
	m := b div a;

	while S[1][1] eq r and #S[2] gt 0 and S[2][1] eq 0 do
		S := LusztigSymbolBShift(S, weight, -1);
	end while;

	return S;

end intrinsic;

//============================================================================
intrinsic LusztigFamiliesB(n::RngIntElt, weight::SeqEnum[RngIntElt] : Display:=true) -> SetEnum
{
The Lusztig families for Bn with weights (b,a) where b=a-a-...-a.
}

	b := weight[1];
	a := weight[2];

	r := b mod a;
	m := b div a;

	N := n+m;	//is sufficiently large I think

	parts := [ [lambda[1], lambda[2]] : lambda in Multipartitions(n,2)];

	symbols := [ LusztigSymbolB(lambda,weight,N) : lambda in parts ];

	contents := [ Content(S) : S in symbols ];

	cts := Sort(SetToSequence(SequenceToSet(contents)));

	if Display then
		str := "<h1>Lusztig families for type B"*Sprint(n)*" and weights "*Sprint(weight)*"</h1>";

		file := CHAMP_GetDir()*"HTML/Outputs/"*Tempname("")*".html";

		str *:= "There are "*Sprint(#cts)*" families.<br><br>";

		count := 0;
		for c in cts do
			count +:= 1;
			str *:= "<h3>Family "*Sprint(count)*"</h3>";
			fam := { i : i in [1..#symbols] | contents[i] eq c};
			S := Minimum({LusztigSymbolB(parts[i],weight) : i in fam});
			C := Content(S);
			str *:= "Symbol:<br>"*Sprint(S[1])*"<br>"*Sprint(S[2])*"<br><br>";
			str *:= "Content: "*HTML(C)*"<br><br>";
			str *:= "Representations:<br>\n";
			//str *:= "<table border=\"0\"><tr>";
			for i in fam do
				lambda := parts[i];
				symbol := symbols[i];
				str *:= Sprint(lambda)*"\t with symbol "*Sprint(symbol)*"<br>";
			end for;
			//str *:= "</tr></table>";
			str *:= "<br><br>";
		end for;

		WriteHTMLStart(file);
		Write(file, str);
		WriteHTMLEnd(file);
		ShowHTML(file);

	end if;

	return 0;

end intrinsic;

//============================================================================
intrinsic LusztigSymbolBOverline(lambda::SeqEnum[SeqEnum[RngIntElt]], t::RngIntElt) -> SeqEnum
{
The overline of a symbol lambda as defined by Lusztig in [Lu03, 22.8]. We then have lambda otimes sgn = overline of lambda.
}

	assert t ge Maximum(SequenceToSet(lambda[1]) join SequenceToSet(lambda[2]));

	lambdabar := [ SetToSequence({0..t} diff {t-lambda[2][i] : i in [1..#lambda[2]]}) , SetToSequence({0..t} diff {t-lambda[1][i] : i in [1..#lambda[1]]}) ];

	return lambdabar;

end intrinsic;
