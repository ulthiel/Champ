//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Attributes for category Grp.
//
//==============================================================================

declare attributes Grp:
	GAP3_Code, //GAP3 code to construct the group
	NumberOfClasses, //number of conjugacy classes
	ClassNames,
	ClassWords,	//sequences of words of representatives of the conjugacy classes
	ClassNames, //names of conjugacy classes
	ClassLengths, //lengths of conjugacy classes
	ClassOrders, //orders of conjugacy classes
	CharacterNames, //Names of the irreducible complex characters.
	Classes, //Sequence containing the conjugacy classes of the group (in Magma format).
	ClassMap, //Maps a group element to the conjugacy class number of this element.
	CharacterTable, //Table of the irreducible complex characters.
	CharacterRing, //Complex character ring (Magma type).
	DBDir,
	DBFile,
	DBName
	;
