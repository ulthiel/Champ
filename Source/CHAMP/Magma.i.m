//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Some Magma specific functions.
//
//==============================================================================

//============================================================================
intrinsic StartProfile()
{Starts profiling.}

	ProfileReset();
	SetProfile(true);

end intrinsic;

//============================================================================
intrinsic StopProfile() -> GrphDir
{Stops profiling and returns the profile graph.}

	SetProfile(false);
	return ProfileGraph();

end intrinsic;

//============================================================================
intrinsic HTMLProfile(P::GrphDir)
{}

	date := ShortDate();
	dir := CHAMP_GetDir()*"Profiles/"*date;
	ret := System("mkdir -p "*dir);
	ProfileHTMLOutput(P, "profile" : Dir:=dir);
	ret := System("cp "*CHAMP_GetDir()*"Profiles/sorttable.js "*dir);
	ret := System("cp "*CHAMP_GetDir()*"Profiles/leaf.gif "*dir);

end intrinsic;


//============================================================================
intrinsic HumanReadableTime(t::RngIntElt) -> MonStgElt
{Prints a time in seconds in human readable format.}

	return HumanReadableTime(t*1.0);

end intrinsic;

//============================================================================
intrinsic HumanReadableTime(t::FldReElt) -> MonStgElt
{Prints a time in seconds in human readable format.}
	F := RealField(2);
	t := F!t;
	if t lt 60 then
		return Sprint(t)*"s";
	elif t lt 3600 then
		return Sprint(t/60)*"m";
	elif t lt 86400 then
		return Sprint(t/3600)*"h";
	elif t lt 604800 then
		return Sprint(t/86400)*"d";
	elif t lt 2629800 then
		return Sprint(t/604800)*"w";
	elif t lt 31557600 then
		return Sprint(t/2629800)*"M";
	else
		return Sprint(t)*"s";
	end if;

end intrinsic;

//===========================================================================
intrinsic PrintAndDelete(str::MonStgElt)
{}

	delstr := "";
	for i:=1 to #str do
		delstr *:="\b";
	end for;

	printf str;
	printf delstr;

end intrinsic;

//===========================================================================
intrinsic GetPercentage(P::RngIntElt, G::RngIntElt) -> MonStgElt
{P/G*100 with four decimal digits as string.}

	p := ChangePrecision(RealField()!(P/G*100), 4);
	return Sprint(p)*"%";

end intrinsic;

//===========================================================================
intrinsic PrintPercentage(P::RngIntElt, G::RngIntElt)
{Prints P/G*100 with four decimal digits and puts back cursor so that percentage can be overwritten.}

	p := ChangePrecision(RealField()!(P/G*100), 4);
	pstr := Sprint(p);
	pstr := pstr[1..Minimum({5,#pstr})];
	printf "%6.2o%%", pstr;
	printf "\b\b\b\b\b\b\b";
end intrinsic;

//============================================================================
intrinsic GetVersionString() -> MonStgElt
{Returns the Magma version as a string.}

	a,b,c := GetVersion();

	return Sprint(a)*"."*Sprint(b)*"-"*Sprint(c);

end intrinsic;


//============================================================================
intrinsic GetMemoryUsageInMB() -> RngIntElt
{Returns the memory used by Magma in mega bytes.}

	return Round( GetMemoryUsage()/1024^2 );

end intrinsic;



//============================================================================
intrinsic PrintAndDelete(s::MonStgElt)
{Print s and set back cursor to the beginning of s so that s can be overwritten.}
	printf s;
	delstring := "";
	for i:=1 to #s do
		delstring *:="\b";
	end for;
	printf delstring;

end intrinsic;
