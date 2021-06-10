//==============================================================================
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see COPYING.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Some operating system functions.
//
//==============================================================================

//============================================================================
intrinsic Sleep(n::RngIntElt)
{Sleeps for n seconds.}

	System("sleep "*Sprint(n));

end intrinsic;

//============================================================================
intrinsic Date( : Unix:=false) -> MonStgElt
{The current date and time.}

	if GetOSType() eq "Win" then
		if Unix eq false then
			pipe := POpen(CHAMP_GetDir()*"Win\\date.exe", "r");
		else
			pipe := POpen(CHAMP_GetDir()*"Win\\date.exe +%s", "r");
		end if;
	else
		if Unix eq false then
			pipe := POpen("date", "r");
		else
			pipe := POpen("date +%s", "r");
		end if;
	end if;
	return Gets(pipe);

end intrinsic;

//============================================================================
intrinsic ShortDate() -> MonStgElt
{}

	if GetOSType() eq "Win" then
		pipe := POpen(CHAMP_GetDir()*"Win\\date.exe +%Y%m%d%H%M%S", "r");
	else
		pipe := POpen("date +%Y%m%d%H%M%S", "r");
	end if;
	return Gets(pipe);

end intrinsic;
