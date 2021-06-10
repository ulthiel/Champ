/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2013-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

declare verbose RngInvar, 5;

//============================================================================
intrinsic CHAMP_GetOS() -> MonStgElt
{}

	return GetEnv("CHAMP_OS");

end intrinsic;

//============================================================================
intrinsic CHAMP_GetDir() -> MonStgElt
{}

    __CHAMP_DIR := GetEnv("CHAMP_DIR");

		//attach trainling slash
    if CHAMP_GetOS() eq "Win" then
    	if __CHAMP_DIR[#__CHAMP_DIR] ne "\\" then
       		__CHAMP_DIR *:= "\\";
    	end if;
    else
   		if __CHAMP_DIR[#__CHAMP_DIR] ne "/" then
       		__CHAMP_DIR *:= "/";
    	end if;
    end if;

    return __CHAMP_DIR;

end intrinsic;

//============================================================================
intrinsic CHAMP_GetVersion() -> MonStgElt
{}

	__CHAMP_DIR := CHAMP_GetDir();

	if CHAMP_GetOS() eq "Unix" then
			ret1 := System("test -d "*__CHAMP_DIR*".git")/256;
			ret2 := System("git --version > /dev/null");
			if ret1 eq 0 and ret2 eq 0 then
				try
					__CHAMP_VER := Pipe("git --git-dir="*__CHAMP_DIR*".git describe", "");
					__CHAMP_VER := __CHAMP_VER[1..#__CHAMP_VER-1];
				catch e
					;
				end try;
			end if;
	end if;
	if assigned __CHAMP_VER then
		return __CHAMP_VER;
	else
		return "Unknown";
	end if;
end intrinsic;
