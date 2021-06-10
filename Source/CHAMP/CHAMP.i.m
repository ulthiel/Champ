//##############################################################################
// CHAMP (CHerednik Algebra Magma Package)
// Copyright (C) 2013–2021 Ulrich Thiel
// Licensed under GNU GPLv3, see License.txt.
// https://github.com/ulthiel/champ
// thiel@mathematik.uni-kl.de
//
//
// Basic CHAMP environment functions.
//
//##############################################################################

//##############################################################################
intrinsic GetOSType() -> MonStgElt
{Returns a string specifying the opering system type Magma is running on (Unix or Win.)}

	return GetEnv("CHAMP_OS_TYPE");

end intrinsic;

//##############################################################################
intrinsic GetOS() -> MonStgElt
{Returns a string specifying the opering system Magma is running on. This is more specific than GetOSType.}

	return GetEnv("CHAMP_OS");

end intrinsic;

//##############################################################################
intrinsic GetOSVersion() -> MonStgElt
{Returns a string specifying the version of the opering system Magma is running on. This is more specific than GetOS.}

	return GetEnv("CHAMP_OS_VER");

end intrinsic;

//##############################################################################
intrinsic GetCPU() -> MonStgElt
{Returns a string specifying the CPU Magma is running on.}

	return GetEnv("CHAMP_CPU");

end intrinsic;

//##############################################################################
intrinsic GetHostname() -> MonStgElt
{Returns a the name of the host Magma is running on.}

	return GetEnv("CHAMP_HOSTNAME");

end intrinsic;

//##############################################################################
intrinsic CHAMP_GetDir() -> MonStgElt
{Returns the base directory of CHAMP.}

	return GetEnv("CHAMP_DIR");

end intrinsic;

//##############################################################################
intrinsic CHAMP_GetVersion() -> MonStgElt
{Returns the CHAMP version.}

	return GetEnv("CHAMP_VER");

end intrinsic;

//##############################################################################
intrinsic GetOSArch() -> MonStgElt
{The operating system architecture.}

  return GetEnv("CHAMP_OS_ARCH");

end intrinsic;
