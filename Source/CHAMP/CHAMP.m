/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/


print "#########################################################";
print "#  CHAMP (CHerednik Algebra Magma Package)              #";
__CHAMP_VER := CHAMP_GetVersion();
__CHAMP_VER *:= &*[" " : i in [1..45-#__CHAMP_VER] ];
print "#  Version "*__CHAMP_VER*"#";
print "#  Copyright (C) 2013-2021 Ulrich Thiel                 #";
print "#  Licensed under GNU GPLv3                             #";
print "#  Please cite                                          #";
print "#    * LMS J. Comput. Math. 18 (2015), no. 1, 266-307   #";
print "#  Contributions by:                                    #";
print "#    * Cedric Bonnafe (Montpellier)                     #";
print "#    * Monika Truong (Stuttgart)                        #";
print "#  thiel@mathematik.uni-kl.de                           #";
print "#  https://github.com/ulthiel/champ                     #";
print "#########################################################";

//============================================================================
/*
	Environment
*/
//SetAutoColumns(true);

SetColumns(0);
