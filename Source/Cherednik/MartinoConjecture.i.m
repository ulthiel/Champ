/*
	CHAMP (CHerednik Algebra Magma Package)
	Copyright (C) 2010-2021 Ulrich Thiel
	Licensed under GNU GPLv3, see COPYING.
	thiel@mathematik.uni-kl.de
	https://ulthiel.com/math
*/

intrinsic IsRefinement(fam1::Setq, fam2::Setq) -> BoolElt
{Returns true if fam1 is a refinement of fam2.}

  for f in fam2 do
    X := { g : g in fam1 | g subset f };
    Y := {};
    for g in X do
      Y join:=g;
    end for;
    if Y ne f then
      return false;
    end if;
  end for;

  return true;

end intrinsic;


intrinsic MartinoConjecture(G::GrpMat) -> BoolElt
{}

  cm := CalogeroMoserFamilies(G);
  rou := RouquierFamilies(G);

  P := Universe(Keys(cm));
  Q := Universe(Keys(rou));

  cmhyp := SetToSequence(Keys(cm));
  cmhypideals := [ ideal<P|h> : h in cmhyp ];
  rouhyp := [ P!h : h in Keys(rou) ];
  rouhypideals := [ ideal<P|h> : h in rouhyp ];

  //First, check the generic families
  if cm[1] eq rou[1] then
    print "Generic Rouquier families are equal to generic Calogero-Moser families";
  else
    if IsRefinement(rou[1], cm[1]) then
      print "Generic Rouquier families refine the generic Calogero-Moser families";
    else
      print "Generic Rouquier families do not refine the generic Calogero-Moser families!";
      return false;
    end if;
  end if;

  //For each CM hyperplane determine the Rouquier families for the
  //corresponding sharp parameter and check if Rouquier families refine
  //CM families.
  //To find the correct hyperplane we need to use ideals.
  for i:=1 to #cmhyp do
    h := cmhyp[i];
    cmfams := cm[h];
    hsharpideal := ideal<P | P!G`MartinoSharp(h)>;
    j := Position(rouhypideals, hsharpideal);
    if j eq 0 then
      hrou := P!1;
    else
      hrou := rouhyp[j];
    end if;
    roufams := rou[hrou];
    if not IsRefinement(roufams, cmfams) then
      return false;
    end if;
  end for;

  return true;

  //
  // rouex := { P!NormalizeRationalHyperplaneEquation(P!G`MartinoSharp(H)) : H in Keys(rou) };
  // if rouex ne cmhyp then
  //   rouexeqcmgen := false;
  // end if;
  // if not rouex subset cmhyp then
  //   rouexcontainedincmgen := false;
  // end if;
  // for H in Keys(rou) do
  //   roufams := rou[H];
  //   Hsharp := P!NormalizeRationalHyperplaneEquation(P!G`MartinoSharp(H));
  //   if Hsharp notin cmhyp then
  //     cmfams := cm[1]`CMFamilies;
  //     rouexcontainedincmgen := false;
  //     rouexeqcmgen := false;
  //   else
  //     cmfams := cm[Hsharp]`CMFamilies;
  //   end if;
  //   if roufams ne cmfams then
  //     roueqcm := false;
  //     //check if each cm family is union of rouquier families (martino's conjecture)
  //     for cmfam in cmfams do
  //       subroufams := { F : F in roufams | F subset cmfam };
  //       if not { i : i in F, F in subroufams} eq cmfam then
  //         roufinercm := false;
  //       end if;
  //     end for;
  //   end if;
  // end for;
  // for H in cmhyp do
  //   cmfams := cm[H]`CMFamilies;
  //   Hsharp := NormalizeRationalHyperplaneEquation(Q!G`MartinoSharp(H));
  //   if Hsharp notin Keys(rou) then
  //     roufams := rou[1];
  //   else
  //     roufams := rou[Hsharp];
  //   end if;
  //   if roufams ne cmfams then
  //     roueqcm := false;
  //     //check if each cm family is union of rouquier families (martino's conjecture)
  //     for cmfam in cmfams do
  //       subroufams := { F : F in roufams | F subset cmfam };
  //       if not { i : i in F, F in subroufams} eq cmfam then
  //         roufinercm := false;
  //       end if;
  //     end for;
  //   end if;
  // end for;
  //
  // return roufinercm;

end intrinsic;
