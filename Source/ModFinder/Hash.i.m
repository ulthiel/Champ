/*
  CHAMP (CHerednik Algebra Magma Package)
  Copyright (C) 2013–2020 Ulrich Thiel
  Licensed under GNU GPLv3, see COPYING.
  https://github.com/ulthiel/champ
  thiel@mathematik.uni-kl.de
*/

//==============================================================================
intrinsic LongHash(X::.) -> MonStgElt
{}
    hash := Sprint(ExtendedType(X));
    T := Type(X);
    if T eq AlgMatElt then
        hash *:= Sprint(BaseRing(X), "Magma");
        hash *:= Sprint(Nrows(X));
        hash *:= Sprint(Ncols(X));
        for i:=1 to Nrows(X) do
            hash *:= Sprint(Hash(X[i]));
        end for;
    elif T eq MtrxSprs then
        hash *:= LongHash(Matrix(X));
    elif T eq ModMatFldElt then
        hash *:= Sprint(BaseRing(X), "Magma");
        hash *:= Sprint(Nrows(X));
        hash *:= Sprint(Ncols(X));
        for i:=1 to Nrows(X) do
            hash *:= Sprint(Hash(X[i]));
        end for;
    elif T eq ModGr then
        hash *:= Sprint(BaseRing(X), "Magma");
        hash *:= Sprint(Dimension(X));
        hash *:= Sprint(X`MatrixDegrees);
        hash *:= Sprint(X`RowDegrees);
        hash *:= Sprint(X`Rep);
        for i:=1 to #X`Matrices do
            hash *:= LongHash(X`Matrices[i]);
        end for;
    elif T eq ModGrp then
        hash *:= Sprint(BaseRing(X), "Magma");
        hash *:= Sprint(Dimension(X));
        hash *:= Sprint(Group(X));
        for i:=1 to #ActionGenerators(X) do
            hash *:= LongHash(ActionGenerator(X,i));
        end for;
    elif T eq GrpMat then
        hash *:= Sprint(BaseRing(X), "Magma");
        hash *:= Sprint(Dimension(X));
        for i:=1 to Ngens(X) do
            hash *:= Sprint(LongHash(Matrix(X.i)));
        end for;
    else
        error "Not implemented for this type of object.";
    end if;

    return hash;

end intrinsic;

//==============================================================================
intrinsic MD5(str::MonStgElt) -> MonStgElt
{}
    os := GetOS();
    if os eq "Darwin" then
        cmd := "md5";
    elif os eq "Linux" then
        cmd := "md5sum";
    end if;
    return Pipe(cmd, str)[1..32]*"-MD5V1";

end intrinsic;

//==============================================================================
intrinsic MD5(X::.) -> MonStgElt
{}

    return MD5(LongHash(X));

end intrinsic;
