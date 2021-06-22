# 	CHAMP

A Cherednik Algebra Magma Package. By [Ulrich Thiel](https://ulthiel.com/math), 2013–2021.

## Scope

With this package you can:
* compute in rational Cherednik algebras (as introduced by [Etingof–Ginzburg](https://arxiv.org/abs/math/0011114));
* compute generators and a presentation of the center of the rational Cherednik algebra at t=0 (the coordinate algebra of the Calogero-Moser space);
* compute Poisson brackets;
* compute decomposition matrices of baby Verma modules and graded characters of simple modules for restricted rational Cherednik algebras (as introduced by [Gordon](https://arxiv.org/abs/math/0202301));
* compute Calogero–Moser families and hyperplanes;
* compute Calogero–Moser cellular characters (as introduced by [Bonnafé–Rouquier](https://arxiv.org/abs/1708.09764)).

The underlying reflection groups can be arbitrary and also the parameters can be arbitrary, including t≠0 and generic parameters valued in polynomial rings or rational function fields. An accompanying database contains many computational results. This document contains a complete overview of the functionality with many examples. The theory and algorithms are discussed in the following publications:
* Thiel, U. (2015). Champ: a Cherednik algebra Magma package. *LMS J. Comput. Math., 18*(1), 266–307. [https://doi.org/10.1112/S1461157015000054](https://doi.org/10.1112/S1461157015000054)
* Bonnafé, C. & Thiel, U. (2021). Calogero–Moser families and cellular characters: computational aspects.
* Thiel, U. (2014). A counter-example to Martino's conjecture about generic Calogero-Moser families. *Algebr. Represent. Theory, 17*(5), 1323–1348. [https://doi.org/10.1007/s10468-013-9449-4](https://doi.org/10.1007/s10468-013-9449-4)
* Thiel, U. (2017). Restricted rational Cherednik algebras. In: *Representation theory—current trends and perspectives* (pp. 681–745). Eur. Math. Soc., Zürich.

### Acknowledgements

I would like to thank Cédric Bonnafé for his contributions and endurance. Furthermore, I would like to thank Dario Mathiä for testing and feedback. 

## Contents

[1. Downloading an running](#downloading)  
[2. Complex reflection groups](#reflgroups)  
[3. Rational Cherednik algebras](#che)  
&nbsp;&nbsp;&nbsp;&nbsp;[3.1 Parameters](#params)  
&nbsp;&nbsp;&nbsp;&nbsp;[3.2 Rational Cherednik algebras at t=0 and Calogero–Moser spaces](#cmspaces)  
&nbsp;&nbsp;&nbsp;&nbsp;[3.3 Poisson brackets](#poisson-brackets)  
[4. Restricted rational Cherednik algebras](#rrca)  
[5 Representation theory of restricted rational Cherednik algebras](#rrca-rep)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.1 Conventions](#rrca-conv)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.2 Working with modules](#rrca-verma)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.3 Computing multiplicities](#rrca-mults)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.4 Calogero–Moser hyperplanes and families](#rrca-cmhyperplanes)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.5 Database](#rrca-db)  
[5. Calogero–Moser cellular characters](#rrca-cellular)

<a name="downloading"></a>

## Downloading and running

You need a [Magma](http://magma.maths.usyd.edu.au/magma/) version of at least 2.19 (current version is 2.26 and CHAMP is tested with 2.25) and the Magma executable needs to be in your system PATH so that you can start Magma by calling ``magma`` in the terminal. CHAMP has an accompanying database containing data about complex reflection groups and rational Cherednik algebras that is necessary for full functionality and that moreover contains many computational results. This database is stored as a [separate repository](https://gitlab.rhrk.uni-kl.de/ulthiel/champ-db) with [Git Large File Storage](https://git-lfs.github.com) on my university's Git server to prevent Github's bandwidth limit on LFS storage.

So, to get CHAMP and the database you can do *one* of the following:

1. More stable: Download the [latest release](https://github.com/ulthiel/champ/releases/latest) of CHAMP together with the database from the release assets.
2. More up to date: Download the [latest source of CHAMP](https://github.com/ulthiel/Champ/archive/refs/heads/master.zip) and the [latest source of the database](https://gitlab.rhrk.uni-kl.de/ulthiel/champ-db/-/archive/master/champ-db-master.zip) and place the database as the folder "Champ-DB" inside the CHAMP directory.
3. Best for developers: Clone the CHAMP repository via ``git clone https://github.com/ulthiel/champ`` and then inside the CHAMP directory clone the Git LFS database via ``git clone https://gitlab.rhrk.uni-kl.de/ulthiel/champ-db Champ-DB``. To this end, you need to install the Git LFS extension first as described [here](https://git-lfs.github.com).

You can then run CHAMP via ```./champ```:

```c++
#########################################################
#  CHAMP (CHerednik Algebra Magma Package)              #
#  Version v1.6.0-beta-2-g009e527                       #
#  Copyright (C) 2013-2021 Ulrich Thiel                 #
#  Licensed under GNU GPLv3                             #
#  Please cite                                          #
#    * LMS J. Comput. Math. 18 (2015), no. 1, 266-307   #
#  Contributions by:                                    #
#    * Cedric Bonnafe (Montpellier)                     #
#    * Monika Truong (Stuttgart)                        #
#  thiel@mathematik.uni-kl.de                           #
#  https://ulthiel.com/math                             #
#########################################################
> 
```

I advise to once run ```./selfcheck``` in the directory ```SelfCheck```. The ReflectionGroups and G5_Verma selfcheck will take a bit of time.

<a name="reflgroups"></a>

## Complex reflection groups

Models for several complex reflection groups, their character tables, character names, models for irreducible representations, etc. is stored in the ReflectionGroups database. The data is taken from (and compatible with) J. Michel's [CHEVIE](https://webusers.imj-prg.fr/~jean.michel/chevie/chevie.html) package from 2014 (there were some character label changes afterwards but this is not dramatic; everything in CHAMP is consistent). The reason for using a database is that we need consistent labelings (of e.g. characters) that allow us to compare results with the literature. A general philosophy in CHAMP is that most objects (like groups) will have attributes (like CharacterTable) which are set by a similarly named procedure operating on the object (using the ~ operator). Usually, it is first checked whether the data exists in the database; if not, it will be computed in a consistent way.

The following examples demonstrate how to use all functions around complex reflection groups:

```C++
//Load the Weyl group B2 in a reflection representation
> W := TypeBReflectionGroup(2);
> W;
MatrixGroup(2, Rational Field)
Generators:
    [-1  2]
    [ 0  1]

    [ 1  0]
    [ 1 -1]

//The database location for this group is stored in the DBDir attribute
> W`DBDir;
ReflectionGroups/B2_CHEVIE/

//Character tables and standard character names are stored in the database.
> CharacterTable(~W);
> W`CharacterTable;
[
    ( 1, 1, 1, -1, -1 ),
    ( 2, -2, 0, 0, 0 ),
    ( 1, 1, -1, -1, 1 ),
    ( 1, 1, 1, 1, 1 ),
    ( 1, 1, -1, 1, -1 )
]
> W`CharacterNames;
[ 11., 1.1, .11, 2., .2 ] //notation for bi-partitions

//IMPORTANT: CharacterTable(W) without the ~ will use Magma's algorithm to
//compute the character table; we won't get a labeling! Hence, always use the
//procedure with the ~ operator.

//Load models for the irreducible representations. Their numbering will match
//the one from the database.
> Representations(~W);
> W`Representations[0]; //I wanted to use positive characteristic
                        //representations one day, hence the 0.

//Fake degrees (graded W-character of the coinvariant algebra)
> FakeDegrees(~W);
> W`FakeDegrees;
[
  q^2,
  q^3 + q,
  q^4,
  1,
  q^2
]

//Other types of reflection groups (with connection to data from the database
//and/or natural choices) can be created with the functions listed below.
//You can also load some special models directly from the database as in the
//following example which loads a particular model of B2 used by Bonnafé-
//Rouquier in some computation
> W := CHAMP_GetFromDB("ReflectionGroups/B2_BR", "GrpMat");
> W;
MatrixGroup(2, Rational Field)
Generators:
    [0 1]
    [1 0]

    [-1  0]
    [ 0  1]
```

You can load all sorts of reflections groups with the following commands:

* ExceptionalComplexReflectionGroup (groups G4 to G37 in Shephard–Todd notation)
* SymmetricReflectionGroup
* TypeBReflectionGroup
* TypeDReflectionGroup
* DihedralReflectionGroup
* CyclicReflectionGroup
* ImprimitiveReflectionGroup (groups G(m,p,n) in Shephard–Todd notation)

Note: I have not imported all the data for all possible groups; my main focus were the exceptional groups. More data can always be added to the database of course.

<a name="che"></a>

## Rational Cherednik algebras

The definitition of rational Cherednik algebras used in CHAMP is exactly the one by [Etingof–Ginzburg](https://arxiv.org/abs/math/0011114). It's best to begin with an example straightaway.

```C++
//Create the rational Cherednik algebra for t and c generic (valued in a
//polynomial ring)
> W := TypeBReflectionGroup(2); //Weyl group of type B2 as above
> H := RationalCherednikAlgebra(W : Type:="EG"); //I will explain the EG below
> H;
Rational Cherednik algebra
Generators:
    w1, w2, y1, y2, x1, x2
Generator degrees:
    0, 0, -1, -1, 1, 1
Base ring:
    Polynomial ring of rank 3 over Rational Field
    Order: Lexicographical
    Variables: t, c1, c2
Group:
    MatrixGroup(2, Rational Field) of order 2^3
    Generators:
    [-1  2]
    [ 0  1]

    [ 1  0]
    [ 1 -1]
t-parameter:
    t
c-parameter:
    Mapping from: { 1 .. 2 } to Polynomial ring of rank 3 over Rational Field
    <1, c1>
    <2, c2>

//There is quite a bit to discuss now but let's start playing directly.
//As you can see in the output, there are generators w1, w2, y1, y2, x1, x2.
//These refer to the generators of the group (the w's), the basis of the space
//W is acting on (the y's) and its dual space (the x's). You can access the i-th
//generator in this numbering with H.i.
> H.3;
[1 0]
[0 1]*(y1)

//As a module, the Cherednik algebra is the group ring of W with coefficients
//in R[V \oplus V^*], where R is the base ring of the parameters. This is how
//algebra elements are represented also in CHAMP. Let's do some arithmetic.
> H.5*H.2;
[ 1  0]
[ 1 -1]*(x1 + x2)
> H.5*H.3;
[-1  0]
[-1  1]*(c2)
+
[-1  2]
[ 0  1]*(c1)
+
[1 0]
[0 1]*(y1*x1 + t)

//IMPORTANT: In Magma, matrices are acting from the *right* on vectors. Hence,
//to keep everything consistent, I have implemented the *opposite* of the
//rational Cherednik alebra as usually written on paper. This may be a bit
//confusing, but in the end it's less confusing than trying to artifically make 
//Magma act on the left.
```

The database contains (some) generators of the center of the rational Cherednik algebra at t=0 (see below). These elements are huge and the computation takes a lot of time, so one really wants to store them. I have implemented a function ```Rprint``` (for "reversible print") that returns program code for an element of an arbitrary Cherednik algebra allowing to reconstruct this element. Here's an example:

```
> W := TypeBReflectionGroup(2);
> H := RationalCherednikAlgebra(W : Type:="EG");
> Rprint(H.1);
/*
	Code for a Cherednik algebra element
	Version: v1.6.0-alpha-66-g295f039
	Date: Fri Jun 11 08:28:16 CEST 2021
*/
//base ring of the group
K := RationalField();
//the group
W := MatrixGroup<2, RationalField() | [ -1, 2, 0, 1 ] , [ 1, 0, 1, -1 ]>;
//the parameters of the Cherednik algebra
R := PolynomialRing(RationalField(), 3);
t := R.1;
c1 := R.2;
c2 := R.3;
t := t;
c := map<{ 1 .. 2 }-> R |[<1,c1>,<2,c2>]>;
//the Cherednik algebra
H := RationalCherednikAlgebra(W,<t,c>);
A := H`GroupAlgebra;
y1 := H`yxAlgebra.1;
y2 := H`yxAlgebra.2;
x1 := H`yxAlgebra.3;
x2 := H`yxAlgebra.4;
//the support of the Cherednik algebra element (as words in the generators)
supp := [[ 1 ]];
//the coefficients of the Cherednik algebra element (these are elements of S)
coeff1 := 1;
hA := Zero(A);
hA +:= coeff1*(A!WordToElement(W,supp[1]));
h := Zero(H); h`Element := hA;
return h
```

<a name="params"></a>

### Parameters

This topic is a bit technical but important. There are two kinds of parameters involved in the relations for the rational Cherednik algebra: a *t-parameter* and a *c-parameter*. Let's take a commutative ring R as base ring. The t-parameter is some fixed element of R; the c-parameter is a function c:Refl(W)/W → R from the conjugacy classes of reflections of W to R. For example, we can let R be a polynomial ring K[t,c<sub>1</sub>,...,c<sub>r</sub>] and define the parameters t and c in the obvious way. In this case we say the parameters are *generic*. If I is an ideal of R, we can also consider R/I as new base ring and get parameters with are *generic for the subscheme* defined by I. For example, we could take a polynomial ring R=K[t,c] and set c(s)=c for all c. This would be the generic *equal* parameter case.

For the construction of the rational Cherednik algebra in CHAMP you can take as base ring R any K-algebra that can be defined in Magma, where K is the base field of the reflection group W, and as parameters you can take any t and maps c with values in R. In particular, you can work with generic parameters, generic parameters on a, say, hyperplane, or special parameters taking values in your base field K. You have complete freedom.

[Ginzburg-Guay-Opdam-Rouquier](https://arxiv.org/abs/math/0212036) considered a Fourier transform on the c-parameter space which makes some expressions in the parameters much simpler (such as equations for the Calogero–Moser hyperplanes). I will refer to these as *k-parameters*. While the c-parameters by Etingof-Ginzburg are indexed by conjugacy classes of reflections, the k-parameters have a double index: the first indexes an orbit [H] of reflection hyperplanes, the second is an index between 0 and |W<sub>H</sub>|-1, where W<sub>H</sub> is the stabilizer of a representative of [H]. Of course, in the end the number of parameters is the same. By default, CHAMP uses k-parameters.

The following examples should make all of the above discussion clear.

```C++
//First, some shortcuts for creating generic rational Cherednik algebras:
> W:=TypeBReflectionGroup(2);
> H:=RationalCherednikAlgebra(W); //generic t and generic k-parameter
> H:=RationalCherednikAlgebra(W : Type:="EG"); //generic t and generic c
> H:=RationalCherednikAlgebra(W,0); //t=0 and generic k-parameter
> H:=RationalCherednikAlgebra(W,0 : Type:="EG"); //t=0 and generic c

//Now, let's have a closer look at parameters. Let's create a generic
//c-parameter.
> CherednikParameter(W : Type:="EG");
Mapping from: { 1 .. 2 } to Polynomial ring of rank 2 over Rational Field
    <1, c1>
    <2, c2>

//This is a map from (labels of) conjugacy classes of reflections of W to the
//polynomial ring in that many variables. Representatives of the conjugacy
//classes of reflections can be obtained as follows:
> W`ReflectionClasses;
[
    [-1  2]
    [ 0  1],

    [ 1  0]
    [ 1 -1]
]

//Let's construct the rational Cherednik algebra of W over the rational numbers
//with t=0 and a c-parameter with values c(1) = -1 and c(2) = 1:
> c := map<{1,2} -> Rationals() | [<1,-1>, <2,1>] >;
> H:=RationalCherednikAlgebra(W,0,c);

//Let's create a c-parameter which is generic for the hyperplane c_1 - c_2
//(this is the generic equal parameter case):
> c := CherednikParameter(W : Type:="EG");
> R:=Codomain(c);
> cH:= SpecializeCherednikParameterInHyperplane(c, R.1-R.2);
> cH;
Mapping from: { 1 .. 2 } to Multivariate rational function field of rank 1 over
Rational Field
    <1, c2>
    <2, c2>
> H := RationalCherednikAlgebra(W,0,cH);

//You can create a generic *rational* c-parameter as follows:
> CherednikParameter(W : Type:="EG", Rational:=true);
Mapping from: { 1 .. 2 } to Multivariate rational function field of rank 2 over
Rational Field
    <1, c1>
    <2, c2>

//Now, let's look at k-parameters (the default):
> k := CherednikParameter(W);
> k;
Mapping from: { 1 .. 2 } to Polynomial ring of rank 2 over Rational Field
    <1, 2*k1_1>
    <2, 2*k2_1>

//The labeling of orbits of reflection hyperplanes is consistent with what is
//stored in
> W`ReflectionLibrary;
//This is an array indexed by orbits of reflection hyperplanes. Each entry is
//again an array indexed by reflection hyperplanes in this orbit. The entries
//of this array are the reflections for the corresponding hyperplane.
//You can work with k-parameters exactly as with the c-parameters above.
```

<a name="cmspaces"></a>

## Rational Cherednik algebras at t=0 and Calogero–Moser spaces

The rational Cherednik algebra H<sub>t=0,c</sub> has a big center Z<sub>c</sub>. The center is a Poisson deformation of the symplectic singularity (V ⊕ V<sup>&ast;</sup>)/W, where W acts on V. The associated variety is called the *Calogero–Moser space* of W at parameter c. CHAMP can compute algebra generators of Z<sub>c</sub> and also a presentation of this algebra (the former works even for large groups like F<sub>4</sub>, the latter involves rather complicated invariant theory computations which are even for small dihedral groups too much; but you can still get some ideas).

The database contains generators of Z<sub>0</sub> (undeformed case) and Z<sub>k</sub> (k generic) for several cases. Some of the elements are extremely large (for G<sub>11</sub> there is one taking up >100MB compressed and >500MB uncompressed)! By default, all functions check the database first and load the data from there if available.


```C++
> W := TypeBReflectionGroup(2);
> H := RationalCherednikAlgebra(W,0);
> CenterGenerators(H); //this needs generic parameters!
[*
    [1 0]
    [0 1]*(y1^2 - 2*y1*y2 + 2*y2^2),
    [-1  0]
    [-1  1]*(k2_1)
    +
    [-1  2]
    [ 0  1]*(k1_1)
    +
    [ 1  0]
    [ 1 -1]*(k2_1)
    +
    [1 0]
    [0 1]*(y1*x1 + y2*x2)
    +
    [ 1 -2]
    [ 0 -1]*(k1_1),
    [1 0]
    [0 1]*(x1^2 + x1*x2 + 1/2*x2^2),
    [1 0]
    [0 1]*(y1^4 - 4*y1^3*y2 + 6*y1^2*y2^2 - 4*y1*y2^3 + 2*y2^4),
    [-1  0]
    [-1  1]*(k2_1*y1^2 - 4*k2_1*y1*y2 + 4*k2_1*y2^2)
    +
    [-1  2]
    [ 0  1]*(k1_1*y1^2 - 2*k1_1*y1*y2)
    +
    [ 1  0]
    [ 1 -1]*(k2_1*y1^2)
    +
    [1 0]
    [0 1]*(y1^3*x1 - 4*y1^2*y2*x1 - y1^2*y2*x2 + 4*y1*y2^2*x1 + 2*y1*y2^2*x2)
    +
    [ 1 -2]
    [ 0 -1]*(-k1_1*y1^2 + 2*k1_1*y1*y2),
    [-1  2]
    [ 0  1]*(-2*k1_1*y2*x1 - 2*k1_1*y2*x2)
    +
    [ 1 -2]
    [ 1 -1]*(-2*k1_1*k2_1)
    +
    [1 0]
    [0 1]*(y1^2*x1^2 + y1^2*x1*x2 + 1/2*y1^2*x2^2 - 4*y1*y2*x1^2 - 4*y1*y2*x1*x2
    - y1*y2*x2^2 + 4*y2^2*x1^2 + 4*y2^2*x1*x2 + y2^2*x2^2)
    +
    [-1  0]
    [ 0 -1]*(-2*k1_1^2)
    +
    [ 1 -2]
    [ 0 -1]*(-2*k1_1*y1*x1 + 2*k1_1*y2*x1)
    +
    [-1  2]
    [-1  1]*(-2*k1_1*k2_1),
    [-1  0]
    [-1  1]*(k2_1*x1^2 + k2_1*x1*x2 + 1/4*k2_1*x2^2)
    +
    [-1  2]
    [ 0  1]*(k1_1*x1^2 + 3/2*k1_1*x1*x2 + 3/4*k1_1*x2^2)
    +
    [ 1  0]
    [ 1 -1]*(1/4*k2_1*x2^2)
    +
    [1 0]
    [0 1]*(y1*x1^3 + 3/2*y1*x1^2*x2 + 3/4*y1*x1*x2^2 + 1/4*y2*x2^3)
    +
    [ 1 -2]
    [ 0 -1]*(k1_1*x1^2 + 1/2*k1_1*x1*x2 + 1/4*k1_1*x2^2),
    [1 0]
    [0 1]*(x1^4 + 2*x1^3*x2 + 3/2*x1^2*x2^2 + 1/2*x1*x2^3 + 1/8*x2^4)
*]
> #CenterGenerators(H);
8
//The computation of the center generators inductively deforms fundamental
//invariants of Z_0 = K[V \oplus V^*]^W. You can compute and acccess these
//fundamental invariants as follows:
> SymplecticDoublingFundamentalInvariants(W);
[
    y1^2 - 2*y1*y2 + 2*y2^2,
    y1*x1 + y2*x2,
    x1^2 + x1*x2 + 1/2*x2^2,
    y1^4 - 4*y1^3*y2 + 6*y1^2*y2^2 - 4*y1*y2^3 + 2*y2^4,
    y1^3*x1 - 4*y1^2*y2*x1 - y1^2*y2*x2 + 4*y1*y2^2*x1 + 2*y1*y2^2*x2,
    y1^2*x1^2 + y1^2*x1*x2 + 1/2*y1^2*x2^2 - 4*y1*y2*x1^2 - 4*y1*y2*x1*x2 -
        y1*y2*x2^2 + 4*y2^2*x1^2 + 4*y2^2*x1*x2 + y2^2*x2^2,
    y1*x1^3 + 3/2*y1*x1^2*x2 + 3/4*y1*x1*x2^2 + 1/4*y2*x2^3,
    x1^4 + 2*x1^3*x2 + 3/2*x1^2*x2^2 + 1/2*x1*x2^3 + 1/8*x2^4
]
//The deformation of an element of Z_0 to an element of Z_c is done with the
//function TruncationInverse which you can also call directly if you are only
//interested in special elements:
> TruncationInverse(H, W`SymplecticDoublingFundamentalInvariants[1]);
[1 0]
[0 1]*(y1^2 - 2*y1*y2 + 2*y2^2)
//On V \oplus V^* we have a natural N^2-grading. We are especially interested in
//algebra generators of N^2-degree (d,d), i.e. of Z-degree 0.
> [ Bidegree(f) : f in W`SymplecticDoublingFundamentalInvariants ];
[ <0, 2>, <1, 1>, <2, 0>, <0, 4>, <1, 3>, <2, 2>, <3, 1>, <4, 0> ]
//We see there are only 2 generators of Z-degree 0.
//You can also directly compute only the degree-0 generators of Z_c as follows
> CenterGeneratorsOfDegreeZero(H);
[*
    [-1  0]
    [-1  1]*(k2_1)
    +
    [-1  2]
    [ 0  1]*(k1_1)
    +
    [ 1  0]
    [ 1 -1]*(k2_1)
    +
    [1 0]
    [0 1]*(y1*x1 + y2*x2)
    +
    [ 1 -2]
    [ 0 -1]*(k1_1),
    [-1  2]
    [ 0  1]*(-2*k1_1*y2*x1 - 2*k1_1*y2*x2)
    +
    [ 1 -2]
    [ 1 -1]*(-2*k1_1*k2_1)
    +
    [1 0]
    [0 1]*(y1^2*x1^2 + y1^2*x1*x2 + 1/2*y1^2*x2^2 - 4*y1*y2*x1^2 - 4*y1*y2*x1*x2
    - y1*y2*x2^2 + 4*y2^2*x1^2 + 4*y2^2*x1*x2 + y2^2*x2^2)
    +
    [-1  0]
    [ 0 -1]*(-2*k1_1^2)
    +
    [ 1 -2]
    [ 0 -1]*(-2*k1_1*y1*x1 + 2*k1_1*y2*x1)
    +
    [-1  2]
    [-1  1]*(-2*k1_1*k2_1)
*]
//We can even compute a presentation of the center of H
> CenterPresentation(H);
[
    3*z1^2*z3 - z1*z2^2 - z1*z6 + 2*k1_1^2*z1 + z2*z5 - 2*z3*z4,
    -4*z1*z2*z3 + 2*z1*z7 + z2^3 + 2*z2*z6 - 4*k2_1^2*z2 - z3*z5,
    2*z1*z8 + z2^2*z3 - 2*z2*z7 - z3*z6 + 2*k1_1^2*z3,
    8*z1^3*z3 - 3*z1^2*z2^2 - 4*z1^2*z6 + (4*k1_1^2 + 8*k2_1^2)*z1^2 +
        2*z1*z2*z5 - 8*z1*z3*z4 + 4*z1*z3*z8 + 2*z2^2*z3^2 + 2*z2^2*z4 -
        4*z2*z3*z7 - 2*z3^2*z6 + 4*k1_1^2*z3^2 + 4*z4*z6 - 8*k2_1^2*z4 - z5^2,
    -7*z1^2*z2*z3 + 6*z1^2*z7 + z1*z2^3 + 3*z1*z2*z6 + (2*k1_1^2 -
        4*k2_1^2)*z1*z2 + 2*z2*z3*z4 - 4*z4*z7 - z5*z6 + 2*k1_1^2*z5,
    8*z1^2*z3^2 - 8*z1^2*z8 - 10*z1*z2^2*z3 + 6*z1*z2*z7 + (8*k1_1^2 -
        4*k2_1^2)*z1*z3 + 2*z2^4 + 3*z2^2*z6 + (-6*k1_1^2 - 8*k2_1^2)*z2^2 +
        z2*z3*z5 - 8*z3^2*z4 + 8*z4*z8 - 2*z5*z7 + (-4*k1_1^2 + 4*k2_1^2)*z6 +
        8*k1_1^4 - 8*k1_1^2*k2_1^2,
    -6*z1^2*z3^2 + 10*z1^2*z8 + 8*z1*z2^2*z3 - 8*z1*z2*z7 - z2^4 - 2*z2^2*z6 +
        (4*k1_1^2 + 4*k2_1^2)*z2^2 + 4*z3^2*z4 - 4*z4*z8 - z6^2 + 4*k1_1^2*z6 -
        4*k1_1^4,
    -4*z1*z2*z3^2 + 2*z1*z2*z8 + 4*z1*z3*z7 + 3*z2^3*z3 - 4*z2^2*z7 + z2*z3*z6 +
        (-2*k1_1^2 - 4*k2_1^2)*z2*z3 - 2*z3^2*z5 + 2*z5*z8 - 2*z6*z7 +
        4*k1_1^2*z7,
    -4*z1*z3^3 + 4*z1*z3*z8 - 2*z2^2*z3^2 - 2*z2^2*z8 + 8*z2*z3*z7 + 4*z3^2*z6 -
        4*k2_1^2*z3^2 - 4*z6*z8 - 4*z7^2 + 8*k2_1^2*z8
]
```

<a name="poisson-brackets"></a> 

The degree-0 center generators for several exceptional complex reflection groups are stored in the database and are loaded automatically when requested. Here's an example:

```
> W:=ExceptionalComplexReflectionGroup(28);
> H := RationalCherednikAlgebra(W,0);
> CenterGeneratorsOfDegreeZero(~H);
Fetched from DB
Deforming center generator 1 of 6
Found center generator in DB.
Deforming center generator 2 of 6
Found center generator in DB.
Deforming center generator 3 of 6
Found center generator in DB.
Deforming center generator 4 of 6
Found center generator in DB.
Deforming center generator 5 of 6
Found center generator in DB.
Deforming center generator 6 of 6
Found center generator in DB.
```

### Poisson brackets

You can compute Poisson brackets between elements in the Cherednik algebra.

```
> W := TypeBReflectionGroup(2);
> H := RationalCherednikAlgebra(W,0);
> PoissonBracket(H.5,H.3);
[-1  0]
[-1  1]*(2*k2_1)
+
[-1  2]
[ 0  1]*(2*k1_1)
+
[1 0]
[0 1]*(1)
```

<a name="rrca"></a>
## Restricted rational Cherednik algebras

The *restricted* rational Cherednik algebra is an important finite-dimensional quotient of the rational Cherednik algebra at t=0. See the paper by [Gordon](https://arxiv.org/abs/math/0202301) or [my paper](https://arxiv.org/abs/1603.05230). Computation in the restricted algebra can be done in CHAMP in the same way as with the uncrestricted algebra.

```C++
> W := TypeBReflectionGroup(2);
> H := RestrictedRationalCherednikAlgebra(W); //generic k-parameter
> H;
Restricted rational Cherednik algebra
Generators:
    w1, w2, y1, y2, x1, x2
Generator degrees:
    0, 0, -1, -1, 1, 1
Base ring:
    Multivariate rational function field of rank 2 over Rational Field
    Variables: k1_1, k2_1
Group:
    MatrixGroup(2, Rational Field) of order 2^3
    Generators:
    [-1  2]
    [ 0  1]

    [ 1  0]
    [ 1 -1]
c-parameter:
    Mapping from: { 1 .. 2 } to Multivariate rational function field of rank 2
    over Rational Field
    <1, 2*k1_1>
    <2, 2*k2_1>
//Here's one caveat: before you can do actual computations in the RRCA, you need to initialize it, which means here that the coinvariant algebra etc. is computed. This can be quite complex and not all of this is necessary when you are just interested in the representation theory, that's why I added an initialize function.
> Initialize(~H);
//Now, we're ready to do computations
> H.5*H.2;
[ 1  0]
[ 1 -1]*(x1 + x2)

//You can convert H into a matrix algebra
> A:=MatrixAlgebra(H);
> A;
Matrix Algebra of degree 512 with 6 generators over Multivariate rational
function field of rank 2 over Rational Field

//We compute the Jacobson radical for the equal parameter case k=[1,1]:
> k := CherednikParameter(W,[1,1]);
> H := RestrictedRationalCherednikAlgebra(W,k);
> Initialize(~H);
> A := MatrixAlgebra(H);
> time J := JacobsonRadical(A); J;
Time: 182.370
Matrix Algebra [ideal of A] of degree 512 and dimension 346 over Rational Field
```

<a name="rrca-rep"></a>

## Representation theory of restricted rational Cherednik algebras

In CHAMP you can compute baby Verma modules (also called standard modules) for restricted rational Cherednik algebras (as defined by [Gordon](https://arxiv.org/abs/math/0202301)). Using modular lifting techniques I introduced in [my paper](https://arxiv.org/abs/1403.6686) you can compute the heads of standard modules (which then give all the simples of the restricted rational Cherednik algebra) as graded modules (also giving the graded W-character) and the (graded) decomposition matrix of standard modules into simples. It works surprisingly well even in huge and complicated examples, and for generic parameters as well.

<a name="rrca-conv"></a>

### Conventions

Let W be a complex reflection group acting on a vector space V over a field K. Let K[V] be the symmetric algebra of V<sup>&ast;</sup>. In the (restricted) rational Cherednik algebra I am putting V<sup>*</sup> in degree +1, V in degree -1, and W in degree 0. This yields a triangular decomposition H = H<sup>-</sup> ⊗ KW ⊗ H<sup>+</sup>. The standard module Δ(λ) of an irreducible W-module λ is obtained by inflating λ to a (H<sup>-</sup> ⊗ KW)-module (i.e. V acting trivial) and then inducing it to an H-module. So, as a vector space, Δ(λ) = K[V]<sub>W</sub> ⊗ λ, where K[V]<sub>W</sub> is the coinvariant algebra. With my grading convention, Δ(λ) lives in *positive* degree.

Note that there are two choices: 1) to put V<sup>&ast;</sup> in degree +1; 2) to inflate λ to a (H^<sup>-</sup> ⊗ KW)-module. You could also put V<sup>&ast;</sup> in degree -1 and/or inflate λ to an (H<sup>+</sup> ⊗ KW)-module. Here is an overview of what is used in the literature:

| Paper | deg V<sup>&ast;</sup> | Δ(λ) |
| ----- | ----------------- | ---- |
| [CHAMP](https://arxiv.org/abs/1403.6686) | +1 | H ⊗<sub>H<sup>-</sup></sub> λ = K[V]<sub>W</sub> ⊗ λ |
| [Bonnafé-Roquier](https://arxiv.org/pdf/1708.09764.pdf) | +1 | H ⊗<sub>H<sup>-</sup></sub> λ = K[V]<sub>W</sub> ⊗ λ|
| [Bellamy-Thiel](https://arxiv.org/abs/1705.08024) | -1 | H ⊗<sub>H<sup>+</sup></sub> λ = K[V]<sub>W</sub> ⊗ λ|
| [Gordon](https://arxiv.org/abs/math/0202301) | -1 | H ⊗<sub>H<sup>-</sup></sub> λ = K[V<sup>&ast;</sup>]<sub>W</sub> ⊗ λ |

So, CHAMP and Bonnafé-Rouquier use the *same* conventions. The difference between Bonnafé-Roquier and Bellamy-Thiel is only an *opposite grading* on the Δ(λ) (up to the grading the modules are the same!). To make this more precise, consider a ℤ-graded algebra A with *triangular decomposition*, i.e. a triple (A<sup>l</sup>, A<sup>0</sup>, A<sup>r</sup>) of graded subalgebras such that the multiplication map A<sup>l</sup> ⊗ A<sup>0</sup> ⊗ A<sup>r</sup> → A is an isomorphism of vector spaces, and moreover the following holds: A<sup>0</sup> is in degree 0, and A<sup>l</sup> is either in positive or in negative degree, and A<sup>r</sup> is in the opposite degree of A<sup>l</sup>. In any case one can define the standard module Δ(λ) = A ⊗<sub>A<sup>r</sup></sub> λ. The inflation is always through the *right* part of the decomposition, so it is up to the grading independent of the aforementioned choice. In Bellamy-Thiel we assumed that A<sup>l</sup> is in negative degree, Bonnafé-Rouquier assume that it is in positive degree. But we both assume that A<sup>l</sup> = K[V]<sub>W</sub>. The Bonnafé-Rouquier assumption is nicer in the sense that the standard modules live in positive degree, which seems more natural (but it doesn't make much of a difference as explained). Only in Gordon the parts of the triangular decomposition are opposite, i.e. A<sup>l</sup> = K[V<sup>&ast;</sup>]<sub>W</sub>.


<a name="rrca-verma"></a>

### Working with modules

```C++
> W := TypeBReflectionGroup(2);
> Representations(~W);
> W`CharacterNames;
[ 11., 1.1, .11, 2., .2 ]

//Construct rational Cherednik algebra for W and generic GGOR parameter
> H:=RestrictedRationalCherednikAlgebra(W);

//We compute the baby Verma module for the W-representation the 2-dimensional
//representation 1.1:
> rho := W`Representations[0][2];
> M:=StandardModules(H, rho);
Graded module of dimension 16 over an algebra with generator degrees [ 0, 0, -1,
-1, 1, 1 ] over Multivariate rational function field of rank 2 over Rational
Field.

//I have implemented an own structure for garded modules that is used. 
//It's called ModGrOld (I started implementing a new type but this isn't fully 
//integrated right now.)
//Recall that as a vector space, M is isomorphic to K[V]_W \otimes \lambda.
//For each algebra generator of H (in this case w1, w2, y1, y2, x1, x2)
//the action is encoded by a matrix. The chosen basis for the coinvariant
//algebra can be viewed with
> W`CoinvariantAlgebra`Basis;
{@
    1,
    x2,
    x1,
    x2^2,
    x1*x2,
    x2^3,
    x1*x2^2,
    x1*x2^3
@}
//and the matrices of the generator actions can be viewed with
> M`Matrices;
[
    Sparse matrix with 16 rows and 16 columns over Multivariate rational
    function field of rank 2 over Rational Field,
    Sparse matrix with 16 rows and 16 columns over Multivariate rational
    function field of rank 2 over Rational Field,
    Sparse matrix with 16 rows and 16 columns over Multivariate rational
    function field of rank 2 over Rational Field,
    Sparse matrix with 16 rows and 16 columns over Multivariate rational
    function field of rank 2 over Rational Field,
    Sparse matrix with 16 rows and 16 columns over Multivariate rational
    function field of rank 2 over Rational Field,
    Sparse matrix with 16 rows and 16 columns over Multivariate rational
    function field of rank 2 over Rational Field
]
//So, the action of y1 is:
> Matrix(M`Matrices[3]);
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[-4*k1_1   4*k1_1   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   4*k1_1   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[2*k1_1 - 2*k2_1   -2*k1_1   0   0   0   0   0   0   0   0   0   0   0   0   0
    0]
[-4*k2_1   -2*k1_1 + 2*k2_1   0   0   0   0   0   0   0   0   0   0   0   0   0
    0]
[0   0   -4*k1_1   4*k1_1   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   4*k1_1   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   2*k1_1 - 2*k2_1   -2*k1_1   0   0   0   0   0   0   0   0   0   0   0
    0]
[0   0   -4*k2_1   -2*k1_1 + 2*k2_1   0   0   0   0   0   0   0   0   0   0   0
    0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   -2*k1_1 - 2*k2_1   2*k1_1   -4*k1_1   4*k1_1   0   0
    0   0   0   0]
[0   0   0   0   0   0   -4*k2_1   2*k1_1 + 2*k2_1   0   4*k1_1   0   0   0   0
    0   0]
[0   0   0   0   0   0   0   0   0   0   -2*k1_1 - 2*k2_1   2*k1_1   -4*k1_1
    4*k1_1   0   0]
[0   0   0   0   0   0   0   0   0   0   -4*k2_1   2*k1_1 + 2*k2_1   0   4*k1_1
    0   0]
//The degrees of the basis vectors of M are:
> M`RowDegrees;
[ 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4 ]

//Let's check if M is really a module for H (check all defining relations):
> IsModule(H,M);
true

//Let's compute a basis of the submodule of M spanned by the 16-th basis
//vector of M (which is x1*x2^3 \otimes e2), where e2 is the second basis
//vector of the W-representation rho):
> Spin(M, M.16);
[*
    (0   0   0   0   0   0   1   0   (2*k1_1^2 + 2*k1_1*k2_1)/(k1_1^2 + k2_1^2)
    -2*k1_1*k2_1/(k1_1^2 + k2_1^2)   0   0   0   0   0   0),
    (0   0   0   0   0   0   0   1   4*k1_1*k2_1/(k1_1^2 + k2_1^2)   (2*k1_1^2 -
    2*k1_1*k2_1)/(k1_1^2 + k2_1^2)   0   0   0   0   0   0),
    (0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0),
    (0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0),
    (0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0),
    (0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0),
    (0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0),
    (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1)
*]
//Hence, M.16 spans a non-trivial submodule.

//Let's try to compute the head of M. This will use my modular technique
//described in the CHAMP paper: specialize parameters, reduce to a finite field,
//use the MeatAxe, and lift everything back. This methods does not have to work,
//but it works surprisingly often. It's impossible to predict, however.K;
> res,L,J,P:=HeadOfLocalModule(M);
//The computation was successful. L is the head and J the radical of M.
//P describes the finite field specialization that was used.
//The function HeadOfLocalModule has many parameters to fine-tune the
//computation.
> L;
Graded module of dimension 8 over an algebra with generator degrees [ 0, 0, -1,
-1, 1, 1 ] over Multivariate rational function field of rank 2 over Rational
Field.
> L`Matrices[3];
[0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0]
[-4*k1_1   4*k1_1   0   0   0   0   0   0]
[0   4*k1_1   0   0   0   0   0   0]
[2*k1_1 - 2*k2_1   -2*k1_1   0   0   0   0   0   0]
[-4*k2_1   -2*k1_1 + 2*k2_1   0   0   0   0   0   0]
[0   0   2*k1_1 - 2*k2_1   -2*k1_1   0   0   0   0]
[0   0   -4*k2_1   -2*k1_1 + 2*k2_1   0   0   0   0]
> IsModule(H,L);
true

//Let's compute the Poincaré series and the graded W-character of L:
> PoincareSeries(L);
2 + 4*t + 2*t^2
// Every H-module is also a W-module. We can compute the corresponding (graded)
// decomposition as follows:
> InGroupSimples(H,M);
(      t t^2 + 1       t       t       t)
//Hence, L = t*(11.) + (t^2+1)*(1.1) + t*(.11) + t*(2.) + t*(.2)
//Note that in degree 0 of L there's a unique W-module, namely the 1.1=rho that
//we started with. This is a general fact and can be used to identify simple
//modules.
> IdentifyModule(H,L);
2   //the second irreducible W-module, i.e. 1.1=rho
```
<a name="rrca-mults"></a>

### Computing multiplicities

The standard module theory of the restricted rational Cherednik algebra leads to the following multiplicity problems:

* [P(λ) : Δ(μ)], functions ProjectivesInSimples and ProjectivesInSimplesQuantum
* [Δ(λ) : L(μ)], functions StandardsInSimples and StandardsInSimplesQuantum
* [L(λ) : μ], functions SimplesInGroupSimples and SimplesInGroupSimplesQuantum
* [Δ(λ) : μ], functions StandardsInGroupSimples and StandardsInGroupSimplesQuantum

In all cases, you can ask for both *graded* and *ungraded* multiplicities. I'm primarily targeting the *graded* multiplicities—from which you can of course immediately obtain the ungraded ones—and this what the above mentioned functions are doing. To represent the graded multiplicities, we can fix a system of representatives of the simples *up to grading shift* and then there are *two* ways to represent the graded multiplicities:

1. You collect for each simple of your system of representatives with which degree shift this occurs. We encode this information as a vector of size the number of simples in the system, and the entries are are (Laurent) polynomials in q. This is what the first named functions above are returning.
2. You fix a grading shift [n] and collect all simples of your system of representatives occuring with this grading shift. We encode this as a (Laurent) polynomial in q with coefficients a polynomial (actually just a linear expression) in the numner of simples in the system. This is what the "Quantum" functions are returning.  

So, these two ways of representing multiplicities is just about what to put first: simple module or grading shift. All this becomes clear in the examples below.

Before going to examples, I want to note that there are some relations between the multiplicities. Brauer reciprocity (combined with the [Δ(λ)] = [∇(λ)] result by Bellamy and myself) says that [P(λ) : Δ(μ)] = [Δ(λ) : L(μ)]. The multiplicities [Δ(λ) : μ] can be computed directly with a fake degree formula by Gordon. If you compute all the [Δ(λ) : μ] and manage to compute *all* the graded modules L(λ), then you know the [L(λ) : μ] and (by a result by Bellamy and myself) you can directly compute the [Δ(λ) : L(μ)] from the *matrix* formula ([Δ(λ) : μ]) = ([Δ(λ) : L(μ)])([L(λ) : μ]). I've built in many convenience functions that allow all these computations automatically. As I will explain allow, it's not always possible to get everything automatically because it's very complicated. For this reason, I'm not using matrices to store the multiplicities but associative arrays which are allowed to have undefined entries. If a method fails, you could still try to build the (irreducible) module in another way and if you succeed you can attach it to the corresponding array and keep computing. 

The ideal and simplest use case is illustrated in the following example:

```c++
> W := TypeBReflectionGroup(2);
> H:=RestrictedRationalCherednikAlgebra(W);
> StandardModules(~H); //computes all the standard modules
> H`StandardModules; //carries all the standard modules; numbering as in W`Representations[0]
Associative Array with index universe { 1 .. 5 }
> SimpleModules(~H); //(tries!) to compute all the simple modules by the method
//as described above
> H`SimpleModules;
Associative Array with index universe { 1 .. 5 }
> SimplesInGroupSimplesQuantum(~H);
> H`SimplesInGroupSimplesQuantum; 
Associative Array with index universe { 1 .. 5 }
// Let's look how these multiplicities are encoded
> f := H`SimplesInGroupSimplesQuantum[5]; f;
11.*q^4 + 1.1*q^3 + (.11 + 2.)*q^2 + 1.1*q + .2
//This means in L(5) we have the W-module 11. in with grading shift 4, the 
//W-module .11 + 2. with grading shift 2 etc.
> Parent(f);
Multivariate rational function field of rank 1 over Polynomial ring of rank 5 over Integer Ring
Variables: q
// Let's look at the other way to represent multiplicities.
> SimplesInGroupSimples(~H);
> H`SimplesInGroupSimples;
Associative Array with index universe { 1 .. 5 }
> H`SimplesInGroupSimples[5];
(    q^4 q^3 + q     q^2     q^2       1)
//The 2nd simple W-module occurs with multiplicity 1 in degrees 3 and 1 in L(5). 
//When we have all the information, we can also determine which standard module occurs 
//at the bottom of a projective (this gives the tilting permutation introduced by Bellamy
//and myself)
> StandardsAtBottomOfProjectives(~H);
> H`StandardsAtBottomOfProjectives;
Associative Array with index universe { 1 .. 5 }
> H`StandardsAtBottomOfProjectives[5];
<5,0>
//This means that Delta(5)[0] is at the bottom of P(5)
```

You can also can all of the above multiplicity functions with an additional integer argument (standing for a simple W-module λ in the fixed ordering) so that you just compute/get the information for the module corresponding to λ. 

I have implemented a function that produces MediaWiki code of all the representation-theoretic information.

```C++
MediaWiki(H);
```

#### Things that can go wrong

The multiplicity computations are extremely complicated. There are some things that can go wrong generically and that will require manual fiddling.

1. The function HeadOfLocalModule computes the unique irreducible quotient of a standard (which, just to remind you, is a huge module over a multivariate rational function field in characteristic zero). This uses a Las Vegas algorithm that I've presented in my original Champ paper. For some reason, it performs exceptionally well. But sometimes, you're just not lucky (like playing in Las Vegas). In this case, you could try to run the function manually a few more times or tweak its (complicated and unpredictable) parameters or you could try other things. 

2. The base field of the simple W-modules is not always the same as the base field of the group. I've simply taken the models from CHEVIE and took the minimal cyclotomic field containing all the entries of the matrices. Now, when you mix several representations—e.g. when you compute decompositon matrices with the above automatic methods—these varying base fields will cause problems (mathematically this is all trivial but the computer complains). So, *before* you do any kind of mixing computations, I advise changing all base rings to a common base ring. Here's an example:

   ~~~c++
   > W1:=ExceptionalComplexReflectionGroup(5); 
   //This model is defined over CyclotomicField(3). 
   //But there are representations defined over CyclotomicField(12). 
   //So, we'll change base rings to CyclotomicField(12) everywhere.
   > W := ChangeRing(W1, CyclotomicField(12));
   > W`DBDir := W1`DBDir; //Needed for loading reps (and everything else) from the database.
   > Representations(~W,0);
   > LiftRepresentationsToCommonBaseField(~W);
   ~~~

<a name="rrca-cmhyperplanes"></a>

### Calogero–Moser hyperplanes and families

The locus of parameters where the number of blocks of the restricted rational Cherednik algebra is less than the number of blocks for the generic algebra is known to be a union of hyperplanes. This locus is moreover known to be contained in the *Euler locus*, which is given by the pairwise differences of the values of the central characters of the simple modules for the generic algebra at the Euler element. To compute the Calogero–Moser hyperplanes, one can either compute the decomposition matrix on each Euler hyperplane and check whether the blocks are generic or not; or one can evaluate the central characters at the degree zero generators of the generic (non-restricted) rational Cherednik algebra. The database contains the Calogero–Moser hyperplanes for exceptional complex reflection groups whenever known:

```c++
> W:=TypeBReflectionGroup(2);
//We will compute the CM hyperplanes and families via central characters
> H := RationalCherednikAlgebra(W,0);
> CalogeroMoserFamilies(~H);
//The attribute H`CalogeroMoserFamilies is then an associative array with keys being the CM
//hyperplanes and entries being the CM families.
> H`CalogeroMoserFamilies;
Associative Array with index universe Polynomial ring of rank 2 over Rational Field
> Keys(H`CalogeroMoserFamilies);
{
k2_1,
k1_1,
1,
k1_1 + k2_1,
k1_1 - k2_1
}
//Here's a shortcut:
> CalogeroMoserHyperplanes(H);
{
k2_1,
k1_1,
1,
k1_1 + k2_1,
k1_1 - k2_1
}
//Get the base ring
> R := Universe(Keys( H`CalogeroMoserFamilies));
> R;
Polynomial ring of rank 2 over Rational Field
Order: Lexicographical
Variables: k1_1, k2_1
//Let's look at the CM families at k1_1-k2_1
> H`CalogeroMoserFamilies[R.1-R.2];
{
{ 1, 2, 5 },
{ 3 },
{ 4 }
}
```

For several complex reflection groups, the Calogero–Moser hyperplanes are stored in the database:

```
> W := ExceptionalComplexReflectionGroup(28);
> CalogeroMoserHyperplanes(W);
[
k2_1,
k1_1,
k1_1 - 2*k2_1,
k1_1 - k2_1,
k1_1 + k2_1,
k1_1 + 2*k2_1,
2*k1_1 - k2_1,
2*k1_1 + k2_1
]
```

In some exceptional cases, one can determine the CM families by using the *Euler families* (families one gets by evaluating the central characters only at the Euler element) combined with *supersingularity* and *rigidity* (see my papers "A Counter-Example to Martino’s Conjecture About Generic Calogero–Moser Families" and "Restricted Rational Cherednik Algebras"). Here's an example:

```c++
> W:=ExceptionalComplexReflectionGroup(23); //H3
> c:=CherednikParameter(W);
//The following gives the Euler families together with the value of the central character
{@
<{@ 9, 10 @}, 0>,
<{@ 7, 8 @}, 5*k1_1>,
<{@ 2 @}, 15*k1_1>,
<{@ 3 @}, -3*k1_1>,
<{@ 1 @}, -15*k1_1>,
<{@ 4 @}, 3*k1_1>,
<{@ 5, 6 @}, -5*k1_1>
@}
//Rigid representations
> RigidRepresentations(W,c);
{}
//Supersingular representations
> SupersingularRepresentations(W);
[ 5, 7, 9, 10 ]
//Try to determine the CM families by using Euler families and supersingularity.
//In this case we're lucky!
> CalogeroMoserFamiliesTry(W,c);
The Euler families are:
{@ 1 @}, {@ 2 @}, {@ 3 @}, {@ 4 @}, {@ 5, 6 @}, {@ 7, 8 @}, {@ 9, 10 @}

Singleton Euler families are CM families, so the following are already CM families:
{@ 1 @}, {@ 2 @}, {@ 3 @}, {@ 4 @}

The supersingular characters are:
5, 7, 9, 10

The following Euler families are CM families due to supersingularity:
{@ 5, 6 @}, {@ 7, 8 @}, {@ 9, 10 @}

Sucessfully determined the CM families. They are:
{@ 1 @}, {@ 2 @}, {@ 3 @}, {@ 4 @}, {@ 5, 6 @}, {@ 7, 8 @}, {@ 9, 10 @}
{@
{@ 1 @},
{@ 2 @},
{@ 3 @},
{@ 4 @},
{@ 5, 6 @},
{@ 7, 8 @},
{@ 9, 10 @}
@}
//For H4 we can determine many CM families but there's something left open:
> W:=ExceptionalComplexReflectionGroup(30);
> c:=CherednikParameter(W);
> c;
Mapping from: { 1 } to Polynomial ring of rank 1 over Cyclotomic Field of order 5 and degree 4
<1, 2*k1_1>
> CalogeroMoserFamiliesTry(W,c);
The Euler families are:
{@ 1 @}, {@ 2 @}, {@ 3, 5 @}, {@ 4, 6 @}, {@ 7, 8, 9, 10, 15, 16, 17, 22, 23, 24, 25, 26, 29, 30, 33, 34 @}, {@ 11, 13 @}, {@ 12, 14 @}, {@ 18, 20 @}, {@ 19, 21 @}, {@ 27 @}, {@ 28 @}, {@ 31 @}, {@ 32 @}

Singleton Euler families are CM families, so the following are already CM families:
{@ 1 @}, {@ 2 @}, {@ 27 @}, {@ 28 @}, {@ 31 @}, {@ 32 @}

The supersingular characters are:
3, 4, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 29, 30, 33, 34

The following Euler families are CM families due to supersingularity:
{@ 3, 5 @}, {@ 4, 6 @}, {@ 11, 13 @}, {@ 12, 14 @}, {@ 18, 20 @}, {@ 19, 21 @}
{@
{@ 1 @},
{@ 2 @},
{@ 27 @},
{@ 28 @},
{@ 31 @},
{@ 32 @},
{@ 3, 5 @},
{@ 4, 6 @},
{@ 11, 13 @},
{@ 12, 14 @},
{@ 18, 20 @},
{@ 19, 21 @}
@}
//Only the big Euler family {@ 7, 8, 9, 10, 15, 16, 17, 22, 23, 24, 25, 26, 29, 30, 33, 34 @} remains.
//I believe it's a CM family as well but I'm not sure.
```



<a name="rrca-db"></a>

### Database

In the database I have stored a lot of data (but certainly not all!) about the representation theory of restricted rational Cherednik algebras for exceptional complex reflection groups, especially multiplicities and Calogero–Moser hyperplanes and families. There could be even more data in the database and one could deduce more data from combining data but things were getting too complex and I didn't pursue this—you're free to expand this. Here's an example:

```c++
> W := ExceptionalComplexReflectionGroup(4);
> rec := RestrictedRationalCherednikAlgebraRepresentationTheory(W); 
> rec;
rec<recformat<ParameterRing, BlGen, DecGenStratification, Data> |
ParameterRing := Polynomial ring of rank 2 over Rational Field
Order: Lexicographical
Variables: k1_1, k1_2,
BlGen := [
k1_2,
k1_1,
k1_1 - 2*k1_2,
k1_1 - k1_2,
k1_1 + k1_2,
2*k1_1 - k1_2
],
Data := Associative Array with index universe Set of subsets of Polynomial ring of rank 2 over Rational Field>
//The entry BlGen gives the CM hyperplanes.
//The entry Data is an associative array indexed by (subsets of) the CM hyperplanes and giving
//information about the representation theory on the (intersection of) the hyperplanes.
> Keys(rec`Data);
{
{
k1_2
},
{
k1_1 - k1_2
},
{
k1_1 + k1_2
},
{
k1_1
},
{
2*k1_1 - k1_2
},
{
k1_1 - 2*k1_2
},
{
1
}
}
//Let's look at the data on the hyperplane 2*k1_1 - k2_1
> R := rec`ParameterRing;
> R;
Polynomial ring of rank 2 over Rational Field
Order: Lexicographical
Variables: k1_1, k2_1
> rec`Data[{2*R.1 - R.2}];
rec<recformat<SimpleDims, SimplePSeries, SimpleGModStruct, SimpleGradedGModStruct, VermaDecomposition, CMFamilies, CuspidalCMFamilies, VermaGradedDecomposition> |
SimpleDims := [ 24, 24, 3, 24, 24, 3, 18 ],
SimplePSeries := [
q^8 + 2*q^7 + 3*q^6 + 4*q^5 + 4*q^4 + 4*q^3 + 3*q^2 + 2*q + 1,
q^8 + 2*q^7 + 3*q^6 + 4*q^5 + 4*q^4 + 4*q^3 + 3*q^2 + 2*q + 1,
2*q + 1,
2*q^6 + 4*q^5 + 4*q^4 + 4*q^3 + 4*q^2 + 4*q + 2,
2*q^6 + 4*q^5 + 4*q^4 + 4*q^3 + 4*q^2 + 4*q + 2,
q + 2,
3*q^4 + 4*q^3 + 4*q^2 + 4*q + 3
],
SimpleGModStruct := [
(1 1 1 2 2 2 3),
(1 1 1 2 2 2 3),
(0 0 1 1 0 0 0),
(1 1 1 2 2 2 3),
(1 1 1 2 2 2 3),
(1 0 0 0 0 1 0),
(0 1 0 1 2 1 3)
],
SimpleGradedGModStruct := [
(1   q^8   q^4   q^7 + q^5   q^3 + q   q^5 + q^3   q^6 + q^4 + q^2),
(q^4   1   q^8   q^5 + q^3   q^7 + q^5   q^3 + q   q^6 + q^4 + q^2),
(0 0 1 q 0 0 0),
(          q^5             q           q^3       q^4 + 1     q^6 + q^2     q^4 + q^2 q^5 + q^3 + q),
(          q^3           q^5             q     q^4 + q^2       q^4 + 1     q^6 + q^2 q^5 + q^3 + q),
(q 0 0 0 0 1 0),
(            0           q^2             0             q       q^3 + q           q^3 q^4 + q^2 + 1)
],
VermaDecomposition := [
(1 0 0 0 0 0 0),
(0 1 0 0 0 0 0),
(0 0 1 0 0 1 1),
(0 0 0 2 0 0 0),
(0 0 0 0 2 0 0),
(0 0 2 0 0 2 2),
(0 0 3 0 0 3 3)
],
CMFamilies := {
{ 1 },
{ 3, 6, 7 },
{ 2 },
{ 4 },
{ 5 }
},
VermaGradedDecomposition := [
(1 0 0 0 0 0 0),
(0 1 0 0 0 0 0),
(  0   0   1   0   0 q^7 q^2),
(      0       0       0 q^2 + 1       0       0       0),
(      0       0       0       0 q^2 + 1       0       0),
(        0         0 q^7 + q^5         0         0   q^2 + 1   q^3 + q),
(0   0   q^6 + q^4 + q^2   0   0   q^5 + q^3 + q   q^4 + q^2 + 1)
]>

```

<a name="rrca-cellular"></a>
## Calogero–Moser cellular characters

```C++
//We compute the cellular characters for G4 at equal parameters
> W:=ExceptionalComplexReflectionGroup(4);
> c:=CherednikParameter(W,[1,1]);
//It's more efficient to compute the cellular characters per CM family or,
//more generally, for a union of CM families, like an Euler family.
//So, let's determine the Euler families.
> EulerFamilies(W,c);
{@
<{@ 7 @}, 0>,
<{@ 5, 6 @}, 2>,
<{@ 2, 3, 4 @}, -4>,
<{@ 1 @}, 8>
@}
//Now, let's compute the cellular characters for the Euler family {2,3,4}:
> cellchar := CalogeroMoserCellularCharacters(W,c,{@ 2,3,4 @});
> cellchar;
[1 1 2]
//This means there is one cellular character, and it decomposes as 
//1*chi_2 + 1*chi_3 + 2*chi_4, where chi_i is the character numbered by i
//in W`CharacterTable. In particular, the Euler family {2,3,4} is indeed 
//a CM family.
//We can also compute all cellular characters at once (but it is more efficient
//to do this per family):
> cellchar := CalogeroMoserCellularCharacters(W,c,{@ 1,2,3,4,5,6,7 @});
> cellchar;
[0 0 0 0 0 0 1]
[1 0 0 0 0 0 0]
[0 1 1 2 0 0 0]
[0 0 0 0 1 1 0]
//We can also automatically compute the cellular characters for all the Euler
//families. The result is a list of pairs consisting of an Euler family and the 
//decomposition matrix.
> cellchar := CalogeroMoserCellularCharacters(W,c);
[*
<
{@ 7 @},

[1]
>,

<
{@ 5, 6 @},

[1 1]
>,

<
{@ 2, 3, 4 @},

[1 1 2]
>,

<
{@ 1 @},

[1]
>
*]
//The function CalogeroMoserCellularCharacters will automatically determine a
//regular vector. It will be more efficient to manually provide an explicit 
//regular vector since otherwise Gaudin operators will first be computed 
//generically and not already specialized in a regular vector:
> V:=VectorSpace(W);
> vreg:=V![-3,5]; //a random choice
> IsRegular(W,vreg);
true
//Let's use this particular regular vector for the computation of cellular 
//characters.
> cellchar := CalogeroMoserCellularCharacters(W,c,{@ 2,3,4 @} : vreg:=vreg);

//The computation of cellular characters uses the Gaudin operators. 
//Here are the functions for Gaudin operators:
> GaudinOperator(W,c); //the full Gaudin operator
> GaudinOperator(W,c, W`Representations[0][7]); //Gaudin operator for representation #7
> GaudinOperator(W,c, V![1,1] ); //Gaudin operator specialized at y=(1,1)
> GaudinOperator(W,c, V![1,1], W`Representations[0][7] ); //Gaudin operator specialized at y=(1,1) and for representation #7
```

