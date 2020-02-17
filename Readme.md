# CHAMP

A Cherednik Algebra Magma Package. By [Ulrich Thiel](https://ulthiel.com/math), 2013–2020.

## Scope

With this package you can:
* compute in rational Cherednik algebras (see [Etingof-Ginzburg](https://arxiv.org/abs/math/0011114))
* compute generators and a presentation of the center of the rational Cherednik algebra at t=0 (the coordinate algebra of the Calogero-Moser space)
* compute Poisson brackets on the center (towards symplectic leaves)
* compute decomposition matrices of baby Verma modules and graded characters of simples for restricted rational Cherednik algebras (see [Gordon](https://arxiv.org/abs/math/0202301)).

The parameters can always be arbitrary, including generic parameters valued in polynomial rings or rational function fields. This document contains a complete overview of the functionality with many examples. The theory and algorithms are discussed in the following publications:
* U. Thiel, CHAMP: A Cherednik Algebra Magma Package
LMS J. Comput. Math. 18 (2015), no. 1, 266–307.
* C. Bonnafé and U. Thiel, Calogero–Moser families and cellular characters: computational aspects (with C. Bonnafé). In preparation (2020).

### Contents

[Downloading an running](#downloading)  
[Complex reflection groups](#reflgroups)  
[Rational Cherednik algebras](#che)  
&nbsp;&nbsp;&nbsp;&nbsp;[Parameters](#params)  
&nbsp;&nbsp;&nbsp;&nbsp;[Rational Cherednik algebras at t=0 and Calogero-Moser spaces](#cmspaces)  
[Restricted rational Cherednik algebras](#rrca)  

<a name="downloading"></a>

## Downloading and running

You need a [Magma](http://magma.maths.usyd.edu.au/magma/) version of at least 2.19 (current version is 2.25). You can then download the [latest CHAMP release](https://github.com/ulthiel/champ/releases/latest) and start it by running ```./champ```. **Important:** for full functionality of CHAMP, you have to download the ReflectionGroups database from the release assets as well and extract it in the ```DB``` directory of CHAMP.

Alternatively, you can clone the git repository. **Important**: due to large binary files in the database, it is stored with [Git Large File Storage](https://git-lfs.github.com). You first have to install this extension as described in the link. Then you can do a ```git clone https://ulthiel.github.com/champ/``` as usual and this will also clone the database.

I advise to once run ```./selfcheck``` in the directory ```SelfCheck```. (The ReflectionGroups selfcheck will take a bit of time but if the first few are fine, the rest should be fine as well).

<a name="reflgroups"></a>

## Complex reflection groups

Models for several complex reflection groups, their character tables, character names, models for irreducible representations, etc. is stored in the ReflectionGroups database. The data is taken from (and compatible with) J. Michel's [CHEVIE](https://webusers.imj-prg.fr/~jean.michel/chevie/chevie.html) package. The reason for using a database is that we need consistent labelings (of e.g. characters) that allow us to compare results with the literature. A general philosophy in CHAMP is that most objects (like groups) will have attributes (like CharacterTable) which are set by a similarly named procedure operating on the object (using the ~ operator). Usually, it is first checked whether the data exists in the database; if not, it will be computed in a consistent way.

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
//and/or natural choices) can be created with the following functions:
//ExceptionalComplexReflectionGroup, SymmetricReflectionGroup,
//TypeBReflectionGroup, TypeDReflectionGroup, DihedralReflectionGroup,
//CyclicReflectionGroup, ImprimitiveReflectionGroup.
//
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

<a name="che"></a>

## Rational Cherednik algebras

```C++
//Create the rational Cherednik algebra for t and c generic (valued in a
//polynomial ring)
> W := TypeBReflectionGroup(2); //Weyl group of type B2 as above
> H := RationalCherednikAlgebra(W : Type:="EG"); //I will explain the EG below
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
//confusing, but in the end it's less confusing than trying to make artifically
//make Magma act on the left.
```

<a name="params"></a>

### Parameters

This topic is a bit technical but important. There are two kinds of parameters involved in the relations for the rational Cherednik algebra: a *t-parameter* and a *c-parameter*. Let's take a commutative ring R as base ring. The t-parameter is some fixed element of R; the c-parameter is a function c:Refl(W)/W -> R from the conjugacy classes of reflections of W to R. For example, we can let R be a polynomial ring K[t,c_1,...,c_r] and define the parameters t and c in the obvious way. In this case we say the parameters are *generic*. If I is an ideal of R, we can also consider R/I as new base ring and get parameters with are *generic for the subscheme* defined by I. For example, we could take a polynomial ring R=K[t,c] and set c(s)=c for all c. This would be the generic *equal* parameter case.

For the construction of the rational Cherednik algebra in CHAMP you can take as base ring R any K-algebra that can be defined in Magma, where K is the base field of the reflection group W, and as parameters you can take any t and maps c with values in R. In particular, you can work with generic parameters, generic parameters on a, say, hyperplane, or special parameters taking values in your base field K. You have complete freedom.

[Ginzburg-Guay-Opdam-Rouquier](https://arxiv.org/abs/math/0212036) considered a Fourier transform on the c-parameter space which makes some expressions in the parameters much simpler (such as equations for the Calogero-Moser hyperplanes). I will refer to these as *k-parameters*. While the c-parameters by Etingof-Ginzburg are indexed by conjugacy classes of reflections, the k-parameters have a double index: the first indexes an orbit [H] of reflection hyperplanes, the second is an index between 0 and |W_H|-1, where W_H is the stabilizer of a representative of [H]. Of course, in the end the number of parameters is the same. By default, CHAMP uses k-parameters.

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
> cH:=SpecializeCherednikParameterInHyperplane(c, R.1-R.2);
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
> k:=CherednikParameter(W);
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

### Rational Cherednik algebras at t=0 and Calogero-Moser spaces

The rational Cherednik algebra H_{t=0,c} has a big center Z_c: it is a Poisson deformation of the symplectic singularity (V \oplus V^*)/W, where W acts on V. The associated variety is called the *Calogero-Moser space* of W at parameter c. CHAMP can compute algebra generators of Z_c and also a presentation of this algebra (the former works even for large groups like F4, the latter involves rather complicated invariant theory computations which are even for small dihedral groups too much; but you can still get some ideas).

The ReflectionGroups database contains generators of Z_0 (undeformed case) and Z_c (c generic) for several cases. Some of the elements are extremely large (for G11 there is one taking up >100MB compressed and >500MB uncompressed)! By default, all functions check the database first and load the data from there if available.


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
[ <2, 0>, <1, 1>, <0, 2>, <4, 0>, <3, 1>, <2, 2>, <1, 3>, <0, 4> ]
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
//Poisson brackets of elements of Z_c can be computed as well
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

```Delphi
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
> H.5*H.2;
[ 1  0]
[ 1 -1]*(x1 + x2)

//You can convert H into a matrix algebra
> A:=MatrixAlgebra(H);
> A;
Matrix Algebra of degree 512 with 6 generators over Multivariate rational
function field of rank 2 over Rational Field

//We compute the Jacobson radical for the equal parameter case k=[1,1]:
> c := CherednikParameter(W,[1,1]);
> H := RestrictedRationalCherednikAlgebra(W,c);
> A := MatrixAlgebra(H);
> time J:=JacobsonRadical(A); J;
Time: 182.370
Matrix Algebra [ideal of A] of degree 512 and dimension 346 over Rational Field
```
