# CHAMP

A Cherednik Algebra Magma Package. By [Ulrich Thiel](https://ulthiel.com/math), 2013–2021.

## Scope

With this package you can:
* compute in rational Cherednik algebras (as introduced by [Etingof–Ginzburg](https://arxiv.org/abs/math/0011114)),
* compute generators and a presentation of the center of the rational Cherednik algebra at t=0 (the coordinate algebra of the Calogero-Moser space),
* compute Poisson brackets on the center,
* compute decomposition matrices of baby Verma modules and graded characters of simples for restricted rational Cherednik algebras (as introduced by [Gordon](https://arxiv.org/abs/math/0202301)).

The underlying reflection groups can be arbitrary and also the parameters can be arbitrary, including generic parameters valued in polynomial rings or rational function fields. This document contains a complete overview of the functionality with many examples. The theory and algorithms are discussed in the following publications:
* U. Thiel, [CHAMP: A Cherednik Algebra Magma Package](https://arxiv.org/abs/1403.6686), LMS J. Comput. Math. 18 (2015), no. 1, 266–307.
* C. Bonnafé and U. Thiel, Calogero–Moser families and cellular characters: computational aspects (with C. Bonnafé). In preparation (2021).

## Contents

[1. Downloading an running](#downloading)  
[2. Complex reflection groups](#reflgroups)  
[3. Rational Cherednik algebras](#che)  
&nbsp;&nbsp;&nbsp;&nbsp;[3.1 Parameters](#params)  
&nbsp;&nbsp;&nbsp;&nbsp;[3.2 Rational Cherednik algebras at t=0 and Calogero-Moser spaces](#cmspaces)  
&nbsp;&nbsp;&nbsp;&nbsp;[3.3 Poisson brackets](#poisson-brackets)  
[4. Restricted rational Cherednik algebras](#rrca)  
&nbsp;&nbsp;&nbsp;&nbsp;[4.1 Representation theory](#rrca-rep)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.1.1 Conventions](#rrca-conv)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.1.2 Working with modules](#rrca-verma)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[4.1.3 Computing multiplicities](#rrca-mults)

<a name="downloading"></a>

## Downloading and running

You need a [Magma](http://magma.maths.usyd.edu.au/magma/) version of at least 2.19 (current version is 2.25). You can then download the [latest CHAMP release](https://github.com/ulthiel/champ/releases/latest) and start it by running ```./champ```. **Important:** for full functionality of CHAMP, you have to download the ReflectionGroups database from the [release assets](https://github.com/ulthiel/champ/releases/latest) as well and extract it in the ```DB``` directory of CHAMP.

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

<a name="che"></a>

## Rational Cherednik algebras

The definitition of rational Cherednik algebras used in Champ is exactly the one by [Etingof–Ginzburg](https://arxiv.org/abs/math/0011114). It's best to begin with an example straightaway.

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

This topic is a bit technical but important. There are two kinds of parameters involved in the relations for the rational Cherednik algebra: a *t-parameter* and a *c-parameter*. Let's take a commutative ring R as base ring. The t-parameter is some fixed element of R; the c-parameter is a function c:Refl(W)/W → R from the conjugacy classes of reflections of W to R. For example, we can let R be a polynomial ring K[t,c<sub>1</sub>,...,c<sub>r</sub>] and define the parameters t and c in the obvious way. In this case we say the parameters are *generic*. If I is an ideal of R, we can also consider R/I as new base ring and get parameters with are *generic for the subscheme* defined by I. For example, we could take a polynomial ring R=K[t,c] and set c(s)=c for all c. This would be the generic *equal* parameter case.

For the construction of the rational Cherednik algebra in CHAMP you can take as base ring R any K-algebra that can be defined in Magma, where K is the base field of the reflection group W, and as parameters you can take any t and maps c with values in R. In particular, you can work with generic parameters, generic parameters on a, say, hyperplane, or special parameters taking values in your base field K. You have complete freedom.

[Ginzburg-Guay-Opdam-Rouquier](https://arxiv.org/abs/math/0212036) considered a Fourier transform on the c-parameter space which makes some expressions in the parameters much simpler (such as equations for the Calogero-Moser hyperplanes). I will refer to these as *k-parameters*. While the c-parameters by Etingof-Ginzburg are indexed by conjugacy classes of reflections, the k-parameters have a double index: the first indexes an orbit [H] of reflection hyperplanes, the second is an index between 0 and |W<sub>H</sub>|-1, where W<sub>H</sub> is the stabilizer of a representative of [H]. Of course, in the end the number of parameters is the same. By default, CHAMP uses k-parameters.

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

The rational Cherednik algebra H<sub>t=0,c</sub> has a big center Z<sub>c</sub>: it is a Poisson deformation of the symplectic singularity (V ⊕ V<sup>&ast;</sup>)/W, where W acts on V. The associated variety is called the *Calogero-Moser space* of W at parameter c. CHAMP can compute algebra generators of Z<sub>c</sub> and also a presentation of this algebra (the former works even for large groups like F<sub>4</sub>, the latter involves rather complicated invariant theory computations which are even for small dihedral groups too much; but you can still get some ideas).

The ReflectionGroups database contains generators of Z<sub>0</sub> (undeformed case) and Z<sub>k</sub> (k generic) for several cases. Some of the elements are extremely large (for G<sub>11</sub> there is one taking up >100MB compressed and >500MB uncompressed)! By default, all functions check the database first and load the data from there if available.


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
```

<a name="poisson-brackets"></a> 
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

### Representation theory

In CHAMP you can compute baby Verma modules (also called standard modules) for restricted rational Cherednik algebras (as defined by [Gordon](https://arxiv.org/abs/math/0202301)). Using modular lifting techniques I introduced in [my paper](https://arxiv.org/abs/1403.6686) you can compute the heads of standard modules (which then give all the simples of the restricted rational Cherednik algebra) as graded modules (also giving the graded W-character) and the (graded) decomposition matrix of standard modules into simples. It works surprisingly well even in huge and complicated examples, and for generic parameters as well.

<a name="rrca-conv"></a>

#### Conventions

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

#### Working with modules

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

In all cases, you can ask for both *graded* and *ungraded* multiplicities. I'm primarily targeting the *graded* multiplicities—from which you can of course immediately obtainen the ungraded ones—and this what the above mentioned functions are doing. To represent the graded multiplicities, we can fix a system of representatives of the simples *up to grading shift* and then there are *two* ways to represent the graded multiplicities:

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

1. The multiplicity computations are extremely complicated. They use the function HeadOfLocalModule to compute the unique irreducible quotient of a standard (which, just to remind you, is a huge module over a multivariate rational function field in characteristic zero). This uses a Las Vegas algorithm that I've presented in my original Champ paper. For some reason, it performs exceptionally well. But sometimes, you're just not lucky (like playing in Las Vegas). In this case, you could try to run the function manually a few more times or tweak its (complicated and unpredictable) parameters or you could try other things. 

2. The base field of the simple W-modules is not always the same as the base field of the group. I've simply taken the models from CHEVIE and took the minimal cyclotomic field containing all the entries of the matrices. Now, when you mix several representations—e.g. when you compute decompositon matrices with the above automatic methods—these varying base fields will cause problems (mathematically this is all trivial but the computer complains). So, *before* you do any kind of mixing computations, I advise doing

   ~~~c++
   > LiftRepresentationsToCommonBaseField(~W);
   ~~~

   This changes the base fields of all the simple W-modules to *one* common base field.