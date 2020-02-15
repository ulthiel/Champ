# CHAMP

A Cherednik Algebra Magma Package. By [Ulrich Thiel](https://ulthiel.com/math), 2013–2020.

## Scope

With this package you can:
* compute in rational Cherednik algebras (see [Etingof-Ginzburg](https://arxiv.org/abs/math/0011114))
* compute generators and a presentation of the center of the rational Cherednik algebra at t=0 (the coordinate algebra of the Calogero-Moser space)
* compute Poisson brackets on the center (towards symplectic leaves)
* compute decomposition matrices and graded characters for restricted rational Cherednik algebras (see [Gordon](https://arxiv.org/abs/math/0202301)).

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

<a name="downloading"></a>

## Downloading and running

You need a [Magma](http://magma.maths.usyd.edu.au/magma/) version of at least 2.19 (current version is 2.25). You can then download the [latest CHAMP release](https://github.com/ulthiel/champ/releases/latest) and start it by running ```./champ```. **Important:** for full functionality of CHAMP, you have to download the ReflectionGroups database from the release assets as well and extract it in the ```DB``` directory of CHAMP.

Alternatively, you can clone the git repository. This has a minor complication: due to large binary files in the database, it is stored with [Git Large File Storage](https://github.com/ulthiel/champ/releases/latest). You first have to install this extension as described in the link. Then you can do a ```git clone https://ulthiel.github.com/champ/``` as usual and this will also clone the database.

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
> W:=TypeBReflectionGroup(2);
> CherednikParameter(W);
Mapping from: { 1 .. 2 } to Polynomial ring of rank 2 over Rational Field
    <1, 2*k1_1>
    <2, 2*k2_1>

//The labeling of orbits of reflection hyperplanes is consistent with what is
//stored in
> W`ReflectionLibrary;
//This is an array indexed by orbits of reflection hyperplanes. Each entry is
//again an array indexed by reflection hyperplanes in this orbit. The entries
//of this array are the reflections for the corresponding hyperplane.
```

<a name="cmspaces"></a>

### Rational Cherednik algebras at t=0 and Calogero-Moser spaces
