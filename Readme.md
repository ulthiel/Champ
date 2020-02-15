# CHAMP

A Cherednik Algebra Magma Package. By [Ulrich Thiel](https://ulthiel.com/math), 2013–2020.

## Scope

With this package you can:
* compute in rational Cherednik algebras (see [Etingof-Ginzburg](https://arxiv.org/abs/math/0011114))
* compute generators and a presentation of the center of the rational Cherednik algebra at t=0 (the coordinate algebra of the Calogero-Moser space)
* compute Poisson brackets on the center (towards symplectic leaves)
* compute decomposition matrices and graded characters for restricted rational Cherednik algebras (see [Gordon](https://arxiv.org/abs/math/0202301)).

The parameters can always be arbitrary, including generic parameters valued in polynomial rings or rational function fields. In the following we will give a complete overview of the functionality. The theory and algorithms is discussed in the following publications:
* U. Thiel, CHAMP: A Cherednik Algebra Magma Package
LMS J. Comput. Math. 18 (2015), no. 1, 266–307.
* C. Bonnafé and U. Thiel, Calogero–Moser families and cellular characters: computational aspects (with C. Bonnafé). In preparation (2020).

## Downloading and running

You need a [Magma](http://magma.maths.usyd.edu.au/magma/) version of at least 2.19 (current version is 2.25). It's most convenient to simply download the [latest release](https://github.com/ulthiel/champ/releases/latest). You can start champ by running ```./champ```.

**Important.** For full functionality of CHAMP, you have to download the ReflectionGroups database from the release assets as well and extract it in the ```DB``` directory of CHAMP.

Alternatively, you can clone the git repository. This has a minor complication: due to large binary files in the database, it is stored with [Git Large File Storage](https://github.com/ulthiel/champ/releases/latest). You first have to install this extension as described in the link. Then you can do a ```git clone https://ulthiel.github.com/champ/``` as usual and this will also clone the database.

## Reflection groups

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
```

Other types of reflection groups (with connection to data from the database and/or natural choices) can be created with the following functions: ExceptionalComplexReflectionGroup, SymmetricReflectionGroup, TypeBReflectionGroup, TypeDReflectionGroup, DihedralReflectionGroup, CyclicReflectionGroup, ImprimitiveReflectionGroup.

You can also load some special models directly from the database as in the following example:

```C++
//Particular model of B2 used by Bonnafé-Rouquier in some computation
> W := CHAMP_GetFromDB("ReflectionGroups/B2_BR", "GrpMat");
> W;
MatrixGroup(2, Rational Field)
Generators:
    [0 1]
    [1 0]

    [-1  0]
    [ 0  1]
```

## Rational Cherednik algebras

```C++
//Create the rational Cherednik algebra for t and c generic (valued in a
//polynomial ring)
> W := TypeBReflectionGroup(2); //Weyl group of type B2 as above
> H := RationalCherednikAlgebra(W);
Rational Cherednik algebra
Generators:
    w1, w2, y1, y2, x1, x2
Generator degrees:
    0, 0, -1, -1, 1, 1
Base ring:
    Polynomial ring of rank 3 over Rational Field
    Order: Lexicographical
    Variables: t, k1_1, k2_1
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
    <1, 2*k1_1>
    <2, 2*k2_1>

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
[-1  1]*(2*k2_1)
+
[-1  2]
[ 0  1]*(2*k1_1)
+
[1 0]
[0 1]*(y1*x1 + t)

//IMPORTANT: In Magma, matrices are acting from the *right* on vectors. Hence,
//to keep everything consistent, I have implemented the *opposite* of the
//rational Cherednik alebra as usually written on paper. This may be a bit
//confusing, but in the end it's less confusing than trying to make Magma act on
//the left.
```
