# CHAMP

A Cherednik Algebra Magma Package. By [Ulrich Thiel](https://ulthiel.com/math), 2013–2020.

With this package you can:
* compute in rational Cherednik algebras (see [Etingof-Ginzburg](https://arxiv.org/abs/math/0011114))
* compute generators of the center of the rational Cherednik algebra at t=0 and even a presentation (i.e. the coordinate algebra of the Calogero-Moser space)
* compute decomposition matrices and graded characters for restricted rational Cherednik algebras (see [Gordon](https://arxiv.org/abs/math/0202301)).

The parameters can be arbitrary, including generic parameters in polynomial rings or rational function fields.

The theory and algorithms is discussed in the following publications:
* U. Thiel, CHAMP: A Cherednik Algebra Magma Package
LMS J. Comput. Math. 18 (2015), no. 1, 266–307.
* C. Bonnafé and U. Thiel, Calogero–Moser families and cellular characters: computational aspects (with C. Bonnafé). In preparation (2020).


## Downloading and running

You need a [Magma](http://magma.maths.usyd.edu.au/magma/) version of at least 2.19 (current version is 2.25). It's most convenient to simply downloaded the [latest release](https://github.com/ulthiel/champ/releases/latest). You can start champ by running ```./champ```.

**For full functionality of CHAMP, you have to download the ReflectionGroups database from the assets as well and extract it in the ```DB``` directory of CHAMP.**

Alternatively, you can clone the git repository; but this has a little twist: the database is stored with [Git Large File Storage](https://github.com/ulthiel/champ/releases/latest), and you first have to install this extension. Then you can do a ```git clone https://ulthiel.github.com/champ/``` as usual.

## Reflection groups

Models for several complex reflection groups, their character tables, character names, models for irreducible representations, and further data is stored in the ReflectionGroups database. The data is taken from J. Michel's [CHEVIE](https://webusers.imj-prg.fr/~jean.michel/chevie/chevie.html) package.

The following examples show how to use all functions around complex reflection groups:

```
//Load the Weyl group B2 in a reflection representation
> W := TypeBReflectionGroup(2);

//Character tables and standard character names are stored in the database.
> CharacterTable(~W);
> W`CharacterNames;
[ 11., 1.1, .11, 2., .2 ]

//The above illustrates a general philosophy: most objects (like groups)
//will have attributes (like CharacterTable) which are set by a procedure
//operating on the object (hence the ~ opeator). A call to Magma's internal
//CharacterTable(W) would compute a character table and we don't have a
//labeling anymore. Keep this in mind.

//Load models for the irreducible representations. Their numbering will match
//the one from the database.
> Representations(~W);
> W`Representations[0]; //I wanted to use positive characteristic
                        //representations one day, hence the 0.

//Get the fake degrees (graded W-character of the coinvariant algebra)
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

Other types of reflection groups (with connection to data from the database and/or natural choices) can be created with the following functions: ExceptionalComplexReflectionGroup, SymmetricReflectionGroup, CyclicReflectionGroup, TypeBReflectionGroup, TypeDReflectionGroup, ImprimitiveReflectionGroup, DihedralReflectionGroup.

## Rational Cherednik algebras

```
//Create the rational Cherednik algebra for t and c generic (valued in a
//polynomial ring)
> W := TypeBReflectionGroup(2);
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

```
