## ReflectionGroups

Contains particular models of (complex) reflection groups, mostly imported from CHEVIE. Recall that these groups are labeled G(m,p,n) and G<sub>*n*</sub> by Shephard and Todd. For reflection groups which are Coxeter groups we prefer the Coxeter name. Also, for cyclic groups we use "Cyc".

## From CHEVIE

| Name          | Description                                                  |
| ------------- | ------------------------------------------------------------ |
| A*n*_CHEVIE   | Irreducible complex reflection group of type A<sub>*n*</sub> (i.e. symmetric group S<sub>n+1</sub>) from CHEVIE (*n*=1..25). Is loaded by ```ComplexReflectionGroup(1,1,n+1)```,  ```TypeAReflectionGroup(n)``` , ```SymmetricReflectionGroup(n+1)```. |
| B*n*_CHEVIE   | Irreducible complex reflection group of type B<sub>*n*</sub> from CHEVIE (*n*=2..25). Is loaded by ```ComplexReflectionGroup(2,1,n)```, ```TypeBReflectionGroup(n)```. |
| Cyc*m*_CHEVIE | Irreducible reflection representation of cyclic group of order *m* from CHEVIE (*m*=2..25). Is loaded by ```ComplexReflectionGroup(m,1,1)```, ```CyclicReflectionGroup(m)```. |
| D*n*_CHEVIE   | Irreducible complex reflection group of type D<sub>*n*</sub> from CHEVIE (*n*=2..25). Is loaded by ```ComplexReflectionGroup(2,2,n)```, ```TypeDReflectionGroup(n)```. |
| Dih*m*_CHEVIE | Irreducible reflection representation of dihedral group of order 2*m* from CHEVIE (*m*=5..25). Is loaded by ```ComplexReflectionGroup(m,m,2)```, ```DihedralReflectionGroup(m)```. |
| E*n*_CHEVIE   | Exceptional irreducible complex reflection group of type E<sub>*n*</sub> from CHEVIE (*n*=6,7,8). Is loaded by ```ComplexReflectionGroup(35)```, ```ComplexReflectionGroup(36)```, and ```ComplexReflectionGroup(37)```, respectively. |
| F4_CHEVIE     | Exceptional irreducible complex reflection group of type F<sub>4</sub> from CHEVIE (*n*=3,4). Is loaded by ```ComplexReflectionGroup(28)```. |
| G*n*_CHEVIE   | Exceptional irreducible complex reflection group of type G<sub>*n*</sub> from CHEVIE (*n*=4..37, except *n*=23,28,30,35,36,37 because they carry Coxeter names). Is loaded by ```ComplexReflectionGroup(n)```. |
| H*n*_CHEVIE   | Exceptional irreducible complex reflection group of type H<sub>*n*</sub> from CHEVIE (*n*=3,4). Is loaded by ```ComplexReflectionGroup(23)``` and ```ComplexReflectionGroup(30)```. |

## Other models

| Name         | Description                                                  |
| ------------ | ------------------------------------------------------------ |
| B2_BR        | Model of B<sub>2</sub> used by Bonnafé–Rouquier.             |
| F4_Bonnafe   | Model of F<sub>4</sub> by Bonnafé from Apr 14, 2020.         |
| G24_Bonnafe  | Model of G<sub>24</sub> by Bonnafé from Apr 14, 2020.        |
| H*n*_Bonnafe | Model of H<sub>n</sub> by Bonnafé from Apr 14, 2020 (*n*=3,4). |