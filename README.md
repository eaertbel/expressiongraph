/*! \mainpage Expression graph library

# Expression graph library.


This library provides expression trees for representation of geometric expressions and automatic differentiation
of these expressions.  This enables to write down geometric expressions at the position level, and automatically 
compute Jacobians and higher order derivatives efficiently and without loss of precision.
The library is built upon the KDL library and uses the same geometric primitives.  Most operators and functions for
the basic primitives of KDL, have an equivalent in the expression tree library.


author: 

               
```
#!c++

                Erwin Aertbelien
                KU Leuven - Dep. of Mechanical Engineering
                Celestijnenlaan 300B - bus 2480
                3001 HEVERLEE
                Belgium

                E-mail: <first name> dot <last name> at kuleuven.be

```

### EUPL License ###

Copyright 2014 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering

Licensed under the EUPL, Version 1.1 only (the "Licence"); You may not use this work except in compliance with the Licence. You may obtain a copy of the Licence at:

[EUPL license](http://ec.europa.eu/idabc/eupl)

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or im plied.
 See the Licence for the specific language governing permissions and limitations under the Licence.


You can submit patches to this repository to the above original author.
To allow for dual licensing, we will only apply patches for which the authors
have given a release.

This library is first release to the public at November,25th, 2014,

### Libraries used: ###
This library links dynamically  with the KDL library (http://wiki.ros.org/kdl). The KDL library is distributed under LGPL.

### Contributing authors: ###
*