# Structure manipulation

An attempt to write some tools to create molecular structures automatically.
Structures are made of MMT, modifier and polyamide-6.


## construct_composite_10

Constructs compisite that is very similar to old segregated system:
    720 mmt atoms (periodic platelet with 12 substitutions Al -> Mg)
    840 modifier atoms (12 molecules x 70 atoms)
    1920 polymer atoms (10 molecules x 192 atoms)


## construct_composite_100

A structure with long polymer chains, designed to check wether the relaxation 
process depeds on the polymer lenght. Also, as soon as the crystallization 
behavior shall be studied, it is important to check the crystallization of 
longer chains.
    720*9 mmt atoms
    840*9 modifiler atoms
    1902*75 polymer atoms


## Source code:

### datafile_content

A container for storing datafile content (also read and write)


### datafile_content_extractor

A tool to extract some atoms from DatafileContent (but it corrupts the original
DatafileContent! Maybe, it would be changed later).


### datafile_content_modifiers

Some utilities (for getting bounding box, etc.). Some of them may be incorrect.
Maybe, it would be changed later.


### insert_poly_into_mmt

Insert polymer chain (stored as DatafileContent) into the empty space of the 
structure made of MMT. Works not in the best way, since the emptiness is 
defined not in general way. Is is good to improve it, because it causes some 
troubles (avoidable by the specific consecuence of structure preparation steps).


### utils

A vary simple utils (for now, inly translate some atoms without any cheks).



## lmp inputs

### lmp_in

A simple script to relax structures, made for test purpose mainly.


### lmp_replicate

Replicate structure (and also delete some atoms maybe). Maybe, it is good idea
to implement this functionality in my scripts.
