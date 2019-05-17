# Structure manipulation


https://github.com/ansko/md_work

An attempt to write some tools to create molecular structures automatically.
Structures are made of MMT, modifier and polyamide-6.

It seems that 75 chains each 100-monomers long in mmt replicated 3 times in both
x and y directions is a good configuration. It means that there are 7500 monomers
in the cell. Summary of systems:

Length      --- Chains number --- Atoms Count

in monomers --- per mmt_331   ---

 10         --- 750           --- 1560*9 + 192*750       = 158040 +

 25         --- 300           --- 1560*9 + (19*25+2)*300 = 157140 +

 50         --- 150           --- 1560*9 + (19*50+2)*150 = 156840 +

100         ---  75           --- 1560*9 + (19*100+2)*75 = 156690 +


## construct_composite_*

Construct compisites that are similar to each other:

    1560*9 = mmt and mod

    N * M-monomer-long chains (N*M = 7500 = 75*100 = 50*150 = 25*300 = 10*750)

Results are in data_structures/comp_mon*_n*.data



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
