# Genome Tree
Given a list of genomes, will generate a phylogenic tree.

## License
This software is dual-licensed under GPL Version 2 and/or GPL Version 3. You may
use this software according to these licenses as is most appropriate for your
project on a case-by-case basis. Both licenses can be found in the root
directory of the repository.

## Prerequisites
Ensure Python and Cargo (Rust) are both installed.

## Basic Instructions
Download genomes from GenBank, which should be received as a .zip type. Some
genomes may be provided with this repository. Place the .zip files in the
genomes_raw directory.

To automatically unzip and preprocess the zip files, run:

    python preprocess.py

A new directory called genomes should have been created. We are ready to begin.

To generate the tree, simply run:

    cargo run --release
    
A phylogenetic tree should have been exported as an xml to the root directory
in a file called 'phylo_tree.xml'. Thank you for using this software.

## Things to note
This software was developed and tested solely on a Linux machine. Python and
Rust are both cross-platform, and as such this should work on other systems
as well, however I can not guarantee this simply because I have not tested it.
