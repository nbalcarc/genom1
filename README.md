# Genome Tree
Given a list of genomes, will generate a phylogenic tree.

## License
This software is dual-licensed under GPL Version 2 and/or GPL Version 3. You may
use this software according to these licenses as is most appropriate for your
project on a case-by-case basis. Both licenses can be found in the root
directory of the repository.

## Basic Instructions
Download genomes from GenBank, which should be received as a .zip type. Some
genomes may be provided with this repository. Place the .zip files in the
genomes_raw directory.

Run:
    'python preprocess.py'
A new directory called genomes should have been created. We are ready to begin.

To generate the tree, simply run:
    'cargo run'
A phylogenetic tree should have been exported as an xml to the root directory
in a file called 'phylo_tree.xml'. Thank you for using this software.

