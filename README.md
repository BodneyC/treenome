TreeNome
========

Repository containing the source code of a new approach to de novo genome assembly.

Where most current assembly techniques and algorithms make extensive use of graphs (overlap and De Bruijn mostly), TreeNome is a tree-based approach which chooses the *most likely* base to follow the sequence before it.

## Compilation

The project is currently in its early stages and is currently platform independent so a simple `make` will do it.

    make

## Usage

Project is too early to have any real usage as yet.

**Command Format**

The program requires a file (preferably in fastq format), this is supplied with the `-f` switch.

    TreeNome -f ./path/to/file


