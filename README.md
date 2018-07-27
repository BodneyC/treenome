TreeNome
========

Repository containing the source code of a new approach to de novo genome assembly.

Where most current assembly techniques and algorithms make extensive use of graphs (overlap and De Bruijn mostly), TreeNome is a tree-based approach which chooses the *most likely* base to follow the sequence before it.

## Compilation

The project is compiled with the makefile:

    make
    
I should note that the project uses [TCLAP](http://tclap.sourceforge.net/) for command line parsing, it is included in the repo as it isn't that large and saves the hassle of having to acquire it one's self.

## Usage

TreeNome can currently take a fastq file and produce a set of trees, representing all paths staring with each base, and build a sequence using the algorithm; it should be noted that the accuracy of the algorithm cannot be confirmed at this stage. The program can store these trees to a file and read that file back in to construct the trees once again.

**Command Format**

The program requires a file, either fastq (-f) or TreeNome output file (-l).

The TCLAP `--help` output should explain the rest of the switches.

```
$ ./TreeNome --help
```
```
USAGE:

   ./TreeNome  {-f <string>|-l <string>} [--phred <33|64>] [-t <int>] [-s
               <string>] [-o] [--] [--version] [-h]


Where:

   -f <string>,  --fastqfile <string>
     (OR required)  Input file in fastq format
         -- OR --
   -l <string>,  --trenomefile <string>
     (OR required)  Input file which was previously outputted from TreeNome


   --phred <33|64>
     Phred base of qualities in fastq

   -t <int>,  --threads <int>
     Number of threads to use

   -s <string>,  --storefile <string>
     File in which to store the tree data

   -o,  --stdout
     Print to stdout

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Tree based de novo DNA assembler
```
