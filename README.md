TreeNome
========

Repository containing the source code of a new approach to de novo genome assembly.

Where most current assembly techniques and algorithms make extensive use of graphs (overlap and De Bruijn mostly), TreeNome is a tree-based approach which chooses the most likely base to follow the sequence before it.

## Compilation

The project is compiled with the makefile:

    $ make
    
It should be noted that the project requires C++11.

It should also be noted that the project uses [TCLAP](http://tclap.sourceforge.net/) for command line parsing, it is included in the repo as it isn't that large and saves the hassle of having to acquire it one's self.

## Usage

TreeNome can currently take a fastq file and produce a set of trees, representing all paths staring with each base, and build a sequence using the algorithm; it should be noted that the accuracy of the algorithm cannot be confirmed at this stage. The program can store these trees to a file and read that file back in to construct the trees once again.

**Command Format**

The program requires a file, either fastq (-f) or TreeNome output file (-l).

An example command (using the example FASTQ file in `./fastq/`) would be as follows:

    $ ./TreeNome -f ./fastq/example.fastq -s ./data_out/out1.gno 

The TCLAP `--help` output should explain the rest of the switches.

```
$ ./TreeNome --help
```
```
USAGE: 

   ./TreeNome  {-f <string>|-l <string>} [--score-sys <Phred+33|Phred+64
               |Solexa+64>] [--thresh <0 to 1>] [-t <int>] -s <string>
               [--store-tree <string>] [-a] [-o] [--] [--version] [-h]


Where: 

   -f <string>,  --fastqfile <string>
     (OR required)  Input file in fastq format
         -- OR --
   -l <string>,  --load-file <string>
     (OR required)  Input file which was previously outputted from TreeNome


   --score-sys <Phred+33|Phred+64|Solexa+64>
     Scoring system. Phred+33 default (Sanger)

   --thresh <0 to 1>
     Threashold value

   -t <int>,  --threads <int>
     Number of threads to use

   -s <string>,  --store-sequence <string>
     (required)  File in which to store the sequence data

   --store-tree <string>
     File in which to store the tree data

   -a,  --analyse
     Perform tree analysis

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
