#include "tclap/CmdLine.h"

struct CmdArgs {
  std::string f_fn, st_fn, l_fn, ss_fn;
  bool print_to_console, store_file, load_file, analyse;
  int phred_base;

  CmdArgs()
      : f_fn(""), st_fn(""), l_fn(""), ss_fn(""), print_to_console(0),
        store_file(0), load_file(0), analyse(0), phred_base(33) {}
};

signed int in_file_check(std::string);
signed int return_args(int, char **, struct CmdArgs &);
