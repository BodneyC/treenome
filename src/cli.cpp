#include "cli.hpp"
#include "GTree.hpp"

signed int in_file_check(std::string filename) {
  std::ifstream testFile(filename);
  if (testFile.good())
    return 0;
  else
    return IN_FILE_ERROR;
}

signed int return_args(int argc, char **argv, struct CmdArgs &cmd_args) {
  try {
    TCLAP::CmdLine cmd("Tree based de novo DNA assembler", ' ', "1.04");
    TCLAP::ValueArg<std::string> f_fn_arg(
        "f", "fastqfile", "Input file in fastq format", false, "", "string");
    TCLAP::ValueArg<std::string> l_fn_arg(
        "l", "load-file",
        "Input file which was previously outputted from TreeNome", false, "",
        "string");
    TCLAP::ValueArg<std::string> st_fn_arg(
        "", "store-tree", "File in which to store the tree data", false, "",
        "string");
    TCLAP::ValueArg<std::string> ss_fn_arg(
        "s", "store-sequence", "File in which to store the sequence data", true,
        "", "string");
    TCLAP::ValueArg<int> thread_arg("t", "threads", "Number of threads to use",
                                    false, 1, "int");
    TCLAP::ValueArg<double> thresh_arg("", "thresh", "Threashold value", false,
                                       0.3f, "0 to 1");
    std::vector<std::string> p_allowed = {"Phred+33", "Phred+64", "Solexa+64"};
    TCLAP::ValuesConstraint<std::string> p_allowed_vc(p_allowed);
    TCLAP::ValueArg<std::string> score_arg(
        "", "score-sys", "Scoring system. Phred+33 default (Sanger)", false,
        "Phred+33", &p_allowed_vc);
    TCLAP::SwitchArg o_switch("o", "stdout", "Print to stdout", cmd, 0);
    TCLAP::SwitchArg a_switch("a", "analyse", "Perform tree analysis", cmd, 0);

    cmd.xorAdd(f_fn_arg, l_fn_arg);
    cmd.add(st_fn_arg);
    cmd.add(ss_fn_arg);
    cmd.add(thread_arg);
    cmd.add(thresh_arg);
    cmd.add(score_arg);
    cmd.parse(argc, argv);

    if (!f_fn_arg.isSet() && !l_fn_arg.isSet())
      return USAGE_ERROR;

    cmd_args.f_fn = f_fn_arg.getValue();
    if (f_fn_arg.isSet() && in_file_check(cmd_args.f_fn))
      return IN_FILE_ERROR;

    cmd_args.l_fn = l_fn_arg.getValue();
    if (l_fn_arg.isSet()) {
      cmd_args.load_file = 1;
      if (in_file_check(cmd_args.l_fn))
        return IN_FILE_ERROR;
    }

    cmd_args.st_fn = st_fn_arg.getValue();
    if (st_fn_arg.isSet())
      cmd_args.store_file = 1;

    if (ss_fn_arg.isSet()) {
      cmd_args.ss_fn = ss_fn_arg.getValue();
      if (cmd_args.ss_fn.substr(cmd_args.ss_fn.find_last_of(".") + 1) != "gno")
        cmd_args.ss_fn += ".gno";
    }

    NUM_THREADS = thread_arg.getValue();
    if (NUM_THREADS > omp_get_max_threads())
      return THREAD_ERROR;

    GTH::score_sys = 0;
    if (score_arg.getValue() == "Phred+33")
      cmd_args.phred_base = 33;
    else if (score_arg.getValue() == "Phred+64")
      cmd_args.phred_base = 64;
    else if (score_arg.getValue() == "Solexa+64") {
      cmd_args.phred_base = 59;
      GTH::score_sys = 1;
    }

    cmd_args.print_to_console = o_switch.getValue();
    cmd_args.analyse = a_switch.getValue();
    GTH::thresh = thresh_arg.getValue();
  } catch (TCLAP::ArgException &e) {
    std::cerr << "Error: " << e.error() << " for arg " << e.argId()
              << std::endl;
  }
  return 0;
}
