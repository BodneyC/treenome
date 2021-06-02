/********************************************************************
 * Filename: TreeNome.cpp [C++ source code]
 *
 * Description: Entry point for TreeNome de novo DNA assembler
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Details: Usage of TCLAP for command line args;
 *		recording of times/memory usages in analysis()
 *
 *******************************************************************/
#include <sys/sysinfo.h>
#include <sys/types.h>

#include "InputFile.hpp"
#include "TreeTop.hpp"
#include "cli.hpp"

#define OUT

void write_trees_to_file(std::string st_fn, TreeTop *tree_top) {
  std::ofstream out_file(st_fn);
  for (int i = 0; i < NBASES; i++)
    out_file << tree_top->tree_strings[i].c_str() << std::endl;
}

TreeTop *create_tree_from_reads(CmdArgs &cmd_args, OUT double &t2) {
  InputFile in_file(cmd_args.f_fn, cmd_args.phred_base);
  if (!in_file.read_fastq()) {
    std::cout << "[ERR]: Input file could not be opened" << std::endl;
    return nullptr;
  } else {
    std::cout << in_file.n_reads << " records found\n" << std::endl;
  }

  TreeTop *tree_top = new TreeTop;

  double t1 = omp_get_wtime();
  tree_top->process_reads_one();
  t2 = omp_get_wtime() - t1;

  return tree_top;
}

TreeTop *loadTreeFromFile(CmdArgs &argList) {
  TreeTop *treeTop = new TreeTop;
  treeTop->reconstruct_trees(argList.l_fn);

  return treeTop;
}

int64_t parseLine(std::string &line) {
  size_t position_of_k = line.find_last_of('k');
  if (position_of_k == std::string::npos) {
    return -1;
  }

  int j = 0;
  while (line[j] < '0' || line[j] > '9') {
    j++;
  }

  std::string val = line.substr(j, position_of_k - j);

  int64_t ret = std::stoll(val, nullptr, 10);

  return ret;
}

int64_t get_mem_in_use() {
  std::ifstream self("/proc/self/status");
  if (!self.is_open())
    return IN_FILE_ERROR;

  int result = -1;
  std::string line;

  while (getline(self, line, '\n')) {
    if (line.substr(0, 7) == "VmSize:") {
      result = parseLine(line);
      break;
    }
  }

  self.close();
  return result;
}

// https://stackoverflow.com/questions/7276826/c-format-number-with-commas
template <typename T> std::string val_w_commas(T val) {
  std::string st = std::to_string(val);
  int i = st.find_last_of(".");
  if (i == -1) {
    i = st.length() - 3;
  } else {
    i -= 3;
    st.erase(st.find_last_not_of('0') + 1, std::string::npos);
    if (st[st.length() - 1] == '.')
      st = st.substr(0, st.length() - 1);
  }

  while (i > 0) {
    st.insert(i, ",");
    i -= 3;
  }

  return st;
}

signed int analysis(TreeTop *tree_top, double time_to_construct,
                    double time_to_build_seq) {
  std::cout
      << "------------------------------------------\nTiming Information:\n"
      << std::endl
      << "Number of threads in use       : " << NUM_THREADS << std::endl
      << "Time to construct trees  ( s ) : " << time_to_construct << std::endl
      << "Time to build sequence   ( s ) : " << time_to_build_seq << '\n'
      << std::endl;

  struct sysinfo mem_info;
  sysinfo(&mem_info);

  int64_t total_virt_mem = mem_info.totalram;
  total_virt_mem *= mem_info.mem_unit;
  double total_vmem_kb = static_cast<double>(total_virt_mem) / 1024.0;

  int64_t vmem_in_use = get_mem_in_use();
  if (vmem_in_use == -1)
    return vmem_in_use;

  std::cout
      << "------------------------------------------\nMemory Information:\n"
      << std::endl
      << "Virtual memory in use   ( KB ) : " << val_w_commas(vmem_in_use)
      << std::endl
      << "Total available         ( KB ) : " << val_w_commas(total_vmem_kb)
      << std::endl
      << "Percent used             ( % ) : "
      << static_cast<double>(vmem_in_use) / total_vmem_kb << std::endl;

  return 0;
}

int main(int argc, char **argv) {
  CmdArgs cmd_args;
  signed int ret_code = return_args(argc, argv, cmd_args);

  switch (ret_code) {
  case USAGE_ERROR:
    std::cout << "[ERR]: Incorrect command usage" << std::endl;
    break;
  case IN_FILE_ERROR:
    std::cout << "[ERR]: Input file not supplied or not present" << std::endl;
    break;
  case THREAD_ERROR:
    std::cout << "[ERR]: " << NUM_THREADS << " threads requested, "
              << omp_get_max_threads() << " available." << std::endl;
    break;
  }
  if (ret_code)
    return ret_code;

  omp_set_num_threads(NUM_THREADS);

  double time_to_construct;

  TreeTop *tree_top;
  if (cmd_args.load_file)
    tree_top = loadTreeFromFile(cmd_args);
  else
    tree_top = create_tree_from_reads(cmd_args, time_to_construct);

  if (!tree_top)
    return FILE_ERROR;

  if (cmd_args.store_file) {
    ret_code = tree_top->store_trees(cmd_args.st_fn);
    if (ret_code)
      return ret_code;
  }

  // Has to be before buildSequence()
  if (cmd_args.print_to_console)
    tree_top->print_tree();

  // Has to be before buildSequence()
  tree_top->analyse_trees();

  double t1 = omp_get_wtime();
  tree_top->build_seq();
  double time_to_build_seq = omp_get_wtime() - t1;

  ret_code = tree_top->store_seq(cmd_args.ss_fn);
  if (ret_code)
    return ret_code;

  if (cmd_args.print_to_console)
    tree_top->print_seq();

  if (cmd_args.analyse) {
    analysis(tree_top, time_to_construct, time_to_build_seq);
    tree_top->print_analysis();
  }

  delete tree_top;

  return 0;
}
