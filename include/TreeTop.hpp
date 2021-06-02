/********************************************************************
 * Filename: TreeTop.hpp [C++ header code]
 *
 * Description: Header for forest class containing GTrees
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
#ifndef _TREETOP_
#define _TREETOP_

#include "GTree.hpp"

// Structure which holds read data and the four GTrees
class TreeTop {
public:
  std::string tree_strings[NBASES];

  TreeTop() : sequence(""), n_reads(0), read_length(0) {}

  /* Common */
  void build_seq();
  void print_tree();
  signed int store_trees(std::string &s_fn);
  void print_seq();
  signed int store_seq(std::string &o_fn);
  void analyse_trees();
  void print_analysis();
  /* FReads */
  void process_reads_one();
  /* FFile */
  signed int reconstruct_trees(std::string &i_fn);

private:
  std::string sequence;
  int64_t n_reads;
  int read_length;
  GTree trees[NBASES];
  int32_t node_cnt[NBASES], max_depth[NBASES], max_deptch_at_threshold[NBASES];

  short max_path();
  bool root_occs_above_threshold();
  bool root_occs_exist();
  /* FReads */
  void create_roots();
  void f_thread(uint64_t i);
};

#endif /*_TREETOP_*/
