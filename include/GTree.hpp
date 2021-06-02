/********************************************************************
 * Filename: GTree.hpp [C++ header code]
 *
 * Description: Declaration of GTree class and NodeInfo struct
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
#ifndef _GTREE_
#define _GTREE_

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include "Node.hpp"
#include "Phred33.hpp"
#include "SeqRead.hpp"

#define USAGE_ERROR -1
#define IN_FILE_ERROR -2
#define OUT_FILE_ERROR -3
#define THREAD_ERROR -4

#define RES 100000

// GTree Helpers
namespace GTH {
extern std::vector<SeqRead> seq_reads;
extern double thresh;
extern double start_weights[NBASES];
extern int64_t start_occs[NBASES];
extern int32_t read_len;

char return_label(int label);
float get_new_weight(double curWeight, int64_t occs, char qual);
void remove_double_ending(std::string &doubleString);
template <typename T> std::string val_to_string(T &val);
} // namespace GTH

// GTree class used in construction and destruction of tree
class GTree {
public:
  GTree()
      : root(nullptr), head(0), n_nodes(0), base_paths(""), occ_paths(""),
        tree_string("") {}
  ~GTree() {}

  /* Sequence creation functions */
  void add_to_seq(uint64_t offset, std::string &sequence);
  void follow_branch(Node *node, short ind, std::string &sequence);
  Node *follow_seq(int32_t offset, std::string &sequence,
                   std::vector<int32_t> &paths);
  std::string cont_seq(Node *node, std::string &sequence);
  void reduce_occs_and_weight(std::vector<int32_t> &paths);
  /* Usefull non-essential */
  Node *get_root();
  int64_t get_n_nodes();
  void print_all_paths(short label);
  std::string store_tree(short label);
  void write_tree_to_file(std::ofstream &storeFile);
  signed short highest_threshold(Node *node);
  /* Tree creation functions */
  void init(short i);
  void add_read_one(uint64_t readNum, short offset);
  void balance_node(Node *node);
  int32_t max_depth(Node *node);
  int32_t max_depth_at_threshold(Node *node);
  /* Tree Recreation Functions */
  void resize_vec(int64_t tmp64);
  void write_vec(std::ifstream &inFile);
  void process_ss_string(std::stringstream &ss);

private:
  Node *root;
  std::atomic<int64_t> head, n_nodes;
  std::string base_paths, occ_paths, tree_string;
  omp_lock_t lock;
  std::deque<Node> d_nodes;

  /* Helper functions */
  short count_children(Node *node);
  /* Tree creation functions */
  void create_root(short ind);
  void pot_add_node(Node *node);
  void update_weight(Node *node, char qual);
  void update_weight_and_occs(Node *node, char qual);
  void create_node(Node *node, short label, char qual, uint64_t rN, int offset);
  /* Public->private overloads */
  void print_all_paths(Node *node, int len, short label);
  void storeTree(Node *node);
};

#endif /*_GTREE_*/
