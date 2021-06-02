/********************************************************************
 * Filename: GTree.cpp [C++ source code]
 *
 * Description: Implementation of GTree class, templatized for
 *		Node/pNode purposes
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
#include "GTree.hpp"

#ifdef __MINGW64__
namespace mingw_fix {
{
  std::ostringstream ss;
  ss << val;
  return ss.str();
} // namespace mingw_fix
} // namespace mingw_fix
#endif /*__MINGW32__*/

/** --------------- Helper Functions --------------- **/
namespace GTH {
std::vector<SeqRead> seq_reads;
double thresh;
double start_weights[NBASES];
int64_t start_occs[NBASES];
int score_sys;
int32_t read_len;

char return_label(int label) {
  switch (label) {
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'T';
  case 3:
    return 'G';
  default: // 7 but it shouldn't
    return 'P';
  }
}

void remove_double_ending(std::string &doub_str) {
  doub_str.erase(doub_str.find_last_not_of('0') + 1, std::string::npos);
  doub_str.erase(doub_str.find_last_not_of('.') + 1, std::string::npos);
}

template <typename T> std::string val_to_string(T &val) {
#ifdef __MINGW64__
  return mingw_fix::to_string(val);
#else
  return std::to_string(val);
#endif /*__MINGW64__*/
}
} // namespace GTH

/** ------------- - Helper Functions --------------- **/
Node *GTree::get_root() {
  if (root)
    return root;
  else
    return nullptr;
}

int64_t GTree::get_n_nodes() { return n_nodes; }

short GTree::count_children(Node *node) {
  short children = 0;

  for (int i = 0; i < NBASES; i++)
    if (node->subnodes[i])
      children++;

  return children;
}

signed short GTree::highest_threshold(Node *node) {
  short ind = -1;
  double max_ratio = GTH::thresh;

  for (short i = 0; i < NBASES; i++) {
    Node *tmp_node = &d_nodes[node->subnodes[i]];
    if (node->subnodes[i] && tmp_node->occs) {
      double cur_ratio = tmp_node->get_ratio();
      if (cur_ratio > max_ratio) {
        max_ratio = cur_ratio;
        ind = i;
      }
    }
  }

  return ind;
}

/** --------------- Genome Creation ---------------- **/
void GTree::reduce_occs_and_weight(std::vector<int32_t> &paths) {
  for (uint32_t i = 0; i < paths.size(); i++) {
    d_nodes[paths[i]].occs--;
    d_nodes[paths[i]].weight = d_nodes[paths[i]].weight.load() - 1;
  }
}

Node *GTree::follow_seq(int32_t offset, std::string &sequence,
                        std::vector<int32_t> &paths) {
  Node *node = root;
  paths.push_back(0);

  // Reach the end of the existing path
  for (uint32_t i = offset; i < sequence.length(); i++) {
    short ind = BASE_IND(sequence[i]);
    if (node->subnodes[ind]) {
      paths.push_back(node->subnodes[ind]);
      node = &d_nodes[node->subnodes[ind]];
    } else {
      return nullptr;
    }
  }

  if (highest_threshold(node) == -1)
    return nullptr;

  return node;
}

void GTree::follow_branch(Node *node, short ind, std::string &sequence) {
  short children;
  do {
    children = count_children(node);
    sequence += GTH::return_label(ind);
    ind = highest_threshold(node);
    node->occs--;
    node->weight = node->weight.load() - 1;
    if (ind == -1)
      return;
    node = &d_nodes[node->subnodes[ind]];
  } while (children);
}

void GTree::add_to_seq(uint64_t offset, std::string &sequence) {
  Node *node = root;
  short ind = 0;
  std::vector<Node *> n_path;

  // End of current sequence
  for (uint32_t i = offset + 1; i < sequence.length(); i++) {
    ind = BASE_IND(sequence[i]);
    if (node->subnodes[ind] && count_children(&d_nodes[node->subnodes[ind]])) {
      n_path.push_back(node);
      node = &d_nodes[node->subnodes[ind]];
    } else {
      return;
    }
  }

  ind = highest_threshold(node);
  if (ind != -1) {
    follow_branch(&d_nodes[node->subnodes[ind]], ind, sequence);

    // Only if something is contributed to the sequence should the
    // occurences be lowered
    for (uint32_t i = 0; i < n_path.size(); i++) {
      n_path[i]->occs--;
      n_path[i]->weight = n_path[i]->weight.load() - 1;
    }
  }
}

/** --------------- Tree from Reads ---------------- **/
void GTree::init(short i) {
  // A deque has been used to prevent frequent reallocations of existing
  // vector data when dNodes is resized
  omp_init_lock(&lock);
  d_nodes.resize(RES);
  create_root(i);
}

/** ---------------- Tree Creation ----------------- **/
void GTree::update_weight(Node *node, char qual) {
  double new_weight, cur_weight = node->weight;
  double p_b_quality =
      GTH::phred_qualities[GTH::score_sys][static_cast<int>(qual)];
  do {
    new_weight = cur_weight + p_b_quality;
  } while (!(node->weight.compare_exchange_weak(cur_weight, new_weight)));
}

void GTree::update_weight_and_occs(Node *node, char qual) {
  node->occs++;
  update_weight(node, qual);
}

void GTree::create_root(short ind) {
  char qual;
  root = &(d_nodes[0]);
  n_nodes++;

  for (int64_t i = 0; i < GTH::seq_reads.size(); i++) {
    int offset = -1;
    // Find algorithm in bool vec
    for (int64_t j = 0; j < GTH::seq_reads[i].size(); j++) {
      if (GTH::seq_reads[i].get_base_idx(j) == ind) {
        offset = j;
        break;
      }
    }
    if (offset != -1) {
      qual = GTH::seq_reads[i].get_quality(offset);
      root->read_num = i;
      root->offset = offset;
      root->weight =
          GTH::phred_qualities[GTH::score_sys][static_cast<int>(qual)];
      root->occs = 1;
      break;
    }
  }
}

void GTree::create_node(Node *node, short ind, char qual, uint64_t rN,
                        int offset) {
  omp_set_lock(&lock);
  head++;
  if ((unsigned)head == d_nodes.size()) {
    d_nodes.resize(d_nodes.size() + RES);
  }
  node->subnodes[ind] = head;
  omp_unset_lock(&lock);
  n_nodes++;

  // Doesn't set occs as addRead...() will do that
  Node *tmpNode = &d_nodes[node->subnodes[ind]];
  tmpNode->read_num = rN;
  tmpNode->offset = offset;
  tmpNode->weight =
      GTH::phred_qualities[GTH::score_sys][static_cast<int>(qual)];
  tmpNode->occs = 1;
}

/** --------------- Read Processing ---------------- **/
void GTree::add_read_one(uint64_t readNum, short offset) {
  std::vector<Node *> paths;

  Node *node = root;
  SeqRead *read = &GTH::seq_reads[readNum];
  bool ret_bool = 0, update_bool = 0;

  // Edge case
  if (root->read_num == readNum && root->offset == offset)
    return;

  for (int i = offset + 1; i < GTH::seq_reads[readNum].size(); i++) {
    short ind = (*read).get_base_idx(i);
    paths.push_back(node);

    omp_set_lock(&node->lock);
    int64_t next_node_location = node->subnodes[ind];
    if (!next_node_location) {
      // std::cout << ( *read ).getQual( i ) << std::endl;
      create_node(node, ind, (*read).get_quality(i), readNum, i);
      ret_bool = update_bool = 1;
      if (count_children(node) == 1)
        // Even if it doesn't balance it, a node was created and so
        // the path needs updating
        balance_node(node);
    }
    if (i + 1 == GTH::seq_reads[readNum].size()) {
      // If the end is reached, they still count as occurrences of the path...
      update_bool = 1;
      // If a subnode exists and it wasn't created above
      if (next_node_location && !ret_bool)
        paths.push_back(&d_nodes[next_node_location]);
      // For balancing purposes
      if (next_node_location && !count_children(&d_nodes[next_node_location]))
        pot_add_node(&d_nodes[next_node_location]);
    }
    omp_unset_lock(&node->lock);

    // It will always enter updateBool, just needs to know when
    if (update_bool)
      for (unsigned int j = 0, k = offset; j < paths.size(); j++, k++)
        update_weight_and_occs(paths[j], (*read).get_quality(k));
    if (ret_bool)
      break;

    node = &d_nodes[next_node_location];
  }
  paths.clear();
}

void GTree::balance_node(Node *node) {
  // Get the offset and read before overiding/updating
  uint64_t left_read_num = node->read_num;
  SeqRead *left_read = &GTH::seq_reads[left_read_num];
  short left_offset = node->offset + 1;
  bool clear = 0;

  // If there is nothing to balance with:
  // (will obviously cause imbalanced weights/occs)
  if (left_offset == (*left_read).size())
    return;

  short left_idx = (*left_read).get_base_idx(left_offset);
  char left_quality = (*left_read).get_quality(left_offset);

  if (node->subnodes[left_idx] &&
      left_offset == d_nodes[node->subnodes[left_idx]].offset &&
      left_read_num == d_nodes[node->subnodes[left_idx]].read_num)
    return;

  // If the paths are different:
  if (!node->subnodes[left_idx]) {
    create_node(node, left_idx, left_quality, left_read_num, left_offset);
    return;
  }
  node = &d_nodes[node->subnodes[left_idx]];
  update_weight_and_occs(node, left_quality);

  // If the paths follow the same route:
  uint64_t right_read_num = node->read_num;
  SeqRead *right_read = &GTH::seq_reads[right_read_num];
  short right_offset = node->offset + 1;
  left_offset++;

  std::vector<short> left_idx_path;
  std::vector<short> right_idx_path;
  short tmp_left_offset = left_offset, tmp_right_offset = right_offset;

  while (left_offset < (*left_read).size() &&
         right_offset < (*right_read).size()) {
    left_idx_path.push_back((*left_read).get_base_idx(left_offset));
    right_idx_path.push_back((*right_read).get_base_idx(right_offset));
    if ((*right_read).get_base_idx(right_offset) !=
        (*left_read).get_base_idx(left_offset)) {
      clear = 1;
      break;
    }
    left_offset++;
    right_offset++;
  }

  // Update to the end of shared bases
  for (unsigned short i = 0; i < left_idx_path.size() - clear; i++) {
    create_node(node, left_idx_path[i],
                (*left_read).get_quality(tmp_left_offset), left_read_num,
                tmp_left_offset);
    node = &d_nodes[node->subnodes[left_idx_path[i]]];
    update_weight_and_occs(node, (*right_read).get_quality(tmp_right_offset));
    tmp_left_offset++;
    tmp_right_offset++;
  }

  // If they differ, create a node each
  if (clear) {
    create_node(node, left_idx_path.back(),
                (*left_read).get_quality(tmp_left_offset), left_read_num,
                tmp_left_offset);
    create_node(node, right_idx_path.back(),
                (*right_read).get_quality(tmp_right_offset), right_read_num,
                tmp_right_offset);
    return;
  }

  // If one ends, create the relevant node
  if (left_offset < (*left_read).size()) {
    left_idx = (*left_read).get_base_idx(tmp_left_offset);
    left_quality = (*left_read).get_quality(tmp_left_offset);
    create_node(node, left_idx, left_quality, left_read_num, tmp_left_offset);
  }
  if (right_offset < (*right_read).size()) {
    short right_idx = (*right_read).get_base_idx(tmp_right_offset);
    char right_qual = (*right_read).get_quality(tmp_right_offset);
    create_node(node, right_idx, right_qual, right_read_num, tmp_right_offset);
  }

  return;
}

void GTree::pot_add_node(Node *node) {
  int64_t read_num = node->read_num;
  SeqRead *read = &GTH::seq_reads[read_num];
  short offset = node->offset + 1;

  if (offset == (*read).size())
    return;

  short ind = (*read).get_base_idx(offset);
  char qual = (*read).get_quality(offset);

  create_node(node, ind, qual, read_num, offset);
}

/** ---------------- Tree Storage ------------------ **/
void GTree::write_tree_to_file(std::ofstream &store_file) {
  for (uint64_t i = 0; i < d_nodes.size(); i++)
    store_file.write((char *)&d_nodes[i], sizeof(Node));
}

/** --------------- Tree From File ----------------- **/
void GTree::resize_vec(int64_t tmp_64) {
  n_nodes = tmp_64;
  if (tmp_64 % RES)
    tmp_64 = (tmp_64 / RES + 1) * RES;
  // std::cout << tmp64 << std::endl;
  d_nodes.resize(tmp_64);
}

void GTree::write_vec(std::ifstream &in_file) {
  for (int64_t i = 0; i < d_nodes.size(); i++) {
    Node tmp_node;
    in_file.read((char *)&tmp_node, sizeof(Node));
    d_nodes[i] = tmp_node;
  }
  root = &d_nodes[0];
}

/** -------------- Useful Functions ---------------- **/
int32_t GTree::max_depth(Node *node) {
  int32_t depth[NBASES] = {-1};
  int32_t ret = 0;

  for (int i = 0; i < NBASES; i++) {
    if (node->subnodes[i]) {
      ret = -1;
      depth[i] = max_depth(&d_nodes[node->subnodes[i]]);
    }
  }

  if (!ret)
    return 1;

  for (int i = 0; i < NBASES; i++)
    if (depth[i] > ret)
      ret = depth[i];

  return ret + 1;
}

int32_t GTree::max_depth_at_threshold(Node *node) {
  int32_t depth[NBASES] = {-1};
  int32_t ret = 0;

  if (node->get_ratio() < GTH::thresh)
    return 0;

  for (int i = 0; i < NBASES; i++) {
    if (node->subnodes[i]) {
      ret = -1;
      depth[i] = max_depth_at_threshold(&d_nodes[node->subnodes[i]]);
    }
  }

  if (!ret)
    return 1;

  for (int i = 0; i < NBASES; i++)
    if (depth[i] > ret)
      ret = depth[i];

  return ret + 1;
}

void GTree::print_all_paths(short label) {
  print_all_paths(root, 0, label);
  base_paths = occ_paths = "";
}

void GTree::print_all_paths(Node *node, int len, short label) {
  occ_paths.erase(len, occ_paths.length());
  std::string val = GTH::val_to_string(node->occs);
  occ_paths += val;
  occ_paths += "-";
  base_paths.erase(len, base_paths.length());
  base_paths += GTH::return_label(label);
  for (uint32_t i = 0; i < val.length(); i++)
    base_paths += '-';
  len += val.length() + 1;
  bool check = 0;
  for (int i = 0; i < NBASES; i++)
    if (node->subnodes[i])
      check = true;
  // std::cout << "WEIGHT: " << node->weight << '\n';
  if (!check) {
    occ_paths.erase(occ_paths.length() - 1);
    std::cout << occ_paths << ": EOS" << std::endl;
    base_paths.erase(base_paths.length() - 1);
    std::cout << base_paths << ": EOS\n" << std::endl;
    return;
  }
  for (int i = 0; i < NBASES; i++)
    if (node->subnodes[i])
      print_all_paths(&d_nodes[node->subnodes[i]], len, i);
}
