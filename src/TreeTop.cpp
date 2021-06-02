/********************************************************************
 * Filename: TreeTop.cpp [C++ source code]
 *
 * Description: Implementation of TreeTop class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Details: Where read processing and the building of the sequence
 *		takes places; also the analysis printing functions
 *
 *******************************************************************/
#include "TreeTop.hpp"

#define MER_LEN 3

int NUM_THREADS;

/** ------------- Sequence Generation -------------- **/
bool TreeTop::root_occs_exist() {
  short ret = 0;

  for (int i = 0; i < NBASES; i++)
    if (trees[i].get_root()->occs == 0)
      ret++;

  return ret == NBASES ? 0 : 1;
}

bool TreeTop::root_occs_above_threshold() {
  short ret = 0;

  for (int i = 0; i < NBASES; i++) {
    if (GTH::start_occs[i]) {
      if (GTH::start_weights[i] / static_cast<double>(GTH::start_occs[i]) <
          GTH::thresh) {
        ret++;
      }
    } else if (!GTH::start_occs[i]) {
      ret++;
    }
  }

  return ret == NBASES ? 0 : 1;
}

short TreeTop::max_path() {
  int start = 0;
  double max_ratio = GTH::thresh;

  for (short i = 0; i < NBASES; i++) {
    if (GTH::start_occs[i]) {
      double curRat =
          GTH::start_weights[i] / static_cast<double>(GTH::start_occs[i]);
      if (curRat >= max_ratio) {
        max_ratio = curRat;
        start = i;
      }
    }
  }

  GTH::start_weights[start]--;
  GTH::start_occs[start]--;
  trees[start].follow_branch(trees[start].get_root(), start, sequence);

  return start;
}

void TreeTop::build_seq() {
  max_path();

  uint64_t offset = 1;
  while (root_occs_above_threshold()) {

    // Possibly make is a tighter gap as its working from single letters
    // ( this would actually be the k-mer match )
    if (offset >= sequence.length() - MER_LEN) {
      sequence += 'N';
      max_path();
      offset += MER_LEN + 1;
    }

    if (trees[BASE_IND(sequence[offset])].get_root()->get_ratio() >=
        GTH::thresh) {
      trees[BASE_IND(sequence[offset])].add_to_seq(offset, sequence);
    }
    offset++;
  }
}

signed int TreeTop::store_seq(std::string &o_fn) {
  std::ofstream store_file(o_fn);
  if (!store_file.is_open()) {
    return OUT_FILE_ERROR;
  }

  store_file << sequence;
  store_file.close();
  return 0;
}

/** ------------ Tree Reconstruction --------------- **/
signed int TreeTop::store_trees(std::string &sFilename) {
  std::ofstream store_file(sFilename, std::ios::binary);
  if (!store_file.is_open()) {
    return OUT_FILE_ERROR;
  }

  for (int i = 0; i < NBASES; i++) {
    store_file.write((char *)&GTH::start_occs[i], sizeof(int64_t));
    store_file.write((char *)&GTH::start_weights[i], sizeof(double));
  }

  for (int i = 0; i < NBASES; i++) {
    int64_t tmp_64 = trees[i].get_n_nodes();
    store_file.write((char *)&tmp_64, sizeof(int64_t));
  }

  for (int i = 0; i < NBASES; i++) {
    trees[i].write_tree_to_file(store_file);
  }

  store_file.close();

  return 0;
}

signed int TreeTop::reconstruct_trees(std::string &iFilename) {
  std::ifstream in_file(iFilename, std::ios::binary);
  if (!in_file.is_open())
    return IN_FILE_ERROR;

  for (int i = 0; i < NBASES; i++) {
    in_file.read((char *)&GTH::start_occs[i], sizeof(int64_t));
    in_file.read((char *)&GTH::start_weights[i], sizeof(double));
  }

  for (int i = 0; i < NBASES; i++) {
    int64_t tmp64;
    in_file.read((char *)&tmp64, sizeof(int64_t));
    trees[i].resize_vec(tmp64);
  }

  for (int i = 0; i < NBASES; i++) {
    trees[i].write_vec(in_file);
  }

  in_file.close();

  return 0;
}

/** --------------- Read Processing ---------------- **/
void TreeTop::f_thread(uint64_t i) {
  for (short j = 0; j < GTH::seq_reads[i].size(); j++) {
    // if( GTH::seqReads[i].getBaseInd( j ) == 1 )
    trees[GTH::seq_reads[i].get_base_idx(j)].add_read_one(i, j);
  }
}

void TreeTop::process_reads_one() {
  for (int i = 0; i < NBASES; i++) {
    trees[i].init(i);
  }

#pragma omp parallel num_threads(NUM_THREADS)
  {
    for (uint64_t i = 0; i < GTH::seq_reads.size(); i += NUM_THREADS) {
#pragma omp for schedule(static, 1)
      for (int j = 0; j < NUM_THREADS; j++) {
        if (i + j < GTH::seq_reads.size()) {
          f_thread(i + j);
        }
      }
    }
  }
}

/** --------------- Misc Functions ----------------- **/
char get_base(int index) {
  switch (index) {
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'G';
  case 3:
    return 'T';
  default:
    return 'N';
  }
}

void TreeTop::analyse_trees() {
  for (int i = 0; i < NBASES; i++) {
    node_cnt[i] = trees[i].get_n_nodes();
    max_depth[i] = trees[i].max_depth(trees[i].get_root());
    max_deptch_at_threshold[i] =
        trees[i].max_depth_at_threshold(trees[i].get_root());
  }
}

void TreeTop::print_analysis() {
  std::cout
      << "\n------------------------------------------\nTree Information:\n"
      << std::endl;

  for (int i = 0; i < NBASES; i++) {
    std::cout << "- " << get_base(i) << "-tree:" << std::endl
              << "  Node count                   : " << node_cnt[i] << std::endl
              << "  Maximum depth of tree        : " << max_depth[i]
              << std::endl
              << "  Max depth above threshold    : "
              << max_deptch_at_threshold[i] << '\n'
              << std::endl;
  }
}

void TreeTop::print_tree() {
  for (int i = 0; i < NBASES; i++) {
    trees[i].print_all_paths(i);
  }
}

void TreeTop::print_seq() {
  // 120 for terminal width's sake
  unsigned short term_width = 120;
  uint64_t i = 0;

  if (sequence.length() > term_width)
    for (; i < sequence.length() - term_width; i += term_width)
      std::cout << sequence.substr(i, term_width) << std::endl;

  std::cout << sequence.substr(i, sequence.length() - i) << std::endl;
}
