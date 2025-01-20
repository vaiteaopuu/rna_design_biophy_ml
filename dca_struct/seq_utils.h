#ifndef SEQ_UTILS_

#include <stdio.h>
#include <string.h>

// Global variable
typedef enum {true, false} Bool;
typedef enum {ADE, CYT, URA, GUA, GAP} Nucleotide;
/* typedef enum {ALA ,CYS ,THR ,GLU ,ASP ,PHE ,TRP ,ILE ,VAL ,LEU ,LYS , */
/*               MET ,ASN ,GLN ,SER ,ARG ,TYR ,HIS ,PRO ,GLY, GAP */
/* } Aminoacid; */

// Proposed mutation in MCMC
typedef struct {
  int nb_seq;
  int len_seq;
  // proposed types
  Nucleotide** seq;
  // the proposed nucleotide
} msa_el;
typedef msa_el *msa_el_p;

// struct ensemble
typedef struct {
  int nb_struct;
  int len_struct;
  // proposed types
  char** struct_ens;
  // the proposed nucleotide
} struct_el;
typedef struct_el *struct_el_p;

// Functions
msa_el_p read_msa(char* seq_file_str);

double**** create_score(int len_seq);

void free_msa(msa_el_p msa_data);

void free_score(double**** bias, int len_seq);

// reverse a string
int rand_int(int from_int, int to_int);

// random float
float rand_float(int min, int max);

// print one seq
void print_seq(Nucleotide* seq, int len_seq);

// print bias
void print_score(double**** bias, int len_seq, double factor, Nucleotide* ref_seq);

// read the scores
double**** read_score(int len_seq, char* dca_score_file);

// convert to seq
Nucleotide* convert_seq(char* seq);

// write log bias
void write_score(FILE *log_bias, double**** bias, int len_seq, double factor, Nucleotide* ref_seq, int id_step);

char* convert_seq2char(Nucleotide* seq, int len_seq);

float** get_contact_mat(char* struct_mat, int len_seq);

char* read_struct_file(char* struct_file_str);

float extract_alpha_value(const char *filename, char* pattern);

float get_number_mut(Nucleotide* seq_i, Nucleotide* seq_j, int seq_len);

struct_el_p read_struct_list(char* struct_file_str);

void free_struct_s(struct_el_p struct_data);
#endif // SEQ_UTILS_
