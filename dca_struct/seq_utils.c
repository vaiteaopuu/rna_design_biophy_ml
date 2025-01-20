#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "seq_utils.h"
#define MAX_CHAR 100000

// convert_aa string to enum amino acids
Nucleotide convert_aa(char x) {
  char nuc_char[] = {'A', 'C', 'U', 'G', '-'};
  int i;
  for (i = ADE; i <= GAP; i++) {
    if (nuc_char[i] == x)
      return i;
  }
  printf("Error in the AA sequence %c\n", x);
  exit(1);
};

Nucleotide* convert_seq(char* seq) {
  int i;
  int len_seq = strlen(seq);
  Nucleotide* aa_seq = (Nucleotide*)malloc(sizeof(Nucleotide)*len_seq);
  for (i = 0; i < len_seq; i++)
    aa_seq[i] = convert_aa(seq[i]);

  return aa_seq;
}

char* convert_seq2char(Nucleotide* seq, int len_seq) {
  int i;
  char nuc_char[] = {'A', 'C', 'U', 'G', '-'};
  char* aa_seq = (char*)malloc(sizeof(char)*len_seq+1);
  for (i = 0; i < len_seq; i++) {
    aa_seq[i] = nuc_char[seq[i]];
  }
  aa_seq[len_seq] = '\0';
  return aa_seq;
}

// read a fasta file
int get_nb_seq(char* seq_file_str) {
  FILE *seq_file;
  char line[MAX_CHAR];
  int nb_seq = 0;
  // initialize chars

  if ((seq_file = fopen(seq_file_str, "r")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }

  while(fgets(line, MAX_CHAR, seq_file) != NULL) {
    nb_seq++;
  }

  return nb_seq;
}

// read a fasta file
msa_el_p read_msa(char* seq_file_str) {
  FILE *seq_file;
  char line[MAX_CHAR];
  char** msa;
  char* seq[MAX_CHAR];
  int i = 0;
  float weigt = 0.0;
  // initialize chars
  int nb_seq = get_nb_seq(seq_file_str), len_seq;
  msa_el_p msa_data = (msa_el_p)malloc(sizeof(msa_el));

  if ((seq_file = fopen(seq_file_str, "r")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }
  msa = (Nucleotide**)malloc(sizeof(Nucleotide*) * nb_seq);

  while (fgets(line, MAX_CHAR, seq_file) != NULL) {
    if (sscanf(line,"%f %s", &weigt, seq) == 2) {
      /* printf("%s %f\n", seq, weigt); */
    } else if (sscanf(line,"%s", seq) == 1) {
      weigt = 1.0;
      /* printf("%s %f\n", seq, weigt); */
    }
    line[strcspn(line, "\n")] = 0;
    msa[i] = convert_seq(seq);
    i++;
    len_seq = strlen(line);
  }
  msa_data->len_seq = len_seq;
  msa_data->nb_seq = nb_seq;
  msa_data->seq = msa;

  fclose(seq_file);
  return msa_data;
}

void free_msa(msa_el_p msa_data) {
  int i;
  for (i = 0; i < msa_data->nb_seq; i++)
    free(msa_data->seq[i]);
  free(msa_data);
}

double**** create_score(int len_seq) {
  int i, j, ai, aj;
  double**** bias = (double****)malloc(sizeof(double***) * len_seq);
  for (i = 0; i < len_seq; i++) {
    bias[i] = (double***)malloc(sizeof(double**) * len_seq);
    for (j = 0; j < len_seq; j++) {
      bias[i][j] = (double**)malloc(sizeof(double*) * (GAP + 1));
      for (ai = ADE; ai <= GAP; ai++) {
        bias[i][j][ai] = (double*)malloc(sizeof(double) * (GAP + 1));
        for (aj = ADE; aj <= GAP; aj++) {
          bias[i][j][ai][aj] = 0.0;
        }
      }
    }
  }
  return bias;
}

void free_score(double**** bias, int len_seq) {
  int i, j, ai;
  for (i = 0; i < len_seq; i++) {
    for (j = 0; j < len_seq; j++) {
      for (ai = ADE; ai <= GAP; ai++) {
        free(bias[i][j][ai]);
      }
      free(bias[i][j]);
    }
    free(bias[i]);
  }
  free(bias);
}

// otemp is the temperature to output
void write_score(FILE *log_bias, double**** bias, int len_seq, double factor, Nucleotide* ref_seq, int id_step) {
  char nuc_char[] = {'A', 'C', 'U', 'G', '-'};
  int i, j, ai, aj;
  fprintf(log_bias, "# STEP %d\n", id_step);
  for (i = 0; i < len_seq; i++) {
    for (ai = ADE; ai <= GAP; ai++) {
      /* fprintf(log_bias, "%d %d %c %c %f\n", i, i, nuc_char[ai], nuc_char[ai], (bias[i][i][ai][ai] - bias[i][i][ref_seq[i]][ref_seq[i]]) * factor); */
      fprintf(log_bias, "%d %d %c %c %f\n", i, i, nuc_char[ai], nuc_char[ai], bias[i][i][ai][ai] * factor);
    }
    for (j = i+1; j < len_seq; j++) {
      for (ai = ADE; ai <= GAP; ai++) {
        for (aj = ADE; aj <= GAP; aj++) {
          /* fprintf(log_bias, "%d %d %c %c %f\n", i, j, nuc_char[ai], nuc_char[aj], (bias[i][j][ai][aj] - bias[i][j][ref_seq[i]][ref_seq[j]]) * factor); */
          fprintf(log_bias, "%d %d %c %c %f\n", i, j, nuc_char[ai], nuc_char[aj], bias[i][j][ai][aj] * factor);
        }
      }
    }
  }
}

// otemp is the temperature to output
void print_score(double**** bias, int len_seq, double factor, Nucleotide* ref_seq) {
  char nuc_char[] = {'A', 'C', 'U', 'G', '-'};
  int i, j, ai, aj;
  for (i = 0; i < len_seq; i++) {
    for (ai = ADE; ai <= GAP; ai++) {
      printf("%d %d %c %c %f\n", i, i, nuc_char[ai], nuc_char[ai], bias[i][i][ai][ai] * factor);
    }
    for (j = i+1; j < len_seq; j++) {
      for (ai = ADE; ai <= GAP; ai++) {
        for (aj = ADE; aj <= GAP; aj++) {
          printf("%d %d %c %c %f\n", i, j, nuc_char[ai], nuc_char[aj], bias[i][j][ai][aj] * factor);
        }
      }
    }
  }
}

// return a random integer in [from_int; to_int] interval
int rand_int(int from_int, int to_int) {
  if (from_int < to_int)
    return rand() % (to_int - from_int) + from_int;
  else
    fprintf(stderr,"ERROR");
  exit(1);
  return 0;
}

// return a random float in [min; max] interval
float rand_float(int min, int max) {
  if (max == min) return min;
  else if (min < max) return (float)rand() / RAND_MAX;;
  // return 0 if min > max
  return 0.0;
}

void print_seq(Nucleotide* seq, int len_seq) {
  int i;
  char nuc_char[] = {'A', 'C', 'U', 'G', '-'};
  for (i = 0; i < len_seq; i++) {
    if (seq[i] > GAP)
      printf(stderr, "%i %d\n", i, seq[i]);
  }
  for (i = 0; i < len_seq; i++)
    printf("%c", nuc_char[seq[i]]);
}

// Initialize all the matrix score
double**** read_score(int len_seq, char* dca_score_file) {
  double**** all_dca = (double****)malloc(sizeof(double***) * len_seq);
  FILE *dca_file;
  char line[MAX_CHAR];
  char typeiS, typejS;
  float scoreij;
  int posi, posj, typei, typej, i, j, ai, aj;
  if ((dca_file = fopen(dca_score_file, "r")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }

  for (i = 0; i < len_seq; i++) {
    all_dca[i] = (double***)malloc(sizeof(double**) * len_seq);
    for (j = 0; j < len_seq; j++) {
      all_dca[i][j] = (double**)malloc(sizeof(double*) * 21);
      for (ai = ADE; ai <= GAP; ai++) {
        all_dca[i][j][ai] = (double*)malloc(sizeof(double) * 21);
        for (aj = ADE; aj <= GAP; aj++) {
          all_dca[i][j][ai][aj] = 0.0;
        }
      }
    }
  }

  // loop in the bias file
  while(fgets(line, MAX_CHAR, dca_file) != NULL) {
    if (sscanf(line,"%d %d %c %c %f", &posi, &posj, &typeiS, &typejS, &scoreij) == 5) {
      typei = convert_aa(typeiS); typej = convert_aa(typejS);
      all_dca[posi][posj][typei][typej] = (double)scoreij;
      all_dca[posj][posi][typej][typei] = (double)scoreij;
    }
  }

  fclose(dca_file);
  return all_dca;
};


// from a dot bracket notation, compute a contact matrix
float** get_contact_mat(char* struct_mat, int len_seq) {
  int i, j, pi;
  Bool* tmp_pos = (Bool*)malloc(sizeof(Bool*) * len_seq);
  Bool* tmp_pk_pos = (Bool*)malloc(sizeof(Bool*) * len_seq);

  // init output mat
  float** target_mat = (float**)malloc(sizeof(float**) * len_seq);
  for (i = 0; i < len_seq; i++) {
    target_mat[i] = (float*)malloc(sizeof(float*) * len_seq);
    tmp_pos[i] = false;
    tmp_pk_pos[i] = false;
  }
  for (i = 0; i < len_seq; i++) {
    for (j = 0; j < len_seq; j++) {
      target_mat[i][j] = 0.0;
    }
  }


  for (i = 0; i < len_seq; i++) {
    // search for the opening bracket
    if (struct_mat[i] == '(' || struct_mat[i] == '<')
      tmp_pos[i] = true;

    if (struct_mat[i] == ')' || struct_mat[i] == '>') {
      pi = len_seq - 1;
      while (tmp_pos[pi] == false && pi > 0)
        pi--;
      if (pi > 0) {
        target_mat[pi][i] = 1.0;
        tmp_pos[pi] = false;
      }
    }

    if (struct_mat[i] == '[')
      tmp_pk_pos[i] = true;
    if (struct_mat[i] == ']') {
      pi = len_seq -1;
      while (tmp_pk_pos[pi] == false && pi > 0)
        pi--;
      if (pi > 0) {
        target_mat[pi][i] = 1.0;
        tmp_pk_pos[pi] = false;
      }
    }
  }
  free(tmp_pos);
  free(tmp_pk_pos);
  return target_mat;
}

// read struct from file
char* read_struct_file(char* struct_file_str) {
  FILE *struct_file;
  char line[MAX_CHAR];
  char* str_struct = (char*)malloc(sizeof(char)*MAX_CHAR);
  // initialize chars

  if ((struct_file = fopen(struct_file_str, "r")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }
  while(fgets(line, MAX_CHAR, struct_file) != NULL) {
    if (strcmp(str_struct, "") != 0) {
      str_struct[strcspn(str_struct, "\n")] = 0;
      strcat(str_struct, &line[0]);
    } else {
      strcpy(str_struct, line);
    }
  }
  str_struct[strcspn(str_struct, "\n")] = 0;
  fclose(struct_file);
  return str_struct;
}


// extract the alpha parameter from the parameters file
float extract_alpha_value(const char *filename, char *pattern) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error opening file");
    return 1.0;
  }

  char line[1000];
  float alpha_value = 1.0;
  while (fgets(line, sizeof(line), file)) {
    if (strncmp(line, pattern, strlen(pattern)) == 0) {
      char format[100];
      snprintf(format, sizeof(format), "%s %%f", pattern);
      sscanf(line, format, &alpha_value);
      break;
    }
  }

  fclose(file);
  return alpha_value;
}

// compute the number of mutations between 2 sequences
float get_number_mut(Nucleotide* seq_i, Nucleotide* seq_j, int seq_len) {
  float nb_mut = 0;
  int i;
  for (i=0; i < seq_len; i++) {
    /* nb_mut += (seq_i[i] != seq_j[i] && seq_i[i] != GAP && seq_i[i] != GAP) ? 1.0 : 0.0; */
    nb_mut += (seq_i[i] != seq_j[i]) ? 1.0 : 0.0;
  }
  return nb_mut;
}

// read a structures file
struct_el_p read_struct_list(char* struct_file_str) {
  FILE *struct_file;
  char line[MAX_CHAR];
  char** struct_ens;
  char* struct_l[MAX_CHAR];
  int i = 0;
  float weigt = 0.0;
  // initialize chars
  int nb_struct = get_nb_seq(struct_file_str), len_struct;
  struct_el_p struct_data = (struct_el_p)malloc(sizeof(struct_el));

  if ((struct_file = fopen(struct_file_str, "r")) == NULL) {
    printf("Error! opening file");
    exit(1);
  }
  struct_data->struct_ens = (char**)malloc(sizeof(char*) * nb_struct);

  while (fgets(line, MAX_CHAR, struct_file) != NULL) {
    line[strcspn(line, "\n")] = 0;
    len_struct = strlen(line);
    struct_data->struct_ens[i] = malloc((len_struct + 1) * sizeof(char));
    strncpy(struct_data->struct_ens[i], line, len_struct + 1);
    i++;
  }
  struct_data->len_struct = len_struct;
  struct_data->nb_struct = nb_struct;
  /* struct_data->struct_ens = struct_ens; */
  fclose(struct_file);
  return struct_data;
}


void free_struct_s(struct_el_p struct_data) {
  int i;
  for (i = 0; i < struct_data->nb_struct; i++)
    free(struct_data->struct_ens[i]);
  free(struct_data);
}
