#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <argp.h>
/* #include <omp.h> */
#include "seq_utils.h"
#include "mcmc.h"
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils/basic.h>

#define INF 1.0/0.0

// default mode is standard MCMC
Bool VERBOSE = false;
Bool OPT = false;

// mode for computing the structure contribution
Bool STRUCT_MODE = false;
Bool GAP_MODE = 0;

// bias update frequency
int UP_FREQ = 1;

// bias update parameters
float LAMBDA = 1.0;
float PHI = 0.0;
float ALPHA = 1.0;
float ALPHA_MUT = 0.0;
float HSTEP = 0.1;
float OUT_FREQ = 0.001;
float TEMP = 1.0;

char nuc_char[] = {'A','C','T','G', '-'};

// Parse arguments -------------------------------------------------------------
const char *argp_program_version = "Version 1.0";
const char *argp_program_bug_address = "@";

/* Program documentation. */
static char doc[] = "Standard MCMC procedure";

/* A description of the arguments we accept. */
static char args_doc[] = "<msa> <struct>";

/* The options we understand. */
static struct argp_option options[] = {
  {"nb_steps"    ,   'n' , "INT"   , 0 , "Nb of MC steps"                     , 0 } ,
  {"up_steps"    ,   'u' , "INT"   , 0 , "update frequency"                   , 0 } ,
  {"mut_wei"     ,   'w' , "FLOAT"  , 0 , "Mutation weight"                             , 0 } ,
  {"ref_seq"     ,   'r' , "STRING" , 0 , "Reference sequence"                          , 0 } ,
  {"const_str"     ,   'c' , "STRING" , 0 , "Constraints"                          , 0 } ,
  {"temp"        ,   't' , "FLOAT" , 0 , "temperature KT"                     , 0 } ,
  {"id_thres"    ,   'i' , "FLOAT" , 0 , "Identity threshold for pseudocount" , 0 } ,
  {"lambda"      ,   'l' , "FLOAT" , 0 , "sequence regularization <TO DO>"    , 0 } ,
  {"alpha"       ,   'a' , "FLOAT" , 0 , "Strructure weight"                  , 0 } ,
  {"phi"       ,   'p' , "FLOAT" , 0 , "DCA weight"                         , 0 } ,
  {"gap"         ,   'g' , "INT"   , 0 , "gap mode 0 = sample with gap      " , 0 } ,
  {"hstep"       ,   'h' , "FLOAT" , 0 , "step size"                          , 0 } ,
  {"out_freq"    ,   'f' , "FLOAT" , 0 , "output frequency"                   , 0 } ,
  {"bias"        ,   'b' , "FILE"  , 0 , "initial values"                     , 0 } ,
  {"struct_mode" ,   's' , "INT"   , 0 , "Structure energy or probability"    , 0 } ,
  {"verbose"     ,   'v' , "INT"   , 0 , "1 if verbose 0 otherwise"           , 0 } ,
  {"mode"        ,   'm' , "INT"   , 0 , "1 MCMC"                             , 0 } ,
  { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char* args[2];
  char* bias;
  float temp;
  char* ref_seq;
  char* const_str;
  long int nb_steps;
  int up_steps;
  float  lambda;
  float  alpha;
  float  phi;
  float  hstep;
  float  out_freq;
  int  verbose;
  int  mode;
  int  struct_mode;
  int  gap_mode;
  float  mut_wei;
  float  id_thres;
};

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'n': arguments->nb_steps = atol(arg); break;
    case 'u': arguments->up_steps = atoi(arg); break;
    case 'w': arguments->mut_wei = atof(arg); break;
    case 't': arguments->temp = atof(arg); break;
    case 'l': arguments->lambda = atof(arg); break;
    case 'p': arguments->phi = atof(arg); break;
    case 'i': arguments->id_thres = atof(arg); break;
    case 'f': arguments->out_freq = atof(arg); break;
    case 'g': arguments->gap_mode = atoi(arg); break;
    case 'h': arguments->hstep = atof(arg); break;
    case 'b': arguments->bias = arg; break;
    case 'a': arguments->alpha = atof(arg); break;
    case 'v': arguments->verbose = atoi(arg); break;
    case 'r': arguments->ref_seq = arg; break;
    case 'c': arguments->const_str = arg; break;
    case 'm': arguments->mode = atoi(arg); break;
    case 's': arguments->struct_mode = atoi(arg); break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 2)
        /* Too many arguments. */
        argp_usage (state);
      arguments->args[state->arg_num] = arg;
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 1)
        /* Not enough arguments. */
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };
// Parse arguments -------------------------------------------------------------

// How many positions are different
float ident_pos(Nucleotide* seq_x, Nucleotide* seq_y, int len_seq) {
  int i;
  float tmp = 0.0;
  float inc = 1.0/(float)len_seq;
  for (i = 0; i < len_seq; i++) {
    if (seq_x[i] == seq_y[i]) {
      tmp += inc;
    }
  }
  return tmp;
}

double update_nrj(double**** bias, Nucleotide* cur_seq, int cpos, Nucleotide old_type, \
                  Nucleotide new_type, int len_seq) {
  // compute the change of nrj
  double old_sum = 0.0, new_sum = 0.0;
  int ppij, old_c;
  for (ppij = 0; ppij < len_seq; ppij++) {
    old_c = cur_seq[ppij];
    if (ppij != cpos) {
      if (ppij < cpos) {
        old_sum += bias[ppij][cpos][old_c][old_type];
        new_sum += bias[ppij][cpos][old_c][new_type];
      } else {
        old_sum += bias[cpos][ppij][old_type][old_c];
        new_sum += bias[cpos][ppij][new_type][old_c];
      }
    }
  }
  return new_sum - old_sum + bias[cpos][cpos][new_type][new_type] - \
    bias[cpos][cpos][old_type][old_type];
}

// compute the total score nrj
double compute_nrj(Nucleotide* seq, int len_seq, double**** bias) {
  double tot_nrj = 0.0;
  int i, j;
  for (i = 0; i < len_seq; i++) {
    for (j = i; j < len_seq; j++) {
      tot_nrj += bias[i][j][seq[i]][seq[j]];
    }
  }
  return tot_nrj;
}


double ens_defect(Nucleotide *seq_nuc, int len_seq, char *target_struct) {
  char *seq_str = convert_seq2char(seq_nuc, len_seq);
  vrna_fold_compound_t *fc =
    vrna_fold_compound(seq_str, NULL, VRNA_OPTION_DEFAULT);
  float ens = vrna_pf(fc, NULL);
  double res = vrna_ensemble_defect(fc, target_struct);

  free(seq_str);
  free_pf_arrays();
  vrna_fold_compound_free(fc);
  return res;
}


// Function to check if two nucleotides are a valid RNA pair
int is_valid_pair(Nucleotide base1, Nucleotide base2) {
  return (base1 == ADE && base2 == URA) ||
    (base1 == URA && base2 == ADE) ||
    (base1 == CYT && base2 == GUA) ||
    (base1 == GUA && base2 == CYT) ||
    (base1 == GUA && base2 == URA) ||
    (base1 == URA && base2 == GUA);
}

double compatibility(Nucleotide *seq_nuc, int len_seq, float** pair_mat) {
  int i, j;
  double tot = 0.0;
  for (i = 0; i <= len_seq; i++) {
    for (j = i+1; j <= len_seq; j++) {
      if (pair_mat[i][j] == 1 && is_valid_pair(seq_nuc[i], seq_nuc[j]))
        tot -= 1;
    }
  }

  return tot;
}


double structure_energy(Nucleotide *seq_nuc, int len_seq, char *target_struct) {
  int i, j;
  char *seq_str = convert_seq2char(seq_nuc, len_seq);
  vrna_fold_compound_t *fc =
    vrna_fold_compound(seq_str, NULL, VRNA_OPTION_DEFAULT);
  float nrj = vrna_eval_structure(fc, target_struct);
  free(seq_str);
  vrna_fold_compound_free(fc);
  return nrj;
}

double probability_mat(Nucleotide* seq_nuc, int len_seq, float** pair_mat) {
  int i, j;
  char* seq_str = convert_seq2char(seq_nuc, len_seq);
  char *propensity = (char *)vrna_alloc(sizeof(char) * (strlen(seq_str) + 1));
  vrna_fold_compound_t *vc = vrna_fold_compound(seq_str, NULL, VRNA_OPTION_DEFAULT);
  vrna_ep_t *ptr, *pair_probabilities = NULL;
  double tot = 0.0, tmp;

  float en = pf_fold(seq_str, NULL);
  FLT_OR_DBL  *bppm = export_bppm();
  for (i = 1; i <= len_seq; i++) {
    for (j = i+1; j <= len_seq; j++) {
      tmp = (bppm[vc->iindx[i] - j] - pair_mat[i-1][j-1]);
      tot += tmp * tmp;
    }
  }
  free_pf_arrays();
  vrna_fold_compound_free(vc);
  free(pair_probabilities);
  free(propensity);
  free(seq_str);

  return (double) sqrt(tot);
}

// Sign function
int sign(double x) {
  if (x > 0.0) {
    return 1;
  } else if (x < 0.0) {
    return -1;
  } else {
    // Use copysign to distinguish between +0 and -0
    return (copysign(1.0, x) > 0.0) ? 1 : -1;
  }
}

void update_parms(Nucleotide* tmp_seq, int len_seq, double**** bias, \
                 double**** exp_prob, double inc) {
  int i, j, ai, aj;
  double tmp, acc;
  double tot = 0.0;

  for (i = 0; i < len_seq; i++) {
    for (j = i + 1; j < len_seq; j++) {
      acc = 0.0;
      for (ai = 0; ai <= GAP; ai++) {
        for (aj = 0; aj <= GAP; aj++) {
          acc += bias[i][j][ai][aj];
        }
      }
      for (ai = 0; ai <= GAP; ai++) {
        for (aj = 0; aj <= GAP; aj++) {
          if (tmp_seq[i] == ai && tmp_seq[j] == aj) {
            tmp = 1.0 - exp_prob[i][j][ai][aj];
          } else
            tmp = -exp_prob[i][j][ai][aj];
          bias[i][j][ai][aj] += inc * (tmp - LAMBDA * sign(bias[i][j][ai][aj]) - acc/(25));
        }
      }
    }
  }

  for (i = 0; i < len_seq; i++) {
    acc = 0.;
    for (ai = 0; ai <= GAP; ai++) {
      acc += bias[i][i][ai][ai];
    }
    for (ai = 0; ai <= GAP; ai++) {
      if (tmp_seq[i] == ai)
        tmp = 1.0 - exp_prob[i][i][ai][ai];
      else
        tmp = -exp_prob[i][i][ai][ai];
      bias[i][i][ai][ai] += inc * (tmp - LAMBDA * sign(bias[i][i][ai][ai]) - acc/5);
    }
  }
}

double compute_struct_cont(Nucleotide* cur_seq, int len_seq, float** pair_mat, char* target_struct) {
  if (STRUCT_MODE == 1)
    return probability_mat(cur_seq, len_seq, pair_mat);
  else if (STRUCT_MODE == 2)
    return ens_defect(cur_seq, len_seq, target_struct);
  else if (STRUCT_MODE == 3)
    return structure_energy(cur_seq, len_seq, target_struct);
  else if (STRUCT_MODE == 4)
    return compatibility(cur_seq, len_seq, pair_mat);
  else
   return 0.0;
}

float update_mut(ref_type, old_type, new_type) {
  if (ref_type == old_type) {
    return 1.0;
  } else if (ref_type == new_type) {
    return -1.0;
  } else {
    return 0;
  }
}

// MCMC in DNA sequence space
// - For each step, perform a mutation in the DNA sequence, then compute the nrj
// - We start from the sequence optimized with diagonal terms only
void mcmc_run(long int nb_steps, Nucleotide** msa, double**** bias,
              double**** exp_prob, int nb_seq, int len_seq, float** pair_mat, struct_el_p struct_data,
              double mean_struct_nrj, double mean_nb_mut, Nucleotide* ref_seq, int* pos_to_mut, int nb_pos_mut) {
  int i, new_seq_id, cur_seq_id, run, mc_run;
  double old_nrj = 0.0;
  double mbeta = -1.0/TEMP;
  double delta_nrj, delta;
  Nucleotide* cur_seq = (Nucleotide*)malloc(sizeof(Nucleotide) * len_seq);
  Nucleotide* tmp_seq = (Nucleotide*)malloc(sizeof(Nucleotide) * len_seq);
  Bool* bias_pos = (Bool*)malloc(sizeof(Bool) * len_seq);

  double beta_l = pow(pow(10.0, -3.0)/HSTEP, 1.0/(double)(nb_steps));

  int step, cpos, old_type, new_type, rpos, ppi, struct_id;
  double struct_nrj, new_struct_nrj, delta_struct_nrj;
  double inc=(double)HSTEP;

  // Start from one sequence of the alignment
  double rand_val;

  // number of mutation from the ref sequence
  double nb_mut, alpha_mut, delta_nb_mut, tmp_mut;

  // main mcmc loop
  for (run = 1; run <= nb_steps; run++) {
    cur_seq_id = rand_int(0, nb_seq);
    for (ppi = 0; ppi < len_seq; ppi++) {
      cur_seq[ppi] = msa[cur_seq_id][ppi];
      tmp_seq[ppi] = msa[cur_seq_id][ppi];
    }
    nb_mut = (double)get_number_mut(cur_seq, ref_seq, len_seq);
    /* nb_mut = 0.; */

    // recompute the nrj
    old_nrj = compute_nrj(cur_seq, len_seq, bias);
    struct_id = rand_int(0, struct_data->nb_struct);
    struct_nrj = compute_struct_cont(cur_seq, len_seq, pair_mat, struct_data->struct_ens[struct_id]);

    for (step = 0; step < UP_FREQ; step++) {
      /* pick a position to mutate */
      cpos = pos_to_mut[rand_int(0, nb_pos_mut)];

      old_type = cur_seq[cpos];

      new_type = rand_int(0, GAP + GAP_MODE);
      while (new_type == old_type)
        new_type = rand_int(0, GAP + GAP_MODE);

      // XXX: compute the change of nrj
      delta_nrj = PHI * update_nrj(bias, cur_seq, cpos, old_type, new_type, len_seq);

      // XXX: compute the change of structure score
      tmp_seq[cpos] = new_type;
      new_struct_nrj = compute_struct_cont(tmp_seq, len_seq, pair_mat, struct_data->struct_ens[struct_id]);
      tmp_seq[cpos] = old_type;
      delta_struct_nrj =  ALPHA * (new_struct_nrj - struct_nrj);

      // compute the mutations
      tmp_mut = update_mut(ref_seq[cpos], old_type, new_type);
      delta_nb_mut = ALPHA_MUT * tmp_mut;
      delta = delta_struct_nrj + delta_nrj + delta_nb_mut;

      // metropolis test
      rand_val = (double)rand_float(0, 1);
      if (delta < 0.0 || exp(mbeta * delta) >= rand_val) {
        cur_seq[cpos] = new_type;
        tmp_seq[cpos] = new_type;
        old_nrj += delta_nrj;
        struct_nrj = new_struct_nrj;
        nb_mut += tmp_mut;
      }

      if (OPT == false && OUT_FREQ > rand_float(0, 1)) {
        printf("%6.1f %6.1f ", step, delta_nrj, struct_nrj);
        print_seq(cur_seq, len_seq);
        printf("\n");
      }
    }

    // chose a few positions to change at the time
    if (OPT==true) {
      inc = beta_l * inc;
      update_parms(cur_seq, len_seq, bias, exp_prob, inc);
      ALPHA -= inc * (mean_struct_nrj - struct_nrj);
      ALPHA_MUT -= inc * (mean_nb_mut - nb_mut);
      fprintf(stderr, "%5.1f %5.1f %5.1f %5.1f\n", old_nrj, struct_nrj, ALPHA, ALPHA_MUT);
    }

  }
}

int main(int argc, char *argv[]) {
  struct arguments arguments;
  // this contains sequences, nb_sequences, length of sequences
  srand(time(0));
  msa_el_p msa_data;
  arguments.temp = 1.0;
  arguments.nb_steps = 1000;
  arguments.up_steps = 1000;
  arguments.lambda = 0.0;
  arguments.phi = 1.0;
  arguments.hstep = 0.1;
  arguments.id_thres = 1.0;
  arguments.bias = "-";
  arguments.const_str = "-";
  arguments.alpha = 0.0;
  arguments.verbose = 0;
  arguments.mode = 0;
  arguments.gap_mode = 1;
  arguments.out_freq = 1.0;
  argp_parse(&argp, argc, argv, 0, 0, &arguments);
  TEMP = (float)arguments.temp;
  HSTEP = (float)arguments.hstep;
  LAMBDA = (float)arguments.lambda;
  PHI = (float)arguments.phi;
  OUT_FREQ = (float)arguments.out_freq;
  UP_FREQ = (int)arguments.up_steps;

  // the structure contribution mode
  // 0 = no structure
  // 1 = SB
  // 2 = Ensemble defect
  STRUCT_MODE = (int)arguments.struct_mode;
  // 0 = you have gaps
  GAP_MODE = (int)arguments.gap_mode;

  if (arguments.verbose == 1)
    VERBOSE = true;
  else
    VERBOSE = false;

  if (arguments.mode == 0)
    OPT = true;
  else
    OPT = false;

  // read the MSA
  msa_data = read_msa(arguments.args[0]);

  // read the target structure
  /* char* struct_str = read_struct_file(arguments.args[1]); */
  /* float** pair_mat = get_contact_mat(struct_str, msa_data->len_seq); */
  struct_el_p struct_data = read_struct_list(arguments.args[1]);
  float** pair_mat = get_contact_mat(struct_data->struct_ens[0], msa_data->len_seq);

  // If we start from another DCA, read it
  double**** bias;
  if (strcmp(arguments.bias, "-") != 0) {
    bias = read_score(msa_data->len_seq, arguments.bias);
    if (arguments.alpha == 0.)
      ALPHA = extract_alpha_value(arguments.bias, "# ALPHA");
    else
      ALPHA = arguments.alpha;
    if (arguments.mut_wei == 0.)
      ALPHA_MUT = extract_alpha_value(arguments.bias, "# ALPHA_MUT");
    else
      ALPHA_MUT = arguments.mut_wei;
  } else
    bias = create_score(msa_data->len_seq);

  // compute the experimental probabilities with pseudocount
  int is, js, i, j;
  double**** exp_prob = create_score(msa_data->len_seq);
  double tmp_wei = 0.0;
  double tot_wei = 0.0;

  double* seq_weight = (double*)malloc(sizeof(double) * msa_data->nb_seq);

  if (arguments.id_thres >= 1.0) {
    for (is = 0; is < msa_data->nb_seq; is++) {
      seq_weight[is] = 1.0;
      tot_wei += 1.0;
    }
  } else {
    for (is = 0; is < msa_data->nb_seq; is++) {
      tmp_wei = 0.0;
      for (js = 0; js < msa_data->nb_seq; js++) {
        if (arguments.id_thres <= ident_pos(msa_data->seq[is], msa_data->seq[js], msa_data->len_seq))
          tmp_wei += 1.0;
      }
      seq_weight[is] = 1.0/tmp_wei;
      tot_wei += seq_weight[is];
    }
  }

  // compute the histogram
  if (OPT == true) {
    for (i = 0; i < msa_data->len_seq; i++) {
      for (is = 0; is < msa_data->nb_seq; is++) {
        exp_prob[i][i][msa_data->seq[is][i]][msa_data->seq[is][i]] += seq_weight[is]/(float)tot_wei;
      }
      for (j = i+1; j < msa_data->len_seq; j++) {
        for (is = 0; is < msa_data->nb_seq; is++) {
          exp_prob[i][j][msa_data->seq[is][i]][msa_data->seq[is][j]] += seq_weight[is]/(float)tot_wei;
        }
      }
    }
  }

  double mean_struct_nrj = 0., tmp_struct_nrj=0, tmp=0.;
  // compute the average structure score
  for (is = 0; is < msa_data->nb_seq; is++) {
    tmp_struct_nrj = 0;
    for (j = 0; j < struct_data->nb_struct; j++) {
      if (tmp_struct_nrj > compute_struct_cont(msa_data->seq[is], msa_data->len_seq, pair_mat, struct_data->struct_ens[j]))
        tmp_struct_nrj = compute_struct_cont(msa_data->seq[is], msa_data->len_seq, pair_mat, struct_data->struct_ens[j]);
    }
    /* mean_struct_nrj += tmp_struct_nrj/struct_data->nb_struct; */
    mean_struct_nrj += tmp_struct_nrj;
  }
  mean_struct_nrj = mean_struct_nrj / tot_wei;

  // read constraints
  int pos_to_mut[1000];
  int nb_pos_mut = 0;
  if (strcmp(arguments.const_str, "-") != 0) {
    char *const_str_copy = strdup(arguments.const_str);
    char *token = strtok(const_str_copy, ",");
    while (token != NULL) {
      pos_to_mut[nb_pos_mut++] = atoi(token); // Convert and store
      token = strtok(NULL, ",");
    }
  } else {
    nb_pos_mut = msa_data->len_seq;
    for (i = 0; i < msa_data->len_seq; i++)
      pos_to_mut[i] = i;
  }

  // compute the average mutation number
  Nucleotide* ref_seq = convert_seq(arguments.ref_seq);

  double mean_nb_mut = 0., tmp_mean = 0.;

  for (is = 0; is < msa_data->nb_seq; is++) {
    mean_nb_mut += (double)get_number_mut(ref_seq, msa_data->seq[is], msa_data->len_seq);
  }
  mean_nb_mut = mean_nb_mut / msa_data->nb_seq;

  printf("# NB_STEPS %d\n", arguments.nb_steps);
  printf("# UP_FREQ %d\n", UP_FREQ);
  printf("# TEMP %f\n", TEMP);
  printf("# LAMBDA %f\n", LAMBDA);
  printf("# ID_THRES %f\n", arguments.id_thres);

  // main function
  mcmc_run(arguments.nb_steps, msa_data->seq, bias, exp_prob, msa_data->nb_seq,
           msa_data->len_seq, pair_mat, struct_data, mean_struct_nrj, mean_nb_mut,
           ref_seq, pos_to_mut, nb_pos_mut);

  if (OPT == true) {
    printf("# ALPHA %f\n", ALPHA/TEMP);
    printf("# ALPHA_MUT %f\n", ALPHA_MUT/TEMP);
    print_score(bias, msa_data->len_seq, 1.0/TEMP, msa_data->seq[0]);
  }

  free_score(bias, msa_data->len_seq);
  free_score(exp_prob, msa_data->len_seq);
  free(seq_weight);
  /* free(struct_str); */

  for (i = 0; i < msa_data->len_seq; i++)
    free(pair_mat[i]);
  free(pair_mat);

  free(ref_seq);
  free_msa(msa_data);
  free_struct_s(struct_data);
  return 0;
}
