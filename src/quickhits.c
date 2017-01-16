/*********

File        quickhits.c
Author      EPW <estebanpw@uma.es>
Description Computes the in-memory hash dictionary of kmers for a query, then searches for hits through the target database and
            produces a histogram of sequences hits

USAGE       <query>            The .fasta file of the query file
            <genomes_database>      The .fasta database containing the genomes
            <sequence_hits_histogram>   The output binary file containing the amount of hits per sequence
            <fixed k of hashtable size> The size of the kmers that will be stored as hash. Consider k=14 takes 2GB, k=16 takes 32GB ram
            <n_threads>        Number of threads to use to build n trees with 4^(n) tree prefixes. Use n=1 for A,C,G,T tree roots; n=2 for AA, AC...TT root trees and 16 threads.



**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include "structs.h"
#include "binaryTreeFunctions.h"
#include "commonFunctions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 64000000 //half a gigabyte if divided by 8 bytes

//Global so its zero-initialized
typedef struct container{
    uint64_t table[4][4][4][4][4][4][4][4][4][4][4][4][4][4];
} Container;

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, long double * minpvalue);
void produce_database(FILE * database, FILE * out_database, long double * pvalues, char * buffer, long double min_p_value, uint64_t n_sequences);

int VERBOSE_ACTIVE = 1;

int main(int argc, char ** av){
    

    clock_t begin, end;
    
    //query to read kmers from, database to find seeds
    FILE * query = NULL, * database = NULL, * out_database = NULL;
    long double minpvalue = 0.5; //Default
    
    init_args(argc, av, &query, &database, &out_database, &minpvalue);
    

    printf("size of container: %"PRIu64"\n", sizeof(Container));
    Container * c_table = (Container *) calloc(1, sizeof(Container));
    if(c_table == NULL) terror("Could not allocate container");
    

    uint64_t n_sequences = 0;

    uint64_t hits;
    unsigned char char_converter[91];
    char_converter['A'] = 0;
    char_converter['C'] = 1;
    char_converter['G'] = 2;
    char_converter['T'] = 3;

    long double p_random_match = powl(ACGT, FIXED_K);
    long double p_rand_by_sumlambdas = p_random_match*(1+SUMOFLAMBDAS);
    long double transform_to_normal = (1- ACGT)/((ACGT)*(ACGT));
    long double mean;
    long double stdev;
    long double pval;
    long double expected_poisson_comp;

    //Variables to account for positions
    //Print info
    fprintf(stdout, "[INFO] Loading query\n");
    //Variables to read kmers
    char c = 'N'; //Char to read character
    //Current length of array and variables for the buffer
    uint64_t temp_len = 0, idx = 0, r = 0, query_len = 0, query_seqs;
    char * temp_seq_buffer = NULL;
    unsigned char curr_kmer[FIXED_K];
    curr_kmer[0] = '\0';
    uint64_t word_size = 0;


    if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }
    //To force reading from the buffer
    idx = READBUF + 1;

    begin = clock();

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
    while((!feof(query) || (feof(query) && idx < r))){

        if(c == '>'){
            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, query); //Skip id

            while(c != '>' && (!feof(query) || (feof(query) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = c;
                    if(word_size < FIXED_K) word_size++;
                    query_len++; 
                    
                }else{
                    if(c != '\n') word_size = 0;
                }
                if(word_size == FIXED_K){
                    //write to hash table
                    
                    c_table->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]]
                    [char_converter[curr_kmer[12]]][char_converter[curr_kmer[13]]]++;

                    memcpy(curr_kmer, &curr_kmer[1], FIXED_K-1);
                    word_size--;
                }
            }
            word_size = 0;
            temp_len++;
            n_sequences++;
            
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);    
        }
        
    }

    end = clock();

    fprintf(stdout, "[INFO] Query completely loaded and of length %"PRIu64". Hash table built. Took %e s\n", temp_len, (double)(end-begin)/CLOCKS_PER_SEC);

    
    begin = clock();

    //Total number of sequences in the DB and total length

    query_seqs = n_sequences;
    uint64_t query_kmers = (query_len) - query_seqs*(FIXED_K - 1);

    n_sequences = 0;


    
    
    //Variables to account for positions
    //Print info
    fprintf(stdout, "[INFO] Loading database and looking up hash table\n");

    //Current length of array
    uint64_t total_kmers = 0, total_kept = 0;
    temp_len = 0; idx = 0; r = 0;
    idx = READBUF + 1;
    curr_kmer[0] = '\0';


    fprintf(out_database, "ID\tHits\tLen\tMU\t\tstdev\t\tPvalue\n");

    hits = 0;
    c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
    while((!feof(database) || (feof(database) && idx < r))){

        if(c == '>'){
            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, database); //Skip id

            while(c != '>' && (!feof(database) || (feof(database) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = c;
                    if(word_size < FIXED_K) word_size++;
                    temp_len++; 
                    
                }else{
                    if(c != '\n') word_size = 0;
                }
                if(word_size == FIXED_K){
                    //ask hash table
                    
                    hits += c_table->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]]
                    [char_converter[curr_kmer[12]]][char_converter[curr_kmer[13]]];

                    memcpy(curr_kmer, &curr_kmer[1], FIXED_K-1);
                    word_size--;
                }
            }
            word_size = 0;
            
            
            
            total_kmers = 2*(temp_len - FIXED_K + 1);
            total_kmers = total_kmers * query_kmers; //Holds the total number of comparisons
            expected_poisson_comp = p_rand_by_sumlambdas*total_kmers;
            mean = expected_poisson_comp/(ACGT);
            stdev = sqrtl(expected_poisson_comp * transform_to_normal);
            pval = 1 - cdf((long double)hits, mean, stdev);
            if(VERBOSE_ACTIVE==1) fprintf(out_database, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%Le\t%Le\t%Le\n", n_sequences, hits, temp_len, mean, stdev, pval);
            if(pval <= minpvalue) total_kept++;

            hits = 0;
            n_sequences++;
            temp_len = 0;
            
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);    
        }
        
    }
    end = clock();
    fprintf(stdout, "[INFO] Database loaded, hits and CDF computed. Took %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fprintf(stdout, "[INFO] The query is %.2f%% similar.\n", (float)100*total_kept/n_sequences);


    



    fprintf(stdout, "[INFO] Done. Deallocating used heap.\n");

    fclose(query);
    fclose(database);
    free(temp_seq_buffer);
    free(c_table);
    
    return 0;
}

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, long double * minpvalue){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           quickhits -query [query] -db [database] -out [output file]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -pval       [Float: 0>=pval>=1] (default: 0.05)\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            exit(1);
        }
        if(strcmp(av[pNum], "-query") == 0){
            *query = fopen64(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open query file");
        }
        if(strcmp(av[pNum], "-db") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open database file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *out_database = fopen64(av[pNum+1], "wt");
            if(out_database==NULL) terror("Could not open output database file");
        }
        if(strcmp(av[pNum], "-pval") == 0){
            *minpvalue = (long double) atof(av[pNum+1]);
            if(*minpvalue < 0 || *minpvalue > 1) terror("Min-pvalue must be between [0,1]");
        }

        pNum++;
    }
    
    if(*query==NULL || *database==NULL || *out_database==NULL) terror("A query, database and output is required");
}

