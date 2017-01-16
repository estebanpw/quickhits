#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include "structs.h"
#include "binaryTreeFunctions.h"
#include "commonFunctions.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))



//Pre-computed optimal k's for 
    //Box, Hunter and Hunter (1978). Statistics for experimenters. Wiley. p. 130.
    // n > 9 (1-p)/p
static uint64_t redir_optimal_k[12] = {13,15,17,18,20,22,23,25,27,28,30,31};
//One based only on database sequence length starting at 10^3
//static uint64_t redir_optimal_k[18] = {4,5,6,8,10,11,13,15,16,18,20,22,23,25,27,28,30,31};

Node_t * getNewLocation(Mempool_t * mp){
    

    Node_t * new_pointer = mp->base + mp->current; //It will increase 32bytes
    mp->current++; //Increase one for the next location
    return new_pointer;
}

inline void init_mem_pool(Mempool_t * mp){
    mp->base = (Node_t *) malloc(POOL_SIZE * sizeof(Node_t));
    if(mp->base == NULL) terror("Could not request memory pool");
    mp->current = 0;
}

Node_t * insertNode_t(uint64_t hash, Mempool_t * mp){

    Node_t * root = getNewLocation(mp);
    if(root == NULL) terror("Could not allocate node");
    root->left = NULL;
    root->right = NULL;
    root->hash = hash;
    //CUCU: Enable if using reps
    //root->reps = 1;    
    return root;
}

int binarySearchNodes(uint64_t query, Node_t * root, Node_t ** found, int insert_yes_no){
    Node_t * current = root;
    

    while(current != NULL){

        if(query == current->hash){
            *found = current;
            //CUCU enable if using reps
            //if(insert_yes_no == 1) current->reps++;
            return 0;  
        }else{
            
            *found=current;
            if(query < current->hash){
                //Go to left
                
                if(current->left != NULL) current = current->left; else return -1;
            }else{
                //Go to right
                
                if(current->right != NULL) current = current->right; else return 1;
            }
        }
    }
    return -2;
}

int binarySearchNodesRanged(Node_t * root, uint64_t r1, uint64_t r2){
    Node_t * current = root;
    

    while(current != NULL){

        if(current->hash >= r1 && current->hash <= r2){
            return 1; //Ranged kmers exist
        }else if(current->hash < r1){  
            current = current->right;
        }else if(current->hash > r2){    
            current = current->left;
        }
    }
    return 0;
}

void * createPrefixTree(void * a){


    BinaryTreeArgs * bta = (BinaryTreeArgs *) a;
    Node_t * insert;

    //Array of mempools to allocate more for each 4 GB needed, up to MAX_MEM_POOLS
    uint64_t requested_pools = 0;
    init_mem_pool(&bta->mem_pools[requested_pools]);

    char c;
    uint64_t curr_pos = 0, crrSeqL = 0, mers_done = 0, threads_redir, curr_hash;
    char b[bta->k_size], b_aux[bta->k_size];
    int result;

    
    while(curr_pos < bta->t_len){
        
        c = bta->sequence[curr_pos];
        curr_pos++;

        if (c == '*') { // Comment, empty or quality (+) line
            crrSeqL = 0; // Reset buffered sequence length
            continue;
        }

        if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
            b[crrSeqL] = c;
            crrSeqL++;
        }else{
            crrSeqL = 0;
        }

        if (crrSeqL >= bta->k_size) { // Full well formed sequence

            threads_redir = hashOfWord(b, bta->prefix_size);
            
            //Insert first node
            if(threads_redir == bta->hash_id){
                curr_hash = hashOfWord(b, bta->k_size);
                if(mers_done == 0){
                    bta->root = insertNode_t(curr_hash, &bta->mem_pools[requested_pools]);    
                    mers_done++;
                }else{
                    //Insert rest of nodes
                    result = binarySearchNodes(curr_hash, bta->root, &insert, 1);
                    if(result == -1){
                        insert->left = insertNode_t(curr_hash, &bta->mem_pools[requested_pools]);
                    }else if(result == 1){
                        insert->right = insertNode_t(curr_hash, &bta->mem_pools[requested_pools]);
                    }

                    if(bta->mem_pools[requested_pools].current == POOL_SIZE){
                        //Malloc new pool
                        requested_pools++;
                        init_mem_pool(&bta->mem_pools[requested_pools]);
                    }
                }
            }
            
	    memcpy(b_aux, b, bta->k_size);
            memcpy(b, &b_aux[1], bta->k_size-1);
            crrSeqL -= 1;
        }

    }

    bta->n_pools = requested_pools;


    //In order traverse to compute accumulated frequencies
    //uint64_t sum = 0;
    //inOrderTraversalPrefixedKmerTree(bta->root, &sum);

    //printf("Sum of this thread: %"PRIu64"\n", sum);


    return NULL;

}


void * createHashTable(void * a){

    BinaryTreeArgs * bta = (BinaryTreeArgs *) a;

    
    char c;
    uint64_t curr_pos = 0, crrSeqL = 0, mers_done = 0, threads_redir, curr_hash;
    char b[bta->k_size], b_aux[bta->k_size];

    
    while(curr_pos < bta->t_len){
        
        c = bta->sequence[curr_pos];
        curr_pos++;

        if (c == '*') { // Comment, empty or quality (+) line
            crrSeqL = 0; // Reset buffered sequence length
            continue;
        }

        if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
            b[crrSeqL] = c;
            crrSeqL++;
        }else{
            crrSeqL = 0;
        }

        if (crrSeqL >= bta->k_size) { // Full well formed sequence

            threads_redir = hashOfWord(b, bta->prefix_size);
            
            //Insert first node
            if(threads_redir == bta->hash_id){
                curr_hash = hashOfWord(b, bta->k_size);
                bta->table[curr_hash] = 1;
                mers_done++;
            }
            
            memcpy(b_aux, b, bta->k_size);
            memcpy(b, &b_aux[1], bta->k_size-1);
            crrSeqL -= 1;
        }

    }



    return NULL;

}

void * searchHashTable(void * a){
    
    BinaryTreeArgs * bta = (BinaryTreeArgs *) a;


    char c;
    uint64_t curr_pos = 0, crrSeqL = 0, threads_redir, curr_hash, currSeq = 0;
    char b[bta->k_size], brev[bta->k_size], b_aux[bta->k_size];
    b[0] = brev[0] = b_aux[0] = '\0';


    uint64_t curr_kmers, optimal_k;
    long double p_random_match, expected_poisson_comp;
    
    //Compute optimal k according to sequence length
    //The number of kmers in the database sequence must be multiplied by two because of reverse strand
    curr_kmers = 2*(bta->seq_lens[0] - (uint64_t) bta->k_size + 1);
    curr_kmers = curr_kmers * bta->query_t_kmers; //Holds the total number of comparisons

    optimal_k = FIXED_K;

    if(bta->hash_id == 0){ //Only one thread should write
        bta->seq_optimal_k[0] = optimal_k;
        p_random_match = powl(ACGT, optimal_k);
        expected_poisson_comp = p_random_match*curr_kmers*(1+ SUMOFLAMBDAS);
        bta->seq_mean[0] = expected_poisson_comp/(ACGT);
        bta->seq_std_dev[0] = sqrtl(expected_poisson_comp * (1- ACGT)/((ACGT)*(ACGT)));
    }

    while(curr_pos < bta->t_len){
        
        c = bta->sequence[curr_pos];
        curr_pos++;


        if (c == '*' && curr_pos < bta->t_len) { // Comment, empty or quality (+) line
            currSeq++;

            //Compute optimal k according to sequence length
            //The number of kmers in the database sequence must be multiplied by two because of reverse strand
            curr_kmers = 2*(bta->seq_lens[currSeq] - (uint64_t) bta->k_size + 1);
            curr_kmers = curr_kmers * bta->query_t_kmers; //Holds the total number of comparisons
           
            if(bta->hash_id == 0){ //Only one thread should write
                expected_poisson_comp = p_random_match*curr_kmers*(1+ SUMOFLAMBDAS);
                bta->seq_optimal_k[currSeq] = optimal_k;
                bta->seq_mean[currSeq] = (long double) expected_poisson_comp/(ACGT);
                bta->seq_std_dev[currSeq] = sqrtl(expected_poisson_comp * (1- ACGT)/((ACGT)*(ACGT)));
            
            }

            crrSeqL = 0; // Reset buffered sequence length
            continue;
        }

        if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
            b[crrSeqL] = c;
            crrSeqL++;
        }else{
            crrSeqL = 0;
        }

        if (crrSeqL >= optimal_k) { // Full well formed sequence
        
            threads_redir = hashOfWord(b, bta->prefix_size);
            strrev(b, brev, optimal_k);
            b[optimal_k] = '\0';
            brev[optimal_k] = '\0';

            
            if(threads_redir == bta->hash_id){
                curr_hash = hashOfWord(b, optimal_k);
                
                if(bta->table[curr_hash] == 1) bta->hits_seq[currSeq]++;
            }

            //Same for reversed
            threads_redir = hashOfWord(brev, bta->prefix_size);
            
            if(threads_redir == bta->hash_id){
                curr_hash = hashOfWord(brev, optimal_k);
                
                if(bta->table[curr_hash] == 1) bta->hits_seq[currSeq]++;

            }
            //This was causing an extremely rare error 
            memcpy(b_aux, b, optimal_k);
            memcpy(b, &b_aux[1], optimal_k-1);
            //memcpy(b, &b[1], optimal_k-1);
            crrSeqL -= 1;
        }

    }
 

    return NULL;

}






////////////

void * searchPrefixedKmerTree(void * a){
    
    BinaryTreeArgs * bta = (BinaryTreeArgs *) a;
    //Node_t * located = NULL;

    char c;
    uint64_t curr_pos = 0, crrSeqL = 0, threads_redir, curr_hash, currSeq = 0, i, ffi_total, lower_hash, upper_hash;
    char b[bta->k_size], brev[bta->k_size], b_upper[bta->k_size], b_lower[bta->k_size], b_aux[bta->k_size];
    b[0] = brev[0] = b_upper[0] = b_lower[0] = '\0';
    //int result;

    uint64_t curr_kmers, optimal_k;
    long double p_random_match, expected_poisson_comp;
    
    //Compute optimal k according to sequence length
    //The number of kmers in the database sequence must be multiplied by two because of reverse strand
    curr_kmers = 2*(bta->seq_lens[0] - (uint64_t) bta->k_size + 1);


    curr_kmers = curr_kmers * bta->query_t_kmers; //Holds the total number of comparisons


    //Remember that 10^9 -> 0....,10^20 -> 32
    optimal_k = redir_optimal_k[MIN(11, MAX(0, lroundl(log10l(curr_kmers)) - 9))];
    //optimal_k = redir_optimal_k[MIN(17, MAX(0, lroundl(log10l(2*bta->seq_lens[0])) - 3))];
    optimal_k = MIN(optimal_k, bta->k_size);
    //optimal_k = bta->k_size; //CUCU

    if(bta->hash_id == 0){ //Only one thread should write
        bta->seq_optimal_k[0] = optimal_k;
        p_random_match = powl(ACGT, optimal_k);
        expected_poisson_comp = p_random_match*curr_kmers*(1+ SUMOFLAMBDAS);
        bta->seq_mean[0] = expected_poisson_comp/(ACGT);
        bta->seq_std_dev[0] = sqrtl(expected_poisson_comp * (1- ACGT)/((ACGT)*(ACGT)));
    }

    while(curr_pos < bta->t_len){
        
        c = bta->sequence[curr_pos];
        curr_pos++;


        if (c == '*' && curr_pos < bta->t_len) { // Comment, empty or quality (+) line
            currSeq++;

            //Compute optimal k according to sequence length
            //The number of kmers in the database sequence must be multiplied by two because of reverse strand
            curr_kmers = 2*(bta->seq_lens[currSeq] - (uint64_t) bta->k_size + 1);
            curr_kmers = curr_kmers * bta->query_t_kmers; //Holds the total number of comparisons

            //Remember that 10^9 -> 0....,10^20 -> 32
            optimal_k = redir_optimal_k[MIN(11, MAX(0, lroundl(log10l(curr_kmers)) - 9))];
            //optimal_k = redir_optimal_k[MIN(17, MAX(0, lroundl(log10l(2*bta->seq_lens[currSeq])) - 3))];
            optimal_k = MIN(optimal_k, bta->k_size);

	       //optimal_k = bta->k_size; // CUCU
            if(bta->hash_id == 0){ //Only one thread should write
                p_random_match = powl(ACGT, optimal_k);
                expected_poisson_comp = p_random_match*curr_kmers*(1+ SUMOFLAMBDAS);
                bta->seq_optimal_k[currSeq] = optimal_k;
                bta->seq_mean[currSeq] = (long double) expected_poisson_comp/(ACGT);
                bta->seq_std_dev[currSeq] = sqrtl(expected_poisson_comp * (1- ACGT)/((ACGT)*(ACGT)));
            
            }
            

            crrSeqL = 0; // Reset buffered sequence length
            continue;
        }

        if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
            b[crrSeqL] = c;
            crrSeqL++;
        }else{
            crrSeqL = 0;
        }
	//if(strncmp("TATACCATACGGATGGATGC", b, optimal_k)==0){ printf("IT EXISTS 1 @ %"PRIu64"\n", curr_pos);}

        if (crrSeqL >= optimal_k) { // Full well formed sequence
	    
	    //if(strncmp("TATACCATACGGATGGATGC", b, optimal_k)==0){ printf("IT EXISTS 2\n");}
            threads_redir = hashOfWord(b, bta->prefix_size);
            strrev(b, brev, optimal_k);
            b[optimal_k] = '\0';
            brev[optimal_k] = '\0';

            
            if(threads_redir == bta->hash_id){
                curr_hash = hashOfWord(b, optimal_k);
                
                if(optimal_k == bta->k_size){ //Standard procedure
                    //result = binarySearchNodes(curr_hash, bta->root, &located, 0);
                    ffi_total = binarySearchNodesRanged(bta->root, curr_hash, curr_hash);

                    bta->hits_seq[currSeq] += ffi_total;
                }else{
                    //We need to use the accumulated frequencies because of the prefix
                    ffi_total = 0;
                    //First, compute the upper and lower bound hashes I.e. (prefixAAAA...A) and (prefixTTTT...T)
                    for(i=0;i<optimal_k;i++) { b_lower[i] = b[i]; b_upper[i] = b[i]; }
                    for(i=optimal_k;i<bta->k_size;i++) { b_lower[i] = 'A'; b_upper[i] = 'T'; }

                    //Second, look for upper bound (prefixTTTT...T)
                    upper_hash = hashOfWord(b_upper, bta->k_size);	    
                    //Thid, look for lower bound (prefixAAAA...A) 
                    lower_hash = hashOfWord(b_lower, bta->k_size);
                    ffi_total = binarySearchNodesRanged(bta->root, lower_hash, upper_hash);
		
                    bta->hits_seq[currSeq] += ffi_total;

		    //if(bta->hash_id == 0) printf("Found %"PRIu64" hits for hash %s and %s\n", ffi_total, b_lower, b_upper);

                }
            }

            //Same for reversed
            threads_redir = hashOfWord(brev, bta->prefix_size);
            
	    if(threads_redir == bta->hash_id){
                curr_hash = hashOfWord(brev, optimal_k);
                
                if(optimal_k == bta->k_size){ //Standard procedure
                    //result = binarySearchNodes(curr_hash, bta->root, &located, 0);
                    ffi_total = binarySearchNodesRanged(bta->root, curr_hash, curr_hash);

                    bta->hits_seq[currSeq] += ffi_total;
                }else{
                    //We need to use the accumulated frequencies because of the prefix
                    ffi_total = 0;
                    //First, compute the upper and lower bound hashes I.e. (prefixAAAA...A) and (prefixTTTT...T)
                    for(i=0;i<optimal_k;i++) { b_lower[i] = brev[i]; b_upper[i] = brev[i]; }
                    for(i=optimal_k;i<bta->k_size;i++) { b_lower[i] = 'A'; b_upper[i] = 'T'; }

                    //Second, look for upper bound (prefixTTTT...T)
                    upper_hash = hashOfWord(b_upper, bta->k_size);

                    lower_hash = hashOfWord(b_lower, bta->k_size);

                    ffi_total = binarySearchNodesRanged(bta->root, lower_hash, upper_hash);

                    bta->hits_seq[currSeq] += ffi_total;


                }
            }
            //This was causing an extremely rare error 
            memcpy(b_aux, b, optimal_k);
            memcpy(b, &b_aux[1], optimal_k-1);
            //memcpy(b, &b[1], optimal_k-1);
            crrSeqL -= 1;
        }

    }
 

    return NULL;

}

void inOrderTraversalPrefixedKmerTree(Node_t * root, uint64_t * sum){
	if(root == NULL){
		return;
	}else{
		inOrderTraversalPrefixedKmerTree(root->left, sum);
		//CUCU: Enable if using node_t with repetitions
		//root->reps = 1;
		//root->acum = *sum;
		//*sum = *sum + root->reps;
		inOrderTraversalPrefixedKmerTree(root->right, sum);
	}
}

/*
typedef struct {
    char * sequence; 
    uint64_t t_len;
    uint64_t n_prefixes;
    uint64_t k_size;
    uint64_t hash_id;
    Mempool_t mem_pools[MAX_MEM_POOLS];
    uint64_t n_pools;
    Node_t * root;
    uint32_t * hits_seq;
} BinaryTreeArgs;
*/

uint64_t computeDiffFreqs(uint64_t result_upper, uint64_t result_lower, Node_t * located_upper, Node_t * located_lower){
    //CUCU: Enable if using reps or acum
    /*
    if(result_upper == 0 && result_lower == 0){
        return located_upper->reps;
    }
    //upper is exact, lower would be inserted to the right
    else if(result_upper == 0 && result_lower == 1){
        return (located_upper->acum + located_upper->reps) - (located_lower->acum + located_lower->reps);
    }
    //upper is exact, lower would be inserted to the left
    else if(result_upper == 0 && result_lower == -1){
        return (located_upper->acum + located_upper->reps) - located_lower->acum;
    }
    //lower is exact, upper would be inserted to the right
    else if(result_lower == 0 && result_upper == 1){
        return (located_upper->acum + located_upper->reps) - located_lower->acum;
    }
    //lower is exact, upper would be inserted to the left
    else if(result_lower == 0 && result_upper == -1){
        return (located_upper->acum) - located_lower->acum;
    }
    //none is exact, upper would go left, lower would go left
    else if(result_upper == -1 && result_lower == -1){
        return (located_upper->acum) - (located_lower->acum);
    }
    //none is exact, upper would go right, lower would go left
    else if(result_upper == 1 && result_lower == -1){
        return (located_upper->acum + located_upper->reps) - (located_lower->acum);
    }
    //none is exact, upper would go left, lower would go right
    else if(result_upper == -1 && result_lower == 1){
        return (located_upper->acum) - (located_lower->acum + located_lower->reps);
    }
    //If they are the same but do not match the query hashes then its a random one
    else if(located_lower->hash == located_upper->hash) return 0;
    //none is exact, upper would go right, lower would go right
    else{
        return (located_upper->acum + located_upper->reps) - (located_lower->acum + located_lower->reps);
       
    }
    */
    return 0;
}

void * computeCDFsByThreads(void * a){
    BinaryTreeArgs * bta = (BinaryTreeArgs *) a;

    uint64_t from, to;
    uint64_t seqs_per_thread = (uint64_t) (floorl((long double) bta->n_sequences_db/ (long double) bta->n_prefixes));


    if(bta->n_prefixes >= bta->n_sequences_db){
        //There are less sequences than threads
        if(bta->hash_id < bta->n_sequences_db){
            bta->pval[bta->hash_id] = 1 - cdf((long double)bta->total_n_hits[bta->hash_id], bta->seq_mean[bta->hash_id], bta->seq_std_dev[bta->hash_id]);
        }
    }else{
        //each one computes a few
        from = bta->hash_id * seqs_per_thread;
        to = (bta->hash_id + 1) * seqs_per_thread;
        uint64_t i;
        for(i=from;i<to;i++){
            bta->pval[i] = 1 - cdf((long double)bta->total_n_hits[i], bta->seq_mean[i], bta->seq_std_dev[i]);
        }
        if(bta->hash_id == bta->n_prefixes-1){
            //If I am the last thread, compute the remaining ones (due to truncating)
            for(i=to;i<bta->n_sequences_db;i++){
                bta->pval[i] = 1 - cdf((long double)bta->total_n_hits[i], bta->seq_mean[i], bta->seq_std_dev[i]);
            }
        }
    }
    return NULL;
}