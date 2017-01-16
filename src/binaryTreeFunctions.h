#ifndef BINARY_TREE_FUNCTIONS_H
#define BINARY_TREE_FUNCTIONS_H

typedef struct {
	char * sequence; 		//The sequence loaded in RAM
    uint64_t t_len; 		//The total length of the sequence
    uint64_t n_prefixes;	//The number of prefixes (aka n_threads)
    uint64_t prefix_size;	//The size of the prefix (in letters)
    uint64_t k_size;		//The actual k-size word
    uint64_t hash_id;		//The ID of the thread to know which prefix to use
    Mempool_t mem_pools[MAX_MEM_POOLS];	//The memory pool for the thread
    uint64_t n_pools;		//Number of pools used
    Node_t * root;			//The root node for the tree
    unsigned char * table;  //Hash table
    uint64_t * hits_seq;	//The number of hits per sequence in the database
    uint64_t * seq_lens;	//The length of each sequence in the DB to tell which K size to use
    long double * seq_mean;	//The computed mean for the binomial distribution [WARNING: Only one thread should overwrite]
    long double * seq_std_dev;	//The computed standard deviation for the binomial distribution [WARNING: Only one thread should overwrite]
    uint64_t * seq_optimal_k;	//Optimal k used
    long double * pval;	//The attained pvalue
    uint64_t query_t_kmers;	//Number of k-mers in the query
    uint64_t n_sequences_db; //Number of sequences in the database
    uint64_t * total_n_hits; //A pointer to hits_seq[0] (its a table)

} BinaryTreeArgs;

/**
 * Get a new memory address from the pool mp
 * 
 */
Node_t * getNewLocation(Mempool_t * mp);

/**
 * Initialize the memory pool to later retrieve individual memory addresses
 * 
 */
inline void init_mem_pool(Mempool_t * mp);

/**
 * Insert a node in a binary tree of Node_t
 * 
 */
Node_t * insertNode_t(uint64_t hash, Mempool_t * mp);

/**
 * Binary search over a tree of Node_t nodes. Use insert_yes_no: [1/0] 1 for inserting; 0 for just looking for a query
 * The Node_t found will hold the last visited node
 */
int binarySearchNodes(uint64_t query, Node_t * root, Node_t ** found, int insert_yes_no);

/**
 * Returns 1 if there is at least one node whose hash is contained in [r1, r2]
 * Otherwise returns 0
 */
int binarySearchNodesRanged(Node_t * root, uint64_t r1, uint64_t r2);

/**
 * Creates a binary tree from an array of words (returns the root node)
 * 
 */
void * createPrefixTree(void * a);

/**
 * Creates a hash table from kmers
 * 
 */
void * createHashTable(void * a);

/**
 * Search the hash table
 * 
 */
void * searchHashTable(void * a);

/**
 * Search a whole array of prefixed k-mers over a prefixed tree
 * 
 */
void * searchPrefixedKmerTree(void * a);

/*
	Iterative in-order tree traverse
*/

void inOrderTraversalPrefixedKmerTree(Node_t * root, uint64_t * sum);

/*
	Computes the accumulated number of mers between two prexied kmers in a tree
*/

uint64_t computeDiffFreqs(uint64_t result_upper, uint64_t result_lower, Node_t * located_upper, Node_t * located_lower);

/*
	Compute CDFs (each threads computes n_seqs/n_threads)
*/

void * computeCDFsByThreads(void * a);

#endif /* BINARY_TREE_FUNCTIONS_H */
