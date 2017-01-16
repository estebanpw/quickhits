#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>

#pragma pack(push, 1)

//Structs required for the dotplot workflow
#define MAXLID 200
#define READBUF 1000000000 //500MB
//#define READBUF 5000000 //5MB
#define INITSEQS 20000 //Number of starting sequences (in database)
#define REALLOC_FREQ 40000
#define POINT 4

#define FIXED_K 15
#define UNIFORM 4
#define SUMOFLAMBDAS 0.3333333
#define ACGT 0.25

#define MAXLID 200

//For quickhits memory management
#define POOL_SIZE 31250000 // 31,250,000  structs of 32 bytes :=>: 1 GB ram
//#define POOL_SIZE 15000000 // for small machines
#define MAX_MEM_POOLS 256 // each thread can hold 256 structs each of 1 GB 
//Node-stack realloc policy
#define NODE_T_REAL_POL 4096

//Struct for words program
typedef struct {
    //Each letter is stored using 2 bits
    //We have 4 letters per byte and a
    //maximum of 32 in 'b'
    unsigned char b[8];
} word;

//Struct for words and sort program
typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
    //Using this field to store the sme word in reverse
    //in the same iteration. Either 'f' or 'r'
    char strand;
} wentryR;

//Struct for w2hd program
typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
    //Using this field to store the sme word in reverse
    //in the same iteration. Either 'f' or 'r'
    char strand;
} locationR;

//Struct for w2hd program
typedef struct {
    //Word compressed in binary format
    word w;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint64_t num;
    //The ocurrences with position and
    //sequence
    locationR *locs;
} hashentryR;

//Struct for words and sort program
typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} wentryF;

//Struct for w2hd program
typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} locationF;

//Struct for w2hd program
typedef struct {
    //Word compressed in binary format
    word w;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint64_t num;
    //The ocurrences with position and
    //sequence
    locationF *locs;
} hashentryF;

//Struct for hits, sortHits and filterHits programs
typedef struct {
    //Ocurrence position in sequence X
    uint64_t posX;
    //Ocurrence position in sequence Y
    uint64_t posY;
    //For multiple sequence files this var
    //reflects in what sequence of X file
    //occurs the word
    uint64_t seqX;
    //For multiple sequence files this var
    //reflects in what sequence of Y file
    //occurs the word
    uint64_t seqY;
} hit;

//Struct for FragHits, af2png and leeFrag programs
struct FragFile {
    //Diagonal where the frag is located
    //This value is calculated as:
    //posX - posY
    int64_t diag;
    //Start position in sequence X
    uint64_t xStart;
    //Start position in Sequence Y
    uint64_t yStart;
    //End position in Sequence X
    uint64_t xEnd;
    //End position in Sequence Y
    uint64_t yEnd;
    //Fragment Length
    //For ungaped aligment is:
    //xEnd-xStart+1
    uint64_t length;
    //Number of identities in the
    //fragment
    uint64_t ident;
    //Score of the fragment. This
    //depends on the score matrix
    //used
    uint64_t score;
    //Percentage of similarity. This
    //is calculated as score/scoreMax
    //Where score max is the maximum
    //score possible
    float similarity;
    //sequence number in the 'X' file
    uint64_t seqX;
    //sequence number in the 'Y' file
    uint64_t seqY;
    //synteny block id
    int64_t block;
    //'f' for the forward strain and 'r' for the reverse
    char strand;
};

//Struct for leeSeqDB function
struct Sequence {
    char ident[MAXLID + 1];
    char *datos;
};

//Struct for holding nucleotide frequencies
struct tFreqs{
    uint64_t A;
    uint64_t C;
    uint64_t G;
    uint64_t T;
    uint64_t total;
};

//Struct for calculating karlin and lambda parameters for different query/sequence and PAM matrix
struct statsHSP{
    struct tFreqs tf;
    double karlin;
    double lambda;
    
};

//Struct for each node in the tree to hold hashes
typedef struct node_t {
    uint64_t hash;
    struct node_t * left;
    struct node_t * right;
} Node_t;

//Struct for allocating memory at once
typedef struct mempool_t{
    Node_t * base;
    uint64_t current;
} Mempool_t;

//Struct for reads index tuple
struct rIndex2 {
    char id[MAXLID];
    uint64_t  rNumber;
    uint64_t  rLen;
    uint64_t  rLmasked; //Masked positions
    uint64_t  nonACGT;  //N's
    uint64_t  pos;      //Start position of sequence
    uint64_t  Lac;      //Accumulated length
};


#endif
