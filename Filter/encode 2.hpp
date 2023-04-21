#ifndef ROLLING_ENCODE
#define ROLLING_ENCODE

#include <bitset>
#include <string>
#include <cstring>

#define TWO_BITS_MASK (3)
#define BITS_PER_BYTE (8)
#define BIG_ENOUGH (1024)

using namespace std;

string decode(uint64_t encoded, bool rna_flag, size_t kmer_size);
string decode(uint64_t encoded, size_t kmer_size);
// Change the 31mer string to 64 bits integer
uint64_t encode(const char *seq);
void encode_forw(uint64_t & value, const char & base);
void encode_forw(uint64_t & value, const char & base, const size_t & kl);
void encode_reve(uint64_t & value, const char & base, const size_t & kl);
uint64_t revcom(unsigned char x);
uint64_t select_canonical (uint64_t & value1,uint64_t & value2);

#endif
