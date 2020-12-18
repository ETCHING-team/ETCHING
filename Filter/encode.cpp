#include "encode.hpp"

string decode(uint64_t encoded, bool rna_flag, size_t kmer_size)
{
  char decoded[32];

  string output;
  int i = sizeof(uint64_t) * BITS_PER_BYTE / 2;

  for (decoded[i--] = '\0'; i >= 0; i--, encoded >>= 2)
    {
      unsigned char byte = encoded & TWO_BITS_MASK;
      if (byte == 2)
        {
          byte = (rna_flag) ? 'U' : 'T';
        }
      else
        {
          byte = 'A' | (byte << 1);
        }
      decoded[i] = byte;
    }
  output = decoded;
  output = output.substr(32-kmer_size);

  return output;
}


string decode(uint64_t encoded, size_t kmer_size)
{
  char decoded[32];

  string output;
  int i = sizeof(uint64_t) * BITS_PER_BYTE / 2;

  for (decoded[i--] = '\0'; i >= 0; i--, encoded >>= 2)
    {
      unsigned char byte = encoded & TWO_BITS_MASK;
      byte = ( byte == 2 ? 'T' : 'A' | (byte << 1) ) ;
      decoded[i] = byte;
    }
  output = decoded;
  output = output.substr(32-kmer_size);

  return output;
}


// Change the 31mer string to 64 bits integer
uint64_t encode(const char *seq)
{
  size_t length = strlen(seq);
  uint64_t result = 0;

  for (size_t i = 0; i < length; i++)
    {
      result = ( result << 2 ) | ( ( seq[i] >> 1 ) & TWO_BITS_MASK );
    }

  return result;
}



void encode_forw(uint64_t & value, const char & base)
{
  value = ( value << 2 ) | ( ( base >> 1 ) & TWO_BITS_MASK );
}



void encode_forw(uint64_t & value, const char & base, const size_t &kl)
{
  value = ( value << 2 ) | ( ( base >> 1 ) & TWO_BITS_MASK );
  value <<= ( 64 - 2 * kl ) ;
  value >>= ( 64 - 2 * kl ) ;
}



uint64_t revcom(unsigned char x)
{
  static const uint64_t table[] = {
    0x02, 0x03, 0x00, 0x01
  };
  return table[x];
}




void encode_reve(uint64_t & value, const char & base, const size_t & kl){
  value =  ( revcom ( ( base >> 1 ) & TWO_BITS_MASK ) << ( 2 * ( kl - 1 ) ) ) | ( value >> 2 ) ;
}



uint64_t select_canonical (uint64_t & value1,uint64_t & value2){
  return value1 > value2 ? value1 : value2 ;
}
