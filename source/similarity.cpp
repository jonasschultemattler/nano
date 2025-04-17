#include <unordered_set>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

// #include <farmhash.h>
// #include <murmurhash.h>


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};



static inline constexpr uint8_t countleadingzeros(const uint64_t x) {
    return __builtin_clzll(x);
}

static inline constexpr uint8_t leading_zeros(const uint64_t hash, const uint8_t precision) {
    const uint64_t mask = (1ULL << precision) - 1u;
    const uint8_t rank = std::countl_zero((hash << precision) | mask) + 1;
    return rank;
}

static inline constexpr uint8_t register_index(const uint64_t hash, const uint8_t precision) {
    return hash >> (64 - precision);
}

static inline constexpr uint64_t hash_1(const uint64_t x) {
    return x;
}

static inline constexpr uint64_t hash_2(const uint64_t x) {
    std::hash<uint64_t> hasher;
    return hasher(x);
}

static inline constexpr uint64_t hash_3(const uint64_t x) {
    __uint128_t result = x;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}


uint64_t jaccard(const std::filesystem::path &filepath, uint8_t const k)
{
    // TODO: implement naive jaccard here

    return 0;
}


uint64_t fracMinHashing(const std::filesystem::path &filepath, uint8_t const k, uint64_t (*hashFunc)(uint64_t))
{
    // TODO: implement FracMinHashing here\

    return 0;
}


int main(int argc, char** argv)
{
    uint8_t const k = 31;
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/celegans.k31.unitigs.fa.ust.fa.gz";
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/ecoli1_k31_ust.fa.gz";

    

}
