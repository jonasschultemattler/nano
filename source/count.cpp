#include <unordered_set>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include "hashfunctions.hpp"
#include "perf.hpp"


const uint8_t k = 31;

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


static inline constexpr uint8_t leading_zeros(const uint64_t x) {
    return std::countl_zero(x);
}

static inline constexpr uint8_t leading_zeros(const uint64_t hash, const uint8_t precision) {
    const uint64_t mask = (1ULL << precision) - 1u;
    const uint8_t rank = std::countl_zero((hash << precision) | mask);
    return rank;
}

static inline constexpr uint64_t register_index(const uint64_t hash, const uint8_t precision) {
    return hash >> (64 - precision);
}


uint64_t naive_couting(const std::filesystem::path &filepath)
{
    // TODO: implement naive counting here
    std::unordered_set<uint64_t> kmerset;
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            kmerset.insert(kmer);
        }
    }
    return kmerset.size();
}


uint64_t flajolet_martin(const std::filesystem::path &filepath, uint64_t (*hashFunc)(uint64_t))
{
    // TODO: implement Flajolet-Martin’s algorithm here
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    uint8_t l = 0;

    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hashFunc(kmer);
            uint8_t zeros = leading_zeros(hash);
            l = std::max(l, zeros);
        }
    }

    return 1 << l;
}


uint64_t hyperloglog(const std::filesystem::path &filepath, uint64_t (*hashFunc)(uint64_t))
{
    // TODO: implement HyperLogLog here
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    const uint8_t precision = 8u;
    const uint64_t m = (1<<precision);
    // const double alpha = 0.679;
    const double alpha = 0.7213 / (1.0 + 1.079 / m);

    uint8_t registers[m];
    std::memset(registers, 0, m*sizeof(uint8_t));

    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hashFunc(kmer);
            const uint64_t index = register_index(hash, precision);
            const uint8_t zeros = leading_zeros(hash, precision) + 1;
            registers[index] = std::max(registers[index], zeros);
        }
    }
    double sum = 0.0;
    for(uint64_t j=0; j < m; ++j) {
        sum += 1.0 / (1ULL << registers[j]);
    }

    return alpha*m*m/sum;
}


int main(int argc, char** argv)
{
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/celegans.k31.unitigs.fa.ust.fa.gz";
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/ecoli1_k31_ust.fa.gz";

    PerformanceTracker tracker;

    tracker.start();
    uint64_t count_ht = naive_couting(file);
    std::cout << "Hashtable: " << count_ht << '\n';
    tracker.report();

    tracker.start();
    uint64_t count_fm = flajolet_martin(file, hash_3);
    uint64_t diff;
    if(count_fm > count_ht)
        diff = count_fm - count_ht;
    else
        diff = count_ht - count_fm;
    std::cout << "Flajolet Martin: " << count_fm << " off by " << (double) diff/count_ht*100 << "%\n";
    tracker.report();

    tracker.start();
    uint64_t count_hll = hyperloglog(file, hash_3);
    if(count_hll > count_ht)
        diff = count_hll - count_ht;
    else
        diff = count_ht - count_hll;
    std::cout << "HyperLogLog: " << count_hll << " off by " << (double) diff/count_ht*100 << "%\n";
    tracker.report();

}
