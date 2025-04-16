#include <unordered_set>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

// #include <farmhash.h>
// #include <murmurhash.h>


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


static inline constexpr uint64_t hash_1(uint64_t kmer) {
    return kmer;
}

static inline constexpr uint64_t hash_2(uint64_t kmer) {
    std::hash<uint64_t> hasher;
    return hasher(kmer);
}

static inline constexpr uint8_t ls1(uint64_t x) {
    return __builtin_clzll(x)-1;
}


uint64_t naive_couting(const std::filesystem::path &filepath, uint8_t const k)
{
    // TODO: implement naive counting here
    std::unordered_set<uint64_t> kmerset;
    uint64_t distinct_kmers = 0;

    // stream over k-mers in DNA file
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            if(kmerset.find(kmer) == kmerset.end()) {
                kmerset.insert(kmer);
                distinct_kmers++;
            }
        }
    }
    return distinct_kmers;
}


uint64_t flajolet_martin(const std::filesystem::path &filepath, uint8_t const k)
{
    // TODO: implement Flajolet-Martin’s algorithm here
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    uint8_t l = 0;

    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hash_2(kmer);
            l = std::max(l, ls1(hash));
        }
    }

    return 1 << l;
}


uint64_t hyperloglog(const std::filesystem::path &filepath, uint8_t const k)
{
    // TODO: implement HyperLogLog here
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    const uint8_t b = 6;
    const uint64_t m = (1<<b);
    const double alpha = 0.709;

    std::vector<uint8_t> M;
    for(int i=0; i < m; i++)
        M.push_back(0);

    uint8_t j = 0;
    uint64_t w;
    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hash_2(kmer);
            w = hash >> (64-b);
            j = hash << b;
            M[j] = std::max(M[j], ls1(w));
        }
    }
    double Z = 0.0;
    for(j=0; j < m; j++) {
        Z += 1.0 / (1 << M[j]);
    }
    return alpha*m*m/Z;
}


int main(int argc, char** argv)
{
    uint8_t const k = 31;
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/example_dna.fa.gz";
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/ecoli1_k31_ust.fa.gz";

    uint64_t kmers = naive_couting(file, k);
    std::cout << "Hashtable: " << kmers << '\n';

    kmers = flajolet_martin(file, k);
    std::cout << "Flajolet Martin: " << kmers << '\n';

    kmers = hyperloglog(file, k);
    std::cout << "HyperLogLog: " << kmers << '\n';

}

