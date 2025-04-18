#include <unordered_set>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

// #include <farmhash.h>
// #include <murmurhash.h>


const std::vector<std::filesystem::path> files = {
        "/Users/adm_js4718fu/datasets/unitigs/ecoli1_k31_ust.fa.gz",
        "/Users/adm_js4718fu/datasets/unitigs/ecoli2_k31_ust.fa.gz",
        "/Users/adm_js4718fu/datasets/unitigs/ecoli4_k31_ust.fa.gz"};
const int n = 3;
const uint8_t k = 31;


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


static size_t set_intersection_size(const std::unordered_set<uint64_t> &set_a,
                                    const std::unordered_set<uint64_t> &set_b)
{
    if(set_b.size() < set_a.size()) {
        return set_intersection_size(set_b, set_a);
    }
    size_t cardinality = 0;
    for(std::unordered_set<uint64_t>::const_iterator it = set_a.begin(); it != set_a.end(); it++) {
        if(set_b.find(*it) != set_b.end())
            cardinality++;
    }
    return cardinality;
}

static size_t set_union_size(const std::unordered_set<uint64_t> &set_a,
                             const std::unordered_set<uint64_t> &set_b)
{
    std::unordered_set<uint64_t> set_union;
    set_union.insert(set_a.begin(), set_a.end());
    set_union.insert(set_b.begin(), set_b.end());
    return set_union.size();
}


void fill_ht(const std::filesystem::path &filepath,
             std::unordered_set<uint64_t> &kmerset)
{
    // stream over k-mers in DNA file
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            kmerset.insert(kmer);
        }
    }
}

void print_matrix(double matrix[n][n])
{
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << '\n';
    }
}

void similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n],
                  double (*similarityFunc)(const std::filesystem::path&, const std::filesystem::path&))
{
    for(int i = 0; i < n; i++) {
        matrix[i][i] = 0;
        for(int j = i+1; j < n; j++) {
            matrix[i][j] = matrix[j][i] = similarityFunc(filepaths[i], filepaths[j]);
        }
    }
}


double jaccard_similarity(const std::filesystem::path &filepath_a,
                          const std::filesystem::path &filepath_b)
{
    // TODO: implement naive jaccard here
    std::unordered_set<uint64_t> kmerset_a;
    std::unordered_set<uint64_t> kmerset_b;
    fill_ht(filepath_a, kmerset_a);
    fill_ht(filepath_b, kmerset_b);

    return (double) set_intersection_size(kmerset_a, kmerset_b)/set_union_size(kmerset_a, kmerset_b);
}


uint64_t minHash_similarity(const std::filesystem::path &filepath_a,
                            const std::filesystem::path &filepath_b)
{
    // TODO: implement MinHashing here

    return 0;
}


uint64_t fracMinHash_similarity(const std::filesystem::path &filepath_a,
                                const std::filesystem::path &filepath_b)
{
    // TODO: implement FracMinHashing here

    return 0;
}


int main(int argc, char** argv)
{
    double matrix[n][n];
    similarities(files, matrix, jaccard_similarity);
    print_matrix(matrix);

}
