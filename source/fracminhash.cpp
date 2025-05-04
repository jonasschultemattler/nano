#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include "hashfunctions.hpp"


const std::vector<std::filesystem::path> files = {
        "/Users/adm_js4718fu/datasets/unitigs/ecoli1_k31_ust.fa.gz",
        "/Users/adm_js4718fu/datasets/unitigs/ecoli2_k31_ust.fa.gz",
        "/Users/adm_js4718fu/datasets/unitigs/ecoli4_k31_ust.fa.gz"};
const int n = 3;
const uint8_t k = 31;


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

void print_matrix(double matrix[n][n]) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << '\n';
    }
}


void fill_fracminhashs(const std::filesystem::path &filepath, uint64_t n, uint64_t minhashs[],
    uint64_t a[], uint64_t b[], uint64_t prime, const double s,
    uint64_t (*hashFunc)(uint64_t)=wyhash)
{
    for(int i = 0; i < n; i++)
        minhashs[i] = UINT64_MAX;

    auto fin = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : fin) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hashFunc(kmer);
            for(int i = 0; i < n; i++) {
                uint64_t hash_i = (hash*a[i]+b[i])%prime;
                if(hash_i <= UINT64_MAX*s)
                    minhashs[i] = std::min(hash_i, minhashs[i]);
            }
        }
    }
}


double fracMinHash_similarity(const std::filesystem::path &filepath_a, const std::filesystem::path &filepath_b,
                              const int permutations=100, const double s=0.1)
{
    // TODO: implement MinHashing here
    const uint64_t prime = (1 << 61) - 1;
    uint64_t a[permutations];
    uint64_t b[permutations];
    for(int i = 0; i < permutations; i++) {
        a[i] = (1+std::rand()) % prime;
        b[i] = (std::rand()) % prime;
    }
    uint64_t minhashs_a[permutations];
    uint64_t minhashs_b[permutations];
    fill_fracminhashs(filepath_a, permutations, minhashs_a, a, b, prime, s);
    fill_fracminhashs(filepath_b, permutations, minhashs_b, a, b, prime, s);

    int y = 0;
    for(int i = 0; i < permutations; i++)
        y += minhashs_a[i] == minhashs_b[i];

    return (double) y/permutations;
}


void fracminhash_similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n],
    const int permutations=100, const double s=0.1)
{
    for(int i = 0; i < n; i++) {
        matrix[i][i] = 0;
        for(int j = i+1; j < n; j++) {
            matrix[i][j] = matrix[j][i] = fracMinHash_similarity(filepaths[i], filepaths[j], permutations, s);
        }
    }
}


int main(int argc, char** argv)
{
    double matrix[n][n];
    double s = 0.0000000000000001;
    const int permutations=100;

    if(argc == 1) {
        fracminhash_similarities(files, matrix, permutations, s);
    }
    else if(argc == 2) {
        s = std::stod(argv[1]);
        fracminhash_similarities(files, matrix, permutations, s);
    }
    else {
        std::cout << "usage: optionally provide number of permutations\n";
        return -1;
    }

    print_matrix(matrix);
}
