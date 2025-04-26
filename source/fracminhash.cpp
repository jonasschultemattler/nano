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


double fracMinHash_similarity(const std::filesystem::path &filepath_a,
                              const std::filesystem::path &filepath_b)
{
    // TODO: implement FracMinHashing here

    return 0;
}

void print_matrix(double matrix[n][n]) {
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


int main(int argc, char** argv)
{
    double matrix[n][n];

    similarities(files, matrix, fracMinHash_similarity);
    print_matrix(matrix);
}
