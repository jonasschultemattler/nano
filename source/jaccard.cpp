#include <unordered_set>

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

double jaccard_similarity(const std::unordered_set<uint64_t> &kmerset_a,
                          const std::unordered_set<uint64_t> &kmerset_b)
{
    // TODO: implement naive jaccard here
    return (double) set_intersection_size(kmerset_a, kmerset_b)/set_union_size(kmerset_a, kmerset_b);
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

void jaccard_similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n])
{
    for(int i = 0; i < n; i++) {
        matrix[i][i] = 0;
        std::unordered_set<uint64_t> kmerset_i;
        fill_ht(filepaths[i], kmerset_i);
        for(int j = i+1; j < n; j++) {
            std::unordered_set<uint64_t> kmerset_j;
            fill_ht(filepaths[j], kmerset_j);
            matrix[i][j] = matrix[j][i] = jaccard_similarity(kmerset_i, kmerset_j);
        }
    }
}


int main(int argc, char** argv)
{
    double matrix[n][n];

    similarities(files, matrix, jaccard_similarity);
    // jaccard_similarities(files, matrix);
    print_matrix(matrix);

}
