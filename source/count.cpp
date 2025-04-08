#include <unordered_set>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


// inline uint32_t hash_1(dna4 &kmer, uint8_t k) {
//     return 0;
// }


uint64_t naive_couting(const std::filesystem::path &filepath, uint8_t const k)
{
    // TODO: implement naive counting here
    std::unordered_set<uint64_t> ht;
    uint64_t kmers = 0;

    // stream over k-mers in DNA file
    auto filestream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    // for (auto & chunk : filestream | seqan3::views::chunk(3000)) {
    for(auto & chunk : filestream) {
        for(auto && kmer : chunk.sequence() | kmer_view) {
            if(ht.count(kmer) == 0) {
                ht.insert(kmer);
                kmers++;
            }
        }
    }
    return kmers;
}


// uint64_t flajolet_martin(const std::filesystem::path &file, uint8_t const k)
// {
//     // TODO: implement Flajolet-Martin’s algorithm here
//     uint64_t kmers = 0;
//     seqan3::sequence_file_input<my_traits> fin{file};

//     for (auto && records : fin | seqan3::views::chunk(3000)) {
//         auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
//         for (auto && kmer : records | kmer_view) {
            
//         }
//     }
//     return kmers;
// }


// uint64_t hyperloglog(const std::filesystem::path &file, uint8_t const k)
// {
//     // TODO: implement HyperLogLog here
//     uint64_t kmers = 0;
//     seqan3::sequence_file_input<my_traits> fin{file};

//     for (auto && records : fin | seqan3::views::chunk(3000)) {
//         auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
//         for (auto && kmer : records | kmer_view) {
            
//         }
//     }
//     return kmers;
// }


int main(int argc, char** argv)
{
    uint8_t const k = 31;
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/example_dna.fa.gz";
    // const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    const std::filesystem::path file = "/Users/adm_js4718fu/datasets/unitigs/ecoli1_k31_ust.fa.gz";

    uint64_t kmers = naive_couting(file, k);
    std::cout << kmers << '\n';

    // kmers = flajolet_martin(&file, k);
    // std::cout << kmers << '\n';

    // kmers = hyperloglog(&file, k);
    // std::cout << kmers << '\n';

}

