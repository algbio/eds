#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <sstream>

struct AlignedSequence {
    std::string seq;

    void align_to(const std::string& reference) {
        // Pad with gaps if needed
        if (seq.size() < reference.size()) {
            seq.append(reference.size() - seq.size(), '-');
        }
    }
};

std::mt19937 rng;

// Read the first sequence from a FASTA file
std::string read_fasta(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, sequence;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue;
        sequence += line;
    }
    return sequence;
}

// Write sequences to a FASTA file
void write_fasta(const std::string& filename, const std::vector<AlignedSequence>& msa) {
    std::ofstream outFile(filename);
    for (size_t i = 0; i < msa.size(); ++i) {
        outFile << ">seq" << i + 1 << "\n";
        const std::string& seq = msa[i].seq;
        for (size_t j = 0; j < seq.size(); j += 60) {
            outFile << seq.substr(j, 60) << "\n";
        }
    }
}

// Generate a random base (optionally excluding one)
char random_base(char exclude = '\0') {
    static const std::string bases = "ACGT";
    char b;
    do {
        b = bases[rng() % 4];
    } while (exclude != '\0' && b == exclude);
    return b;
}

// Mutate a sequence with substitutions, deletions, and insertions, aligning insertions across all
std::pair<AlignedSequence, AlignedSequence> mutate_with_insertions(const AlignedSequence& parent, double rate) {
    std::string child_seq;
    std::string updated_parent_seq;
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (size_t i = 0; i < parent.seq.size(); ++i) {
        double r = dist(rng);

        if (parent.seq[i] == '-') {
            // Preserve gaps
            child_seq += '-';
            updated_parent_seq += '-';
        } else if (r < rate / 3.0) {
            // Insertion: new base in child, gap in parent
            char inserted = random_base();
            child_seq += inserted;
            updated_parent_seq += '-';

            // Then copy the base
            child_seq += parent.seq[i];
            updated_parent_seq += parent.seq[i];
        } else if (r < 2 * rate / 3.0) {
            // Deletion
            child_seq += '-';
            updated_parent_seq += parent.seq[i];
        } else if (r < rate) {
            // Substitution
            child_seq += random_base(parent.seq[i]);
            updated_parent_seq += parent.seq[i];
        } else {
            // No mutation
            child_seq += parent.seq[i];
            updated_parent_seq += parent.seq[i];
        }
    }

    return { {child_seq}, {updated_parent_seq} };
}

// Generate mutations in a binary tree structure
void generate_binary_mutations(std::vector<AlignedSequence>& msa, size_t index, size_t& current, size_t max, double rate) {
    if (current >= max || index >= msa.size()) return;

    // Mutate and update parent
    auto [left, updated_parent] = mutate_with_insertions(msa[index], rate);
    auto [right, _] = mutate_with_insertions(updated_parent, rate);

    // Update parent in MSA
    msa[index] = updated_parent;

    // Align all existing sequences to new parent length
    size_t new_len = updated_parent.seq.size();
    for (auto& s : msa) {
        if (s.seq.size() < new_len) {
            s.seq.append(new_len - s.seq.size(), '-');
        }
    }

    // Align new children to parent
    left.align_to(updated_parent.seq);
    right.align_to(updated_parent.seq);

    msa.push_back(left);
    ++current;
    if (current >= max) return;

    msa.push_back(right);
    ++current;
    if (current >= max) return;

    generate_binary_mutations(msa, msa.size() - 2, current, max, rate);
    generate_binary_mutations(msa, msa.size() - 1, current, max, rate);
}

int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cerr << "Usage: " << argv[0] << " input.fasta output.fasta num_sequences mutation_rate [seed]\n";
        return 1;
    }
    unsigned int seed;
    if (argc >= 6) {
        seed = std::stoul(argv[5]);
        std::cout << "Using provided seed: " << seed << "\n";
    } else {
        seed = std::random_device{}();
        std::cout << "Generated seed: " << seed << "\n";
    }
    std::mt19937 rng(seed);

    std::string input_file = argv[1];
    std::string output_file = argv[2];
    size_t num_sequences = std::stoi(argv[3]);
    double mutation_rate = std::stod(argv[4]);

    std::string root_seq = read_fasta(input_file);
    if (root_seq.empty()) {
        std::cerr << "Error: Failed to read input FASTA.\n";
        return 1;
    }

    std::vector<AlignedSequence> msa;
    msa.push_back({ root_seq });

    size_t current = 1;
    generate_binary_mutations(msa, 0, current, num_sequences, mutation_rate);

    write_fasta(output_file, msa);
    std::cout << "MSA with " << msa.size() << " sequences written to " << output_file << "\n";

    std::cout << "\nSequence lengths and gap counts:\n";
    for (size_t i = 0; i < msa.size(); ++i) {
        size_t gaps = std::count(msa[i].seq.begin(), msa[i].seq.end(), '-');
        std::cout << "seq" << i + 1 << ": " << msa[i].seq.size() << " bp, " << gaps << " gaps\n";
    }

    std::cout << "\nVerifying alignment...\n";
    size_t ref_len = msa[0].seq.size();
    bool same_len = true;
    for (size_t i = 0; i < msa.size(); ++i) {
        size_t len = msa[i].seq.size();
        if (len != ref_len) {
            std::cout << "Mismatch: seq" << i + 1 << " has length " << len << ", expected " << ref_len << "\n";
            same_len = false;
        }
    }
    if (same_len) {
        std::cout << "All sequences same length: YES\n";
    } else {
        std::cout << "All sequences same length: NO\n";
    }

    return 0;
}
