# set gene dictionary
genes = {}
real_gene = {}
codon_freq = {}
codon_by_amino_acid = {}
# once stored I will have to double-check inside-of each def to make sure that the gene id is real

def read_fasta(fasta_file_path):
    """Read a FASTA file and return a dictionary of gene IDs to sequences."""
    # goes through the file
    current_id = None
    for line in open(fasta_file_path).readlines():
        line = line.rstrip()

        if line.startswith(">"):
            current_id = line[1:]
            # [1:] skips the first line that is a header
            # line 0 is the header
            # ai was used to help with this
        else:
            genes[current_id] = line
    return genes

def read_codon_freqs(path):
    """
    goes through codon file
    :param path:
    :return:
    """
    # used ai to understand how to get through a .csv properly
    codon_by_amino_acid = {}
    for line in open(path).readlines()[1:]:  # read the file skipping the header line
        line = line.rstrip()
        parts = line.split(",")

        # assigning each part of the codon freq
        codon = parts[0]
        amino_acid = parts[1]
        freq = float(parts[2])

        codon_freq[codon] = amino_acid
        codon_by_amino_acid.setdefault(amino_acid, []).append((codon, freq))  # double parentheses for pass a single tuple

    return codon_freq, codon_by_amino_acid


def optimize_gene(gene_id, fasta_file_path, codon_freq_table_file_path):
    """Optimize a gene by replacing codons with the most frequent synonym for each amino acid."""
    genes = read_fasta(fasta_file_path)
    if gene_id not in genes:
        return -1
    # return -1 if there is no id
    codon_freq, codon_by_amino_acid = read_codon_freqs(codon_freq_table_file_path)

    sequence = genes[gene_id]
    codons = []
    i = 0
    while i < len(sequence):
        codon = sequence[i:i+3]
        codons.append(codon)
        i = i+3
    new_codons = []
    for codon in codons:
        amino_acid = codon_freq[codon]
        list_of_codon_freq_pair_options = codon_by_amino_acid[amino_acid]

        best_codon = list_of_codon_freq_pair_options[0][0]
        best_freq = list_of_codon_freq_pair_options[0][1]
        for option in list_of_codon_freq_pair_options:
            this_codon = option[0]
            this_freq = option[1]
            if this_freq > best_freq:
                best_codon = this_codon
                best_freq = this_freq
        new_codons.append(best_codon)
    optimized_codon_sequence = ""
    for j in new_codons:
        optimized_codon_sequence = optimized_codon_sequence + j
    output_file_write = open(gene_id + "_optimized.fasta", "w")
    output_file_write.write(">" + gene_id + "\n")
    output_file_write.write(optimized_codon_sequence + "\n")
    output_file_write.close()

    return optimized_codon_sequence


def deoptimize_gene(gene_id, fasta_file_path, codon_freq_table_file_path):
    """Deoptimize a gene by replacing codons with the least frequent synonym for each amino acid."""
    de_codons = [];
    genes = read_fasta(fasta_file_path)  # I am pretty sure I need to read the files twice to pass the tests
    if gene_id not in genes:
        return -1

    codon_freq, codon_by_amino_acid = read_codon_freqs(codon_freq_table_file_path)
    sequence = genes[gene_id]
    codons = []
    i = 0
    while i < len(sequence):
        codon = sequence[i:i + 3]
        de_codons.append(codon)
        i = i + 3
    new_codons = []
    for codon in de_codons:
        amino_acid = codon_freq[codon]
        list_of_codon_freq_pair_options = codon_by_amino_acid[amino_acid]

        best_codon = list_of_codon_freq_pair_options[0][0]
        best_freq = list_of_codon_freq_pair_options[0][1]
        for option in list_of_codon_freq_pair_options:
            this_codon = option[0]
            this_freq = option[1]
            if this_freq < best_freq:
                best_codon = this_codon
                best_freq = this_freq
        new_codons.append(best_codon)
    deoptimized_codon_sequence = ""
    for j in new_codons:
        deoptimized_codon_sequence = deoptimized_codon_sequence + j
    output_file_write = open(gene_id + "_deoptimized.fasta", "w")
    output_file_write.write(">" + gene_id + "\n")
    output_file_write.write(deoptimized_codon_sequence + "\n")
    output_file_write.close()

    return deoptimized_codon_sequence


def menu():
    """
    A text-based menu that prompts the user to enter file paths to a .fasta file
    and a codon frequency file, a gene id, and whether to optimize or deoptimize
    the gene. Runs in a loop until the user chooses to quit.
    """
    while True:
        fasta_file_path_input = input("Enter fasta file path: ")
        codon_freq_file_path_input = input("Enter codon frequency file path: ")
        gene_id_input = input("Enter gene id: ")
        what_to_do_with_2 = input("Deoptimize or optimize the gene? D or O: ")

        if what_to_do_with_2 == "D" or what_to_do_with_2 == "d":
            result = deoptimize_gene(gene_id_input, fasta_file_path_input, codon_freq_file_path_input)
        elif what_to_do_with_2 == "O" or what_to_do_with_2 == "o":
            result = optimize_gene(gene_id_input, fasta_file_path_input, codon_freq_file_path_input)
        else:
            print("bad input")
            result = None

        if result == -1:
            print("Gene id not found in the given fasta file.")
        elif result is not None:
            print(result)

        again = input("Process another gene? (y/n): ")
        if again == "n" or again == "N":
            break

"""Note - the following lines will call your menu() function if run from this file,
but will not run your menu if the file is imported (e.g., into the tests)"""
if __name__ == "__main__":
    menu()