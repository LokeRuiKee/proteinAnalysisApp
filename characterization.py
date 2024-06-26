import numpy as np
import matplotlib.pyplot as plt

def side_chain_charge(aa, pH=7):

    # Using if conditions to know the charge of an input amino acid at a specific pH
    # The numbers that pH is compared against are taken from pKa tables of amino acids. They differ among websites / references.

    if aa == "A" or aa == "G" or aa == "V" or aa == "L" or aa == "I" or aa == "M" or aa == "P" or aa == "F" or aa == "W" or aa == "S" or aa == "T" or aa == "Q" or aa == "N":
        return 0
    elif aa == "C":
        if pH < 8.4:
            return 0
        elif pH == 8.4:
            return -0.5
        else:
            return -1
    elif aa == "Y":
        if pH < 10.5:
            return 0
        elif pH == 10.5:
            return -0.5
        else:
            return -1

    elif aa == "D" or aa == "E":
        if pH < 4:
            return 0
        elif pH == 4:
            return -0.5
        else:
            return -1

    elif aa == "K":
        if pH < 10.5:
            return 1
        elif pH == 10.5:
            return 0.5
        else:
            return 0

    elif aa == "R":
        if pH < 12.5:
            return 1
        elif pH == 12.5:
            return 0.5
        else:
            return 0

    elif aa == "H":
        if pH < 6:
            return 1
        elif pH == 6:
            return 0.5
        else:
            return 0
    else:
        return 0

def adjustTotal(num, protein_seq, pH):

    # This function is to adjust the total charge, which is passed as a function.
    # The total charge that needs to be adjusted is lacking the charge from the C-terminus (last amino acid in the sequence)
    # This function is a helper for the function "total_charge(protein_seq, pH=7)"
    # It checks the last residue and the pH and adjusts the num, which when used will be the total_charge(protein_seq, pH=7)

    if protein_seq[-1] == "G" or protein_seq[-1] == "L" or protein_seq[-1] == "I" or protein_seq[-1] == "D" or protein_seq[-1] == "V":
        if pH > 9.6:
            num = num - 1
        elif pH == 9.6:
            num = num - 0.5
    elif protein_seq[-1] == "A" or protein_seq == "E":
        if pH > 9.7:
            num = num - 1
        elif pH == 9.7:
            num = num - 0.5
    elif protein_seq[-1] == "M" or protein_seq == "S" or protein_seq == "H":
        if pH > 9.2:
            num = num - 1
        elif pH == 9.2:
            num = num - 0.5
    elif protein_seq[-1] == "P":
        if pH > 10.6:
            num = num - 1
        elif pH == 10.6:
            num = num - 0.5
    elif protein_seq[-1] == "F" or protein_seq[-1] == "Q" or protein_seq[-1] == "T" or protein_seq[-1] == "Y":
        if pH > 9.1:
            num = num - 1
        elif pH == 9.1:
            num = num - 0.5
    elif protein_seq[-1] == "W":
        if pH > 9.4:
            num = num - 1
        elif pH == 9.4:
            num = num - 0.5
    elif protein_seq[-1] == "N":
        if pH > 8.8:
            num = num - 1
        elif pH == 8.8:
            num = num - 0.5
    elif protein_seq[-1] == "C":
        if pH > 8.2:
            num = num - 1
        elif pH == 8.2:
            num = num - 0.5
    elif protein_seq[-1] == "K" or protein_seq[-1] == "R":
        if pH > 9:
            num = num - 1
        elif pH == 9:
            num = num - 0.5
    return num

def pair_combinations():

    # This is to create all possible pair amino acid combinations (400 combinations)
    # It puts the 400 pairs in a dictionary where the count is initialized as 0 for each pair

    amino_acids = "ACDEFGHIKLMNPQRSTVWYU"
    combinations = {}
    for i in range(len(amino_acids)):
        for j in range(len(amino_acids)):
            combination = amino_acids[i] + amino_acids[j]
            combinations[combination] = 0
    return combinations

pair_combs = pair_combinations() # This is the variable that stores the pair combinations

def count_duplicate_pair(protein_seq, duplicate):

    # This function is to handle the case where the pair is duplicate --> "AA", "LL", etc. for the remaining 18 residues
    # This is important because a sequence "AAA" has an "AA" count of 2 and not 1 --> based on literature
    # It takes a duplicate (such as AA or YY) and compares it with all possible 2-mers of protein sequence
    # If they are equal, a counter is incremented.
    # The output is a counter reflecting the count of the duplicate in the protein sequence

    sum = 0
    two_mers = len(protein_seq) - 2 + 1
    for i in range(two_mers):
        kmer = protein_seq[i: i+2]
        if kmer == duplicate:
            sum += 1
    return sum

def pair_occurence(protein_seq):

    # This function fills the initialized dictionary above with the real counts of pairs
    # It handles the case if the pair is duplicate using the function above
    # If the pair is not duplicate, it is counted normally using built-in ".count"

    pair_occurence = pair_combs
    for pair in pair_occurence:
        if pair[0] != pair[1]:
            pair_occurence[pair] = protein_seq.count(pair)
        else:
            pair_occurence[pair] = count_duplicate_pair(protein_seq, pair)
    return pair_occurence

def pair_composition(protein_seq):

    # This function is to normalize the above function by the (length of the protein sequence - 1) to find composition
    # Therefore, occurrence is a count, but composition is a percentage
    # It's (length - 1) and not just (length) based on literature, I don't know really why tho.

    composition = {}
    for pair in pair_occurence(protein_seq):
        composition[pair] = round((pair_occurence(protein_seq)[pair] * 100) / (len(protein_seq) - 1), 2)
    return composition

def average_pair_composition(list_of_sequences):

    # This function makes use of the function above it to calculate pair composition for a list of proteins, rather than an individual protein
    # It loops over every protein and saves its composition dictionary
    # In the average dictionary, the pair composition, for each protein sequence, is incremented

    average_pair_composition = pair_combs
    # This is to increment the pairs
    for sequence in list_of_sequences:
        composition = pair_composition(sequence)
        for pair in composition:
            average_pair_composition[pair] += round(composition[pair], 2)

    # This is to find the average value of each pair
    for pair in average_pair_composition:
        average_pair_composition[pair] = round(average_pair_composition[pair] / len(list_of_sequences),2)

    return average_pair_composition

def isoelectric_point(protein_seq):

    # This function uses a brute-force algorithm to locate the specific pH at which the total charge is 0

    min_positive = 100000
    min_negative = -100000
    pH_min_positive = 7
    pH_min_negative = 7
    zero_charge_list = []
    counter = 0
    # it tests all pHs from 2.0 until 12.6, with a step of 0.1, until it reaches a total charge of 0
    for ph in np.arange(2.0, 12.6, 0.1):
        if total_charge(protein_seq, ph) == 0:
            # Upon finding at least one pH that makes the total charge 0, I increment a counter which is initialized as 0
            counter += 1
            # Since many pHs might make the charge 0, I add all the pHs that make the charge 0 in a list called zero_charge_list
            zero_charge_list.append(ph)
        # To handle any errors, I assume no pH will make the total charge 0.
        # So I check if the total charge is > than 0 but < min_positive, in which case the min_positive is updated
        # ... to update the pH at which the charge is the most small positive charge
        elif total_charge(protein_seq, ph) > 0:
            if total_charge(protein_seq, ph) < min_positive:
                min_positive = total_charge(protein_seq, ph)
                pH_min_positive = ph
        # By the same token for cases where the total charge will be negative...
        else:
            if total_charge(protein_seq, ph) > min_negative:
                min_negative = total_charge(protein_seq, ph)
                pH_min_negative = ph

    # Finally, if the counter (initialized as 0) remains 0, meaning no pH made the total charge 0
    # I return the pH --> min positive charge and pH --> min negative charge (least negative)
    # such that the true pH is a number between those two

    if counter:
        return zero_charge_list
    else:
        return [pH_min_positive, pH_min_negative]

def average_isoelectric_point(list_of_sequences):

    # Utilizes the above function to calculate average isoelectric point for a list of sequences

    summ = 0
    for seq in list_of_sequences:
        # Gets an average of the isoelectric points, since the above functions outputs at least two pHs, if not more
        iso = sum(isoelectric_point(seq)) / len(isoelectric_point(seq))
        # This average is added to a counter, initialized as 0
        summ += iso

    # The average of the entire list is then calculated
    return round(summ / len(list_of_sequences), 1)

def total_charge(protein_seq, pH=7):

    # Utilizes side_chain_charge(aa, pH=7) and adjustTotal(num, protein_seq, pH) to calculate the total charge

    total = 0
    for aa in protein_seq:
        total += side_chain_charge(aa, pH)
    if pH < 2:
        total = total + 1
    elif pH == 2:
        total = total + 0.5
    else:
        total = adjustTotal(total, protein_seq, pH)
    return total

def average_total_charge(list_of_sequences, pH = 7):

    # Utilizes the above function to calculate average total charge for a list of protein sequences

    sum = 0
    for seq in list_of_sequences:
        sum += total_charge(seq, pH)

    return round(sum / len(list_of_sequences), 1)

def num_negative_charge(protein_seq):
    return protein_seq.count("D") + protein_seq.count("E")

def average_num_negative_charge(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_negative_charge(seq)

    return round(sum / len(list_of_sequences),0)

def num_positive_charge(protein_seq):
    return protein_seq.count("R") + protein_seq.count("H") + protein_seq.count("K")

def average_num_positive_charge(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_positive_charge(seq)

    return round(sum / len(list_of_sequences),0)

def num_hydrophobic(protein_seq):
    return protein_seq.count("G") + protein_seq.count("A") + protein_seq.count("V") + protein_seq.count("L") + protein_seq.count("I") + protein_seq.count("P") + protein_seq.count("F") + protein_seq.count("M") + protein_seq.count("W")

def average_num_hydrophobic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_hydrophobic(seq)

    return round(sum / len(list_of_sequences),0)

def num_hydrophilic(protein_seq):
    return protein_seq.count("H") + protein_seq.count("E") + protein_seq.count("D") + protein_seq.count("N") + protein_seq.count("Q") + protein_seq.count("K") + protein_seq.count("R")

def average_num_hydrophilic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_hydrophilic(seq)

    return round(sum / len(list_of_sequences),0)

def num_aliphatic(protein_seq):
    return protein_seq.count("G") + protein_seq.count("A") + protein_seq.count("V") + protein_seq.count("L") + protein_seq.count("I") + protein_seq.count("P")

def average_num_aliphatic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_aliphatic(seq)

    return round(sum / len(list_of_sequences),0)

def num_aromatic(protein_seq):
    return protein_seq.count("F") + protein_seq.count("Y") + protein_seq.count("W")

def average_num_aromatic(list_of_sequences):
    sum = 0
    for seq in list_of_sequences:
        sum += num_aromatic(seq)

    return round(sum / len(list_of_sequences),0)

def average_protein_length(list_of_sequences):
    sum = 0

    for seq in list_of_sequences:
        sum += len(seq)

    return round(sum / len(list_of_sequences),0)

def percent_negative_charge(protein_seq):
    return round(num_negative_charge(protein_seq) / len(protein_seq),2)

def average_percent_negative_charge(list_of_sequences):
    return round(average_num_negative_charge(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_positive_charge(protein_seq):
    return round(num_positive_charge(protein_seq) / len(protein_seq),2)

def average_percent_positive_charge(list_of_sequences):
    return round(average_num_positive_charge(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_hydrophobic(protein_seq):
    return round(num_hydrophobic(protein_seq) / len(protein_seq),2)

def average_percent_hydrophobic(list_of_sequences):
    return round(average_num_hydrophobic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_hydrophilic(protein_seq):
    return round(num_hydrophilic(protein_seq) / len(protein_seq),2)

def average_percent_hydrophilic(list_of_sequences):
    return round(average_num_hydrophilic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_aliphatic(protein_seq):
    return round(num_aliphatic(protein_seq) / len(protein_seq),2)

def average_percent_aliphatic(list_of_sequences):
    return round(average_num_aliphatic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def percent_aromatic(protein_seq):
    return round(num_aromatic(protein_seq) / len(protein_seq),2)

def average_percent_aromatic(list_of_sequences):
    return round(average_num_aromatic(list_of_sequences) / average_protein_length(list_of_sequences),2)

def AA_Occurence(protein_seq):

    # Counts the occurence of each amino acid in a protein sequence

    aa_occurence_dic = {}
    list_of_aa = ["A", "G", "P", "V", "L", "I",
                  "M", "C", "F", "Y", "W", "H",
                  "K", "R", "Q", "N", "E", "D",
                  "S", "T", "U"]
    for residue in list_of_aa:
        aa_occurence_dic[residue] = protein_seq.count(residue)
    return aa_occurence_dic

def Avg_Occurence(list_of_sequences):

    # Utilizes the above function to calculate the average occurence of each amino acid in a list of protein sequences

    list_of_aa= ["A", "G", "P", "V", "L", "I",
                 "M", "C", "F", "Y", "W", "H",
                 "K", "R", "Q", "N", "E", "D",
                 "S", "T", "U"]

    average_occurence = {}
    for letter in list_of_aa:
        average_occurence[letter] = 0
    for seq in list_of_sequences:
        aa_occurence = AA_Occurence(seq)
        for aa in aa_occurence:
            average_occurence[aa] += aa_occurence[aa]

    for aa in average_occurence:
        average_occurence[aa] = average_occurence[aa] / len(list_of_sequences)

    return average_occurence

def AA_Composition(protein_seq):

    # Utilizes AA_Occurence to normalize occurence by protein length --> to calculate composition

    aa_composition_dic = AA_Occurence(protein_seq)
    for aa in aa_composition_dic:
        aa_composition_dic[aa] = round((aa_composition_dic[aa]/len(protein_seq))*100,2)
    return aa_composition_dic

def Avg_Composition(list_of_sequences):

    # Utilizes AA_Composition to find the average composition for each amino acid in a list of protein sequences

    list_of_aa = ["A", "G", "P", "V", "L", "I",
                 "M", "C", "F", "Y", "W", "H",
                 "K", "R", "Q", "N", "E", "D",
                 "S", "T", "U"]

    average_composition = {}
    for letter in list_of_aa:
        average_composition[letter] = 0
    for seq in list_of_sequences:
        aa_composition = AA_Composition(seq)
        for aa in aa_composition:
            average_composition[aa] += aa_composition[aa]

    for aa in average_composition:
        average_composition[aa] = round(average_composition[aa] / len(list_of_sequences),2)

    return average_composition

def Protein_Molecular_Weight(protein_seq):

    # Amino acids molecular weight (Dalton) source: https://worldwide.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/#:~:text=The%20average%20molecular%20weight%20of,(kDa)%20is%201%2C000%20daltons
    # Molecular weight can be found by adding the Daltons of each amino acid to a counter --> protein_moelcular_weight
    # The trick is to subtract 18 * (length - 1) at the end because for every two residues to connect, a water molecule is lost
    # That is why there is the 18 * (len(protein_seq) - 1) part...

    aa_mass_dic = {"A": 89, "G": 75, "P": 115, "V": 117, "L": 131, "I": 131,
                  "M": 149, "C": 121, "F": 165, "Y": 181, "W": 204, "H": 155,
                  "K": 146, "R": 174, "Q": 146, "N": 132, "E": 147, "D": 133,
                  "S": 105, "T": 119, "U": 167}
    protein_molecular_weight = 0
    for aa in protein_seq:
        protein_molecular_weight += aa_mass_dic[aa]

    protein_molecular_weight = round(protein_molecular_weight - (18 * (len(protein_seq) - 1)),2)
    return protein_molecular_weight

def Avg_molecular_weight(list_of_sequences):

    # Utilizes the above function to calculate average molecualr weight in a list of protein sequences

    summ = 0
    for seq in list_of_sequences:
        summ += Protein_Molecular_Weight(seq)
    return round(summ / len(list_of_sequences),2)

def get_total_hydrophobicity(protein_seq):

    # Finds the total hydrophobicity of a protein sequence based on Kyte-Doolittle scale

    aa_hydrophobicity = {'I': 4.50, 'V': 4.20, 'L': 3.80, 'F': 2.80, 'C': 2.50, 'M': 1.90, 'A': 1.80, 'G': -0.40, 'T': -0.70, 'S': -0.80,
                         'W': -0.90,'Y': -1.30, 'P': -1.60, 'H': -3.20, 'E': -3.50, 'N': -3.50, 'Q': -3.50, 'D': -3.50, 'K': -3.90, 'R': -4.50, 'U': 2.50}
    total_hydrophobicity = 0
    for aa in protein_seq:
        total_hydrophobicity += aa_hydrophobicity[aa]
    return round(total_hydrophobicity,2)

def get_average_hydrophobicity(list_of_proteins):

    # Utilizes the above function to calculate average hydrophobicity in a list of protein sequences

    sum_hydrophobicity = 0
    for protein in list_of_proteins:
        sum_hydrophobicity += get_total_hydrophobicity(protein)
    avg_hydrophobicity = sum_hydrophobicity/len(list_of_proteins)
    return round(avg_hydrophobicity,2)

def composition_plot(protein, num=1):

    # This function serves to draw a plot for the composition of the protein sequence
    # X axis has every amino acid of the 20 amino acids
    # Y axis has the frequency of every amino acid
    # For example, "A --> 21%, Y --> "13%, and so on"

    plt.title("Protein Composition")
    if num == 1:
        composition = AA_Composition(protein)
    else:
        composition = Avg_Composition(protein)
    plt.xlabel("Amino Acid")
    plt.ylabel("Frequency (%)")
    amino_acid = list(composition.keys())
    frequency = list(composition.values())
    plt.bar(amino_acid, frequency, color ='orange', width = 0.4)
    for i in range(len(amino_acid)):
        plt.text(i, frequency[i], frequency[i], ha = 'center')

    plt.show()

def export(input, name, file_name, num=1):

    # The main function that is used to export all the information about the sequence to a .txt file
    # input is the protein sequence(s)
    # name is the protein name / list name
    # file_name is the name of the file to which info will be exported
    # num is the number of proteins in the input, 1 by default

    exported_file = open(file_name, "w")
    amino_acids = "ACDEFGHIKLMNPQRSTVWYU"

    # For exporting information of an individual protein
    if num == 1:
        exported_file.writelines("Protein name: " + name + "\n")
        exported_file.writelines("Protein sequence: " + input + "\n")
        exported_file.writelines("Protein length: " + str(len(input)) + "\n")
        exported_file.writelines("Protein molecular weight: " + str(Protein_Molecular_Weight(input)) + "\n")
        exported_file.writelines("Protein net charge at pH = 7: " + str(total_charge(input, 7)) + "\n")
        exported_file.writelines("Protein isoelectric point range: " + str(round(isoelectric_point(input)[0],1)) + " - " + str(round(isoelectric_point(input)[-1],1)) + " | Mean: " + str(round(sum(isoelectric_point(input)) / len(isoelectric_point(input)), 2)) + "\n")
        exported_file.writelines("Protein net hydrophobicity: " + str(get_total_hydrophobicity(input)) + "\n")
        exported_file.writelines("Number of positively charged residues: " + str(num_positive_charge(input))+ " and their percentage: " + str(percent_positive_charge(input)) + "\n")
        exported_file.writelines("Number of negatively charged residues: "+ str(num_negative_charge(input)) + " and their percentage: " + str(percent_negative_charge(input)) + "\n")
        exported_file.writelines("Number of hydrophobic residues: "+ str(num_hydrophobic(input)) + " and their percentage: " + str(percent_hydrophobic(input)) + "\n")
        exported_file.writelines("Number of hydrophilic residues: "+ str(num_hydrophilic(input))+ " and their percentage: " + str(percent_hydrophilic(input)) + "\n")
        exported_file.writelines("Number of aliphatic residues: " + str(num_aliphatic(input)) + " and their percentage: " + str(percent_aliphatic(input)) + "\n")
        exported_file.writelines("Number of aromatic residues: " + str(num_aromatic(input)) + " and their percentage: " + str(percent_aromatic(input)) + "\n")

        exported_file.writelines("\n")
        exported_file.writelines("AA" + "\t" + "Occurence" + "\t" + "Composition" + "\n")
        occurence = AA_Occurence(input)
        composition = AA_Composition(input)
        for i in range(20):
            a = amino_acids[i]
            exported_file.writelines(a + "\t" + str(occurence[a]) + "\t" + str(composition[a]) + "\n")
            exported_file.writelines("\n")

    # For exporting information of a list of protein sequences
    else:
        exported_file.writelines("Protein list name: " + name + "\n")
        exported_file.writelines("Dataset protein count: " + str(len(input)) + "\n")
        exported_file.writelines("Average protein length: " + str(average_protein_length(input)) + "\n")
        exported_file.writelines("Average protein molecular weight: " + str(Avg_molecular_weight(input)) + "\n")
        exported_file.writelines("Average protein net charge at pH = 7: " + str(average_total_charge(input, 7)) + "\n")
        exported_file.writelines("Average protein isoelectric point: " + str(average_isoelectric_point(input)) + "\n")
        exported_file.writelines("Average protein hydrophobicity: " + str(get_average_hydrophobicity(input)) + "\n")
        exported_file.writelines("Average number of positively charged residues: " + str(average_num_positive_charge(input))+ " and their percentage: " + str(average_percent_positive_charge(input)) + "\n")
        exported_file.writelines("Average number of negatively charged residues: "+ str(average_num_negative_charge(input)) + " and their percentage: " + str(average_percent_negative_charge(input)) + "\n")
        exported_file.writelines("Average number of hydrophobic residues: "+ str(average_num_hydrophobic(input)) + " and their percentage: " + str(average_percent_hydrophobic(input)) + "\n")
        exported_file.writelines("Average number of hydrophilic residues: "+ str(average_num_hydrophilic(input))+ " and their percentage: " + str(average_percent_hydrophilic(input)) + "\n")
        exported_file.writelines("Average number of aliphatic residues: " + str(average_num_aliphatic(input)) + " and their percentage: " + str(average_percent_aliphatic(input)) + "\n")
        exported_file.writelines("Average number of aromatic residues: " + str(average_num_aromatic(input)) + " and their percentage: " + str(average_percent_aromatic(input)) + "\n")
        exported_file.writelines("\n")
        exported_file.writelines("AA" + "\t" + "Occurence" + "\t" + "Composition" + "\n")
        avg_occurence = Avg_Occurence(input)
        avg_composition = Avg_Composition(input)
        for i in range(20):
            a = amino_acids[i]
            exported_file.writelines(a + "\t" + str(avg_occurence[a]) + "\t" + str(avg_composition[a]) + "\n")
            exported_file.writelines("\n")