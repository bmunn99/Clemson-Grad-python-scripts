# To extract consensus identity between all species
from Bio import AlignIO

# Load the alignment
alignment = AlignIO.read("/data2/lackey_lab/shalam/PHYL/three/test4/AGO3_short.aln", "clustal")

num_sequences = len(alignment)
num_positions = alignment.get_alignment_length()

consensus = []
overall_identity = 0

# For each position in the alignment
for i in range(num_positions):
    # Get all amino acids at this position
    column = alignment[:, i]
    
    # Count occurrences of each amino acid at this position
    amino_acid_counts = {aa: column.count(aa) for aa in set(column)}

    # Sort amino acids by their occurrences
    sorted_aas = sorted(amino_acid_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Add the most common amino acid to the consensus sequence
    consensus_aa, count = sorted_aas[0]
    consensus.append(consensus_aa)
    
    # Calculate the identity for this position
    identity = count / num_sequences
    overall_identity += identity

# Average the identity scores to get the overall consensus identity
overall_identity /= num_positions

print("Consensus sequence:", ''.join(consensus))
print("Overall consensus identity:", overall_identity)


# To extract consensus identity between Homo sapien and all other species
from Bio import AlignIO

# Load the alignment
alignment = AlignIO.read("/data2/lackey_lab/shalam/PHYL/three/test4/AGO3_short.aln", "clustal")

# Find the index of the Homo species (assuming it's the first sequence)
homo_index = 0

# Initialize a dictionary to store the consensus identity for each species
consensus_identities = {}

# For each species in the alignment
for i in range(len(alignment)):
    if i == homo_index:
        continue  # Skip comparing Homo to itself
    
    species_name = alignment[i].id
    homo_sequence = alignment[homo_index]
    species_sequence = alignment[i]
    
    num_positions = alignment.get_alignment_length()
    matching_positions = sum(1 for h, s in zip(homo_sequence, species_sequence) if h == s)
    
    # Calculate the identity for this species
    identity = matching_positions / num_positions
    
    # Store the identity in the dictionary
    consensus_identities[species_name] = identity

# Print the consensus identities
for species, identity in consensus_identities.items():
    print(f"Consensus identity with {species}: {identity:.2%}")
