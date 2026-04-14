"""
DNA Transcription and Translation Simulator
==========================================
A modular Python implementation for simulating the central dogma of molecular biology:
DNA → RNA → Protein

Author: Biotech Student
Purpose: Educational tool for understanding transcription and translation processes
"""

class DNATranscriptionTranslation:
    """
    A comprehensive class for simulating DNA transcription and RNA translation.
    
    Biological Background:
    - Transcription: DNA template strand is read 3' → 5' to synthesize RNA 5' → 3'
    - Translation: mRNA is read 5' → 3' in triplet codons to synthesize proteins
    """
    
    def __init__(self):
        """Initialize with the standard genetic code (codon table) and properties."""
        # Standard genetic code: 64 codons → 20 amino acids + 3 stop signals
        self.codon_table = {
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',  # * = stop codon
            'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',  # AUG = start codon
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        # Amino acid properties for advanced analysis
        self.aa_properties = {
            'A': 'Nonpolar', 'R': 'Basic', 'N': 'Polar', 'D': 'Acidic',
            'C': 'Polar', 'Q': 'Polar', 'E': 'Acidic', 'G': 'Nonpolar',
            'H': 'Basic', 'I': 'Nonpolar', 'L': 'Nonpolar', 'K': 'Basic',
            'M': 'Nonpolar', 'F': 'Nonpolar', 'P': 'Nonpolar', 'S': 'Polar',
            'T': 'Polar', 'W': 'Nonpolar', 'Y': 'Polar', 'V': 'Nonpolar'
        }

        # Average molecular weights of amino acids (in Daltons) minus water for peptide bonds
        self.aa_weights = {
            'A': 71.08, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.14,
            'E': 129.12, 'Q': 128.13, 'G': 57.05, 'H': 137.14, 'I': 113.16,
            'L': 113.16, 'K': 128.17, 'M': 131.19, 'F': 147.18, 'P': 97.12,
            'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13
        }

    def validate_dna_sequence(self, dna_seq):
        """Validate DNA sequence contains only valid nucleotides."""
        valid_bases = set('ATCG')
        dna_seq = dna_seq.upper().replace(' ', '')
        
        invalid_bases = set(dna_seq) - valid_bases
        if invalid_bases:
            raise ValueError(f"Invalid DNA bases found: {invalid_bases}")
        
        return True

    def transcribe_dna_to_rna(self, dna_sequence, template_strand=True):
        """Transcribe DNA to RNA following biological transcription rules."""
        self.validate_dna_sequence(dna_sequence)
        dna_sequence = dna_sequence.upper().replace(' ', '')
        
        if template_strand:
            transcription_map = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        else:
            transcription_map = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'}
        
        rna_sequence = ''.join(transcription_map[base] for base in dna_sequence)
        return rna_sequence

    def find_open_reading_frames(self, rna_sequence):
        """Find all possible Open Reading Frames (ORFs) without nesting."""
        rna_sequence = rna_sequence.upper().replace(' ', '')
        orfs = []
        start_codon = 'AUG'
        stop_codons = {'UAA', 'UAG', 'UGA'}
        
        for frame in range(3):
            i = frame
            while i < len(rna_sequence) - 2:
                codon = rna_sequence[i:i+3]
                
                if codon == start_codon:
                    start_pos = i
                    j = i + 3
                    found_stop = False
                    
                    while j < len(rna_sequence) - 2:
                        stop_codon = rna_sequence[j:j+3]
                        if stop_codon in stop_codons:
                            orf_seq = rna_sequence[start_pos:j+3]
                            orfs.append((start_pos, j+2, frame+1, orf_seq))
                            # Advance 'i' past this ORF to prevent nested translation
                            i = j + 3 
                            found_stop = True
                            break
                        j += 3
                    
                    if not found_stop:
                        i += 3
                else:
                    i += 3
                    
        return orfs

    def translate_rna_to_protein(self, rna_sequence, find_orfs=True):
        """Translate RNA sequence to protein using the genetic code."""
        rna_sequence = rna_sequence.upper().replace(' ', '')
        
        if find_orfs:
            orfs = self.find_open_reading_frames(rna_sequence)
            results = {'proteins': [], 'orf_info': orfs}
            
            for start_pos, end_pos, frame, orf_seq in orfs:
                protein = self._translate_sequence(orf_seq)
                results['proteins'].append({
                    'sequence': protein,
                    'start_pos': start_pos,
                    'end_pos': end_pos,
                    'frame': frame,
                    'length': len(protein)
                })
            
            return results
        else:
            protein = self._translate_sequence(rna_sequence)
            return {'proteins': [{'sequence': protein, 'length': len(protein)}]}

    def _translate_sequence(self, rna_sequence):
        """Helper function to translate RNA sequence to amino acid sequence."""
        protein = []
        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i+3]
            if len(codon) != 3:
                break
                
            if codon in self.codon_table:
                amino_acid = self.codon_table[codon]
                if amino_acid == '*':  # Stop codon
                    break
                protein.append(amino_acid)
            else:
                protein.append('X')  # X = unknown amino acid
                
        return ''.join(protein)

    def analyze_protein_properties(self, protein_sequence):
        """Analyze basic properties and molecular weight of translated protein."""
        if not protein_sequence:
            return {'error': 'Empty protein sequence'}
        
        aa_count = {}
        property_count = {'Nonpolar': 0, 'Polar': 0, 'Acidic': 0, 'Basic': 0}
        molecular_weight = 18.02 # Add weight of terminal H and OH groups (water)
        
        for aa in protein_sequence:
            if aa in self.aa_properties:
                aa_count[aa] = aa_count.get(aa, 0) + 1
                prop = self.aa_properties[aa]
                property_count[prop] += 1
                molecular_weight += self.aa_weights.get(aa, 0)
        
        total_aa = len(protein_sequence)
        
        return {
            'length': total_aa,
            'molecular_weight_kDa': round(molecular_weight / 1000, 2), # Convert Da to kDa
            'composition': aa_count,
            'properties': property_count,
            'property_percentages': {
                prop: round(count/total_aa * 100, 1) 
                for prop, count in property_count.items()
            }
        }

    def simulate_transcription_translation(self, dna_sequence, verbose=True):
        """Complete simulation: DNA → RNA → Protein with detailed output."""
        if verbose:
            print("=== DNA TRANSCRIPTION & TRANSLATION SIMULATION ===\n")
            print(f"Input DNA sequence: 5'-{dna_sequence}-3'")
        
        rna_sequence = self.transcribe_dna_to_rna(dna_sequence, template_strand=False)
        if verbose:
            print(f"Transcribed RNA:    5'-{rna_sequence}-3'\n")
        
        translation_results = self.translate_rna_to_protein(rna_sequence, find_orfs=True)
        
        if verbose:
            print("=== TRANSLATION RESULTS ===")
            if translation_results['orf_info']:
                for i, (protein_info, orf_info) in enumerate(zip(
                    translation_results['proteins'], 
                    translation_results['orf_info']
                )):
                    start_pos, end_pos, frame, orf_seq = orf_info
                    print(f"\nORF {i+1} (Frame {frame}, positions {start_pos}-{end_pos}):")
                    print(f"  RNA: {orf_seq}")
                    print(f"  Protein: {protein_info['sequence']}")
                    
                    analysis = self.analyze_protein_properties(protein_info['sequence'])
                    print(f"  Length: {analysis['length']} amino acids")
                    print(f"  Estimated Weight: {analysis['molecular_weight_kDa']} kDa")
                    print(f"  Properties: {analysis['property_percentages']}")
            else:
                print("No complete ORFs found (no start/stop codon pairs)")
        
        results = {
            'dna_input': dna_sequence,
            'rna_sequence': rna_sequence,
            'translation_results': translation_results,
            'summary': {
                'total_orfs': len(translation_results['orf_info']),
                'longest_protein': max(
                    [p['length'] for p in translation_results['proteins']], 
                    default=0
                )
            }
        }
        
        return results

    def simulate_point_mutation(self, dna_sequence, position, new_base):
        """Simulate a point mutation and compare original vs mutated protein."""
        original_result = self.simulate_transcription_translation(dna_sequence, verbose=False)
        
        dna_list = list(dna_sequence.upper())
        if 0 <= position < len(dna_list):
            original_base = dna_list[position]
            dna_list[position] = new_base.upper()
            mutated_dna = ''.join(dna_list)
            
            mutated_result = self.simulate_transcription_translation(mutated_dna, verbose=False)
            
            return {
                'original': {
                    'dna': dna_sequence,
                    'proteins': [p['sequence'] for p in original_result['translation_results']['proteins']]
                },
                'mutated': {
                    'dna': mutated_dna,
                    'proteins': [p['sequence'] for p in mutated_result['translation_results']['proteins']],
                    'mutation': f"{original_base}→{new_base.upper()} at position {position}"
                }
            }
        else:
            raise ValueError(f"Position {position} out of range for sequence length {len(dna_sequence)}")

    def reverse_translate_protein(self, protein_sequence):
        """Generate possible DNA sequences that could code for given protein."""
        reverse_codons = {}
        for codon, aa in self.codon_table.items():
            if aa != '*':
                # Convert RNA codon back to DNA
                reverse_codons.setdefault(aa, []).append(codon.replace('U', 'T'))
        
        if not protein_sequence:
            return []
        
        dna_sequence = []
        for aa in protein_sequence.upper():
            if aa in reverse_codons:
                # Use first available codon
                dna_sequence.append(reverse_codons[aa][0])
            else:
                raise ValueError(f"Unknown amino acid: {aa}")
        
        return [''.join(dna_sequence)]


def demonstrate_central_dogma():
    """Demonstration function showing transcription, translation, and advanced features."""
    simulator = DNATranscriptionTranslation()
    
    print("BIOLOGICAL BACKGROUND:")
    print("=" * 50)
    print("TRANSCRIPTION: DNA template strand (3'→5') is read by RNA polymerase")
    print("               to synthesize complementary RNA (5'→3')")
    print("               Key change: DNA Thymine (T) → RNA Uracil (U)")
    print("\nTRANSLATION:   mRNA (5'→3') is read by ribosomes in triplet codons")
    print("               Each codon specifies one amino acid via genetic code")
    print("               Start: AUG (Methionine), Stop: UAA, UAG, UGA")
    print("\n" + "=" * 50 + "\n")
    
    # Example 1
    print("EXAMPLE 1: Simple gene sequence")
    print("-" * 30)
    dna1 = "ATGGCCTTTGACAAGTAG"
    simulator.simulate_transcription_translation(dna1)
    
    # Example 2
    print("\n" + "=" * 50 + "\n")
    print("EXAMPLE 2: Complex sequence with multiple reading frames")
    print("-" * 50)
    dna2 = "CGATGGCCTTTAGCATGAAATTTGGATAG"
    simulator.simulate_transcription_translation(dna2)

    # Example 3
    print("\n" + "=" * 50 + "\n")
    print("EXAMPLE 3: Mutation Simulation")
    print("-" * 30)
    try:
        # Mutating position 5 (C -> A)
        mutation_result = simulator.simulate_point_mutation("ATGGCCTTTGACAAGTAG", 5, 'A')
        print("Original protein:", mutation_result['original']['proteins'])
        print("Mutated protein: ", mutation_result['mutated']['proteins'])
        print("Mutation:        ", mutation_result['mutated']['mutation'])
    except Exception as e:
        print(f"Mutation simulation error: {e}")
    
    # Example 4
    print("\n" + "=" * 50 + "\n")
    print("EXAMPLE 4: Reverse Translation")
    print("-" * 30)
    try:
        possible_dna = simulator.reverse_translate_protein("MAF")  # Met-Ala-Phe
        print(f"Protein sequence: MAF")
        print(f"Possible DNA:     {possible_dna[0]}")
    except Exception as e:
        print(f"Reverse translation error: {e}")


if __name__ == "__main__":
    demonstrate_central_dogma()
