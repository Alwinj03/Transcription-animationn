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
        """Initialize with the standard genetic code (codon table)."""
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

    def validate_dna_sequence(self, dna_seq):
        """
        Validate DNA sequence contains only valid nucleotides.
        
        Args:
            dna_seq (str): DNA sequence to validate
            
        Returns:
            bool: True if valid, False otherwise
            
        Raises:
            ValueError: If sequence contains invalid characters
        """
        valid_bases = set('ATCG')
        dna_seq = dna_seq.upper().replace(' ', '')
        
        invalid_bases = set(dna_seq) - valid_bases
        if invalid_bases:
            raise ValueError(f"Invalid DNA bases found: {invalid_bases}")
        
        return True

    def transcribe_dna_to_rna(self, dna_sequence, template_strand=True):
        """
        Transcribe DNA to RNA following biological transcription rules.
        
        Biological Process:
        - RNA polymerase reads template strand 3' → 5'
        - Synthesizes RNA complementary strand 5' → 3'
        - T in DNA becomes U in RNA
        
        Args:
            dna_sequence (str): Input DNA sequence
            template_strand (bool): If True, treats input as template strand;
                                  if False, treats as coding strand
                                  
        Returns:
            str: Transcribed RNA sequence
            
        Example:
            DNA template: 3'-TACGGCATG-5'
            RNA product:  5'-AUGCCGUAC-3'
        """
        # Validate input
        self.validate_dna_sequence(dna_sequence)
        dna_sequence = dna_sequence.upper().replace(' ', '')
        
        # Transcription mapping
        if template_strand:
            # Template strand: complement and replace T with U
            transcription_map = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        else:
            # Coding strand: just replace T with U (same sequence as RNA)
            transcription_map = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'}
        
        rna_sequence = ''.join(transcription_map[base] for base in dna_sequence)
        
        return rna_sequence

    def find_open_reading_frames(self, rna_sequence):
        """
        Find all possible Open Reading Frames (ORFs) in RNA sequence.
        
        Biological Context:
        - ORF: sequence from start codon (AUG) to stop codon (UAA, UAG, UGA)
        - Reading frames: 3 possible frames depending on where translation starts
        
        Args:
            rna_sequence (str): Input RNA sequence
            
        Returns:
            list: List of tuples (start_pos, end_pos, frame, orf_sequence)
        """
        rna_sequence = rna_sequence.upper().replace(' ', '')
        orfs = []
        start_codon = 'AUG'
        stop_codons = {'UAA', 'UAG', 'UGA'}
        
        # Check all 3 reading frames
        for frame in range(3):
            i = frame
            while i < len(rna_sequence) - 2:
                codon = rna_sequence[i:i+3]
                
                # Found start codon
                if codon == start_codon:
                    start_pos = i
                    # Look for stop codon in same reading frame
                    j = i + 3
                    while j < len(rna_sequence) - 2:
                        stop_codon = rna_sequence[j:j+3]
                        if stop_codon in stop_codons:
                            # Found complete ORF
                            orf_seq = rna_sequence[start_pos:j+3]
                            orfs.append((start_pos, j+2, frame+1, orf_seq))
                            break
                        j += 3
                i += 3
                
        return orfs

    def translate_rna_to_protein(self, rna_sequence, find_orfs=True):
        """
        Translate RNA sequence to protein using the genetic code.
        
        Biological Process:
        - Ribosome reads mRNA 5' → 3' in triplet codons
        - Each codon specifies one amino acid via tRNA
        - Translation starts at AUG (start codon) and ends at stop codons
        
        Args:
            rna_sequence (str): Input RNA sequence
            find_orfs (bool): If True, only translate ORFs; if False, translate entire sequence
            
        Returns:
            dict: Contains 'proteins' list and 'orf_info' if find_orfs=True
        """
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
            # Translate entire sequence from beginning
            protein = self._translate_sequence(rna_sequence)
            return {'proteins': [{'sequence': protein, 'length': len(protein)}]}

    def _translate_sequence(self, rna_sequence):
        """
        Helper function to translate RNA sequence to amino acid sequence.
        
        Args:
            rna_sequence (str): RNA sequence to translate
            
        Returns:
            str: Amino acid sequence (single-letter codes)
        """
        protein = []
        
        # Process sequence in triplets (codons)
        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i+3]
            
            # Skip incomplete codons
            if len(codon) != 3:
                break
                
            # Translate codon to amino acid
            if codon in self.codon_table:
                amino_acid = self.codon_table[codon]
                if amino_acid == '*':  # Stop codon
                    break
                protein.append(amino_acid)
            else:
                # Handle unknown codons (shouldn't happen with valid RNA)
                protein.append('X')  # X = unknown amino acid
        
        return ''.join(protein)

    def analyze_protein_properties(self, protein_sequence):
        """
        Analyze basic properties of translated protein.
        
        Args:
            protein_sequence (str): Amino acid sequence
            
        Returns:
            dict: Protein analysis including composition and properties
        """
        if not protein_sequence:
            return {'error': 'Empty protein sequence'}
        
        # Count amino acids
        aa_count = {}
        property_count = {'Nonpolar': 0, 'Polar': 0, 'Acidic': 0, 'Basic': 0}
        
        for aa in protein_sequence:
            if aa in self.aa_properties:
                aa_count[aa] = aa_count.get(aa, 0) + 1
                prop = self.aa_properties[aa]
                property_count[prop] += 1
        
        total_aa = len(protein_sequence)
        
        return {
            'length': total_aa,
            'composition': aa_count,
            'properties': property_count,
            'property_percentages': {
                prop: round(count/total_aa * 100, 1) 
                for prop, count in property_count.items()
            }
        }

    def simulate_transcription_translation(self, dna_sequence, verbose=True):
        """
        Complete simulation: DNA → RNA → Protein with detailed output.
        
        Args:
            dna_sequence (str): Input DNA sequence
            verbose (bool): If True, print detailed steps
            
        Returns:
            dict: Complete results including RNA, proteins, and analysis
        """
        if verbose:
            print("=== DNA TRANSCRIPTION & TRANSLATION SIMULATION ===\n")
            print(f"Input DNA sequence: 5'-{dna_sequence}-3'")
        
        # Step 1: Transcription
        rna_sequence = self.transcribe_dna_to_rna(dna_sequence, template_strand=False)
        if verbose:
            print(f"Transcribed RNA:    5'-{rna_sequence}-3'\n")
        
        # Step 2: Translation
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
                    
                    # Analyze protein properties
                    analysis = self.analyze_protein_properties(protein_info['sequence'])
                    print(f"  Length: {analysis['length']} amino acids")
                    print(f"  Properties: {analysis['property_percentages']}")
            else:
                print("No complete ORFs found (no start/stop codon pairs)")
        
        # Step 3: Complete analysis
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


def demonstrate_central_dogma():
    """
    Demonstration function showing transcription and translation with examples.
    """
    simulator = DNATranscriptionTranslation()
    
    print("BIOLOGICAL BACKGROUND:")
    print("=" * 50)
    print("TRANSCRIPTION: DNA template strand (3'→5') is read by RNA polymerase")
    print("               to synthesize complementary RNA (5'→3')")
    print("               Key change: DNA Thymine (T) → RNA Uracil (U)")
    print()
    print("TRANSLATION:   mRNA (5'→3') is read by ribosomes in triplet codons")
    print("               Each codon specifies one amino acid via genetic code")
    print("               Start: AUG (Methionine), Stop: UAA, UAG, UGA")
    print("\n" + "=" * 50 + "\n")
    
    # Example 1: Simple sequence with one ORF
    print("EXAMPLE 1: Simple gene sequence")
    print("-" * 30)
    dna1 = "ATGGCCTTTGACAAGTAG"  # Contains start codon ATG and stop codon TAG
    result1 = simulator.simulate_transcription_translation(dna1)
    
    print("\n" + "=" * 50 + "\n")
    
    # Example 2: Longer sequence with multiple potential ORFs
    print("EXAMPLE 2: Complex sequence with multiple reading frames")
    print("-" * 50)
    dna2 = "CGATGGCCTTTAGCATGAAATTTGGATAG"
    result2 = simulator.simulate_transcription_translation(dna2)
    
    return result1, result2


def create_mutation_simulator():
    """
    Extended functionality: Simulate mutations and their effects.
    
    Returns:
        function: Mutation simulation function
    """
    def simulate_point_mutation(dna_sequence, position, new_base):
        """
        Simulate a point mutation and compare original vs mutated protein.
        
        Args:
            dna_sequence (str): Original DNA sequence
            position (int): Position to mutate (0-indexed)
            new_base (str): New base to insert
            
        Returns:
            dict: Comparison of original vs mutated sequences
        """
        simulator = DNATranscriptionTranslation()
        
        # Original sequence
        original_result = simulator.simulate_transcription_translation(dna_sequence, verbose=False)
        
        # Mutated sequence
        dna_list = list(dna_sequence.upper())
        if 0 <= position < len(dna_list):
            original_base = dna_list[position]
            dna_list[position] = new_base.upper()
            mutated_dna = ''.join(dna_list)
            
            mutated_result = simulator.simulate_transcription_translation(mutated_dna, verbose=False)
            
            return {
                'original': {
                    'dna': dna_sequence,
                    'proteins': [p['sequence'] for p in original_result['translation_results']['proteins']]
                },
                'mutated': {
                    'dna': mutated_dna,
                    'proteins': [p['sequence'] for p in mutated_result['translation_results']['proteins']],
                    'mutation': f"{original_base}→{new_base} at position {position}"
                },
                'effect': 'Calculate mutation effects here'
            }
        else:
            raise ValueError(f"Position {position} out of range for sequence length {len(dna_sequence)}")
    
    return simulate_point_mutation


# Advanced feature: Reverse translation (Protein → possible DNA sequences)
def reverse_translate_protein(protein_sequence, codon_usage_bias=None):
    """
    Generate possible DNA sequences that could code for given protein.
    
    Biological Context:
    - Genetic code is degenerate: multiple codons can code for same amino acid
    - Codon usage bias: organisms prefer certain codons over others
    
    Args:
        protein_sequence (str): Amino acid sequence
        codon_usage_bias (dict): Optional codon preferences
        
    Returns:
        list: Possible DNA sequences (showing degeneracy of genetic code)
    """
    # Reverse codon table: amino acid → possible codons
    reverse_codons = {}
    codon_table = DNATranscriptionTranslation().codon_table
    
    for codon, aa in codon_table.items():
        if aa != '*':  # Exclude stop codons
            if aa not in reverse_codons:
                reverse_codons[aa] = []
            # Convert RNA codon back to DNA
            dna_codon = codon.replace('U', 'T')
            reverse_codons[aa].append(dna_codon)
    
    # For demonstration, return first possible sequence
    if not protein_sequence:
        return []
    
    dna_sequence = []
    for aa in protein_sequence.upper():
        if aa in reverse_codons:
            # Use first codon (could implement codon bias here)
            dna_sequence.append(reverse_codons[aa][0])
        else:
            raise ValueError(f"Unknown amino acid: {aa}")
    
    return [''.join(dna_sequence)]


# Main execution and examples
if __name__ == "__main__":
    # Run main demonstration
    result1, result2 = demonstrate_central_dogma()
    
    # Demonstrate mutation simulation
    print("\nEXAMPLE 3: Mutation Simulation")
    print("-" * 30)
    mutation_sim = create_mutation_simulator()
    
    try:
        mutation_result = mutation_sim("ATGGCCTTTGACAAGTAG", 5, 'A')
        print("Original protein:", mutation_result['original']['proteins'])
        print("Mutated protein: ", mutation_result['mutated']['proteins'])
        print("Mutation:", mutation_result['mutated']['mutation'])
    except Exception as e:
        print(f"Mutation simulation error: {e}")
    
    # Demonstrate reverse translation
    print("\nEXAMPLE 4: Reverse Translation")
    print("-" * 30)
    try:
        possible_dna = reverse_translate_protein("MAF")  # Met-Ala-Phe
        print(f"Protein sequence: MAF")
        print(f"Possible DNA:     {possible_dna[0]}")
    except Exception as e:
        print(f"Reverse translation error: {e}")
