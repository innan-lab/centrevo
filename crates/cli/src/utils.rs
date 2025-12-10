use centrevo_sim::base::{Nucleotide, Sequence};

pub fn sequence_from_indices(indices: Vec<u8>) -> Sequence {
    let nucleotides: Vec<Nucleotide> = indices
        .into_iter()
        .map(|i| Nucleotide::from_index(i).unwrap_or(Nucleotide::A))
        .collect();
    Sequence::from_nucleotides(nucleotides)
}
