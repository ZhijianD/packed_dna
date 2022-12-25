//! A general-purpose genomics crate for dealing with DNA.

#![warn(missing_docs)]
use std::{convert::TryFrom, fmt::Display, str::FromStr, slice::Iter};

// TODO: add a packed module with the PackedDna struct
//
// this struct must have the following:
// 1. A representation that is more memory efficient that simply storing a vector of `Nuc`
// 2. A FromStr implementation (should be case insensitive like the `Nuc` impl)
// 3. A `FromIterator` implementation to construct it from an iterator over `Nuc`s
// 4. A `fn get(&self, idx: usize) -> Nuc` getter for a particular nucleotide
//
// Make sure to unit test and document all elements
// Also, the internal representation of the PackedDna struct should be privately scoped
mod packed {
    use super::*;
    use std::iter::FromIterator;
    use crate::{Nuc, ParseNucError};
    #[derive(Debug, thiserror::Error)]
    #[error("failed to generate DNA from string or iterator given")]
    pub struct ParseDnaError();
    
    #[derive(Debug)]
    pub struct PackedDna {
         size: usize,
         list_of_nuc: Vec<Nuc>,
         
    }
    
    impl FromStr for PackedDna {
        type Err = ParseDnaError;
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let mut v: Vec<Nuc> = Vec::new();
            let mut size: usize = 0;
            for c in s.chars() {
                let temp: Result<Nuc, ParseNucError<char>> = Nuc::try_from(c);
                match temp {
                    Ok(val) => v.push(val), 
                    _ => return Err(ParseDnaError())
                };
                size = size + 1;
            }
            return Ok(PackedDna {size: size, list_of_nuc: v })
        }
    }

    impl FromIterator<Nuc> for PackedDna {
        fn from_iter<I: IntoIterator<Item = Nuc>>(iter: I) -> PackedDna {
            let content: Vec<Nuc> = iter.into_iter().collect();
            let size: usize = content.len();
            PackedDna { size: size, list_of_nuc: content }
        }
    }
    impl PackedDna {
        pub fn get(&self, idx: usize) -> Nuc {
            let e = self.list_of_nuc.get(idx);
            *e.unwrap()
        }
        // count number of each Nuc for given dna
        pub fn nuc_counter(&self) {
            for i in &self.list_of_nuc {
                println!("{:?}", i);
            }
        }
    }
}

/// A nucleotide
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nuc {
    /// Adenine//
    A,
    /// Cytosine
    C,
    /// Guanine
    G,
    /// Thymine
    T,
}

/// An error that can occur when parsing a nucleotide.
#[derive(Debug, thiserror::Error)]
#[error("failed to parse nucleotide from {0}")]
pub struct ParseNucError<T: Display>(T);

impl TryFrom<char> for Nuc {
    type Error = ParseNucError<char>;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_uppercase() {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            _ => Err(ParseNucError(value)),
        }
    }
}

impl FromStr for Nuc {
    type Err = ParseNucError<String>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let upper = s.to_ascii_uppercase();
        match upper.as_str() {
            "A" => Ok(Self::A),
            "C" => Ok(Self::C),
            "G" => Ok(Self::G),
            "T" => Ok(Self::T),
            _ => Err(ParseNucError(upper)),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::iter::FromIterator;

    // use structs and functions to be tested
    use super::*;
    use crate::packed::{PackedDna};

    #[test]
    fn tryfrom_char() {
        // test case-insensitive inputs
        let a = Nuc::try_from('a').unwrap();
        let c = Nuc::try_from('C').unwrap();
        let g = Nuc::try_from('g').unwrap();
        let t = Nuc::try_from('T').unwrap();
        assert!(a == Nuc::A);
        assert!(c == Nuc::C);
        assert!(g == Nuc::G);
        assert!(t == Nuc::T);

        // test invalid char inputs
        let invalid1 = Nuc::try_from('W');
        let invalid2 = Nuc::try_from('?');
        let invalid3 = Nuc::try_from(' ');
        assert!(invalid1.is_err());
        assert!(invalid2.is_err());
        assert!(invalid3.is_err());
    }

    #[test]
    fn fromstr_nuc() {
        // test case insensitive inputs
        let a = Nuc::from_str("a").unwrap();
        let c = Nuc::from_str("C").unwrap();
        let g = Nuc::from_str("g").unwrap();
        let t = Nuc::from_str("T").unwrap();
        assert!(a == Nuc::A);
        assert!(c == Nuc::C);
        assert!(g == Nuc::G);
        assert!(t == Nuc::T);

        // test invalid string inputs
        let invalid1 = Nuc::from_str("W");
        let invalid2 = Nuc::from_str("?");
        let invalid3 = Nuc::from_str(" ");
        let invalid4 = Nuc::from_str("AGC");
        let invalid5 = Nuc::from_str(">?:");
        assert!(invalid1.is_err());
        assert!(invalid2.is_err());
        assert!(invalid3.is_err());
        assert!(invalid4.is_err());
        assert!(invalid5.is_err());
    }

    #[test]
    fn fromstr_dna() {
        // test case insensitive valid dna input
        let a = PackedDna::from_str("a").unwrap();
        let c = PackedDna::from_str("C").unwrap();
        let g = PackedDna::from_str("g").unwrap();
        let t = PackedDna::from_str("T").unwrap();
        let acgt = PackedDna::from_str("aCgT").unwrap();
        assert!(a.get(0) == Nuc::A);
        assert!(c.get(0) == Nuc::C);
        assert!(g.get(0) == Nuc::G);
        assert!(t.get(0) == Nuc::T);
        assert!(acgt.get(0) == Nuc::A);
        assert!(acgt.get(1) == Nuc::C);
        assert!(acgt.get(2) == Nuc::G);
        assert!(acgt.get(3) == Nuc::T);

        // test invalid inputs
        let invalid1 = PackedDna::from_str("AeWrc");
        let invalid2 = PackedDna::from_str(" ACGT_");
        let invalid3 = PackedDna::from_str("+%acgt");
        let invalid4 = PackedDna::from_str("cagt!");
        assert!(invalid1.is_err());
        assert!(invalid2.is_err());
        assert!(invalid3.is_err());
        assert!(invalid4.is_err());
    }

    #[test]
    fn fromiter_dna() {
        let acgt:PackedDna = PackedDna::from_iter(vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T].into_iter());
        assert_eq!(acgt.get(0), Nuc::A);
        assert_eq!(acgt.get(1), Nuc::C);
        assert_eq!(acgt.get(2), Nuc::G);
        assert_eq!(acgt.get(3), Nuc::T);
    }
}
