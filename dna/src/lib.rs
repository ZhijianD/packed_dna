//! A general-purpose genomics crate for dealing with DNA.

use std::{convert::TryFrom, fmt::Display, str::FromStr};

/// A packed module with the PackedDna struct
///
/// this struct have the following:
/// 1. A representation of Nuc with each Nuc taking up 2 bits
/// 2. A FromStr implementation (case insensitive)
/// 3. A `FromIterator` implementation to construct it from an iterator over `Nuc`s
/// 4. A `fn get(&self, idx: usize) -> Nuc` getter for a particular nucleotide
/// 5. A `nuc_count_data(&self)` that returns count of each nucleotide as a vector with 4 element

pub mod packed {
    use super::*;
    use std::{convert::TryInto, iter::FromIterator};

    /// An error that can occur when parsing a DNA.
    #[derive(Debug, thiserror::Error)]
    #[error("failed to generate DNA from string or iterator given")]
    pub struct ParseDnaError();

    /// wrapper of individual byte
    /// content contains 4 Nucleotides, unoccupied bits are zeros
    /// Each nucleotide takes 2 bit with following representation:
    /// {00: Nuc::A,
    ///  01: Nuc::C,
    ///  10: Nuc::G,
    ///  11: Nuc::T}
    /// e.g. ACGT => 00011011 (27)
    #[derive(Debug)]
    struct Node {
        content: u8,
    }

    /// A DNA with its size and list of Nodes
    /// e.g. ACTGAT => list[0b0001_1110, 0b0011_0000], size = 6
    #[derive(Debug)]
    pub struct PackedDna {
        size: usize,
        list_of_nuc: Vec<Node>,
    }

    /// Create Dna from given string slice
    impl FromStr for PackedDna {
        type Err = ParseDnaError;
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let length: usize = s.len();
            let mut num_nodes: usize = length / 4; // num of nodes required
            let last_size: usize = length % 4; // size of last node
            if last_size != 0 {
                // offset
                num_nodes += 1;
            }
            let mut res_list: Vec<Node> = vec![]; // result list
            let mut str_index: usize = 0; // pointer to characters in s
            let mut res_index: usize = 0; // pointer to nodes in result
            while res_index < num_nodes {
                let mut content: u8 = 0b00_00_00_00; // content of current node

                // generate next 4 nucleotide, will throw ParseDnaError if cannot parse nucleotide
                let s0 = Nuc::try_from(s.chars().nth(str_index).unwrap());
                let s0 = match s0 {
                    Ok(val) => val,
                    _ => return Err(ParseDnaError()),
                };
                let s1 = s.chars().nth(str_index + 1);
                let s2 = s.chars().nth(str_index + 2);
                let s3 = s.chars().nth(str_index + 3);

                content = PackedDna::encoder(content, s0, 0);

                if s1 != None {
                    let s1 = Nuc::try_from(s1.unwrap());
                    let s1 = match s1 {
                        Ok(val) => val,
                        _ => return Err(ParseDnaError()),
                    };
                    content = PackedDna::encoder(content, s1, 2);
                }

                if s2 != None {
                    let s2 = Nuc::try_from(s2.unwrap());
                    let s2 = match s2 {
                        Ok(val) => val,
                        _ => return Err(ParseDnaError()),
                    };
                    content = PackedDna::encoder(content, s2, 4);
                }
                if s3 != None {
                    let s3 = Nuc::try_from(s3.unwrap());
                    let s3 = match s3 {
                        Ok(val) => val,
                        _ => return Err(ParseDnaError()),
                    };
                    content = PackedDna::encoder(content, s3, 6);
                }

                res_list.push(Node { content });
                str_index += 4; // iterate over next 4 character in s
                res_index += 1;
            }
            Ok(PackedDna {
                size: length,
                list_of_nuc: res_list,
            })
        }
    }

    /// Create Dna from given Nucleotide Iterator
    impl FromIterator<Nuc> for PackedDna {
        fn from_iter<I: IntoIterator<Item = Nuc>>(iter: I) -> PackedDna {
            let mut iter: <I as IntoIterator>::IntoIter = iter.into_iter();
            let mut size: usize = 0;
            let mut node_list: Vec<Node> = vec![];
            loop {
                let mut content: u8 = 0b00_00_00_00;

                // get next 4 Nucleotides if possible
                let s0: Option<Nuc> = iter.next();
                let s1: Option<Nuc> = if s0 != None { iter.next() } else { None };
                let s2: Option<Nuc> = if s1 != None { iter.next() } else { None };
                let s3: Option<Nuc> = if s2 != None { iter.next() } else { None };

                // encode byte with all 4 nucleotide if they exist
                if s0.is_some() {
                    content = PackedDna::encoder(content, s0.unwrap(), 0);
                    size += 1;
                } else {
                    break;
                }
                if s1.is_some() {
                    content = PackedDna::encoder(content, s1.unwrap(), 2);
                    size += 1;
                }
                if s2.is_some() {
                    content = PackedDna::encoder(content, s2.unwrap(), 4);
                    size += 1;
                }
                if s3.is_some() {
                    content = PackedDna::encoder(content, s3.unwrap(), 6);
                    size += 1;
                }
                node_list.push(Node { content });
                if s1.is_none() || s2.is_none() || s3.is_none() {
                    break;
                }
            }
            PackedDna {
                size,
                list_of_nuc: node_list,
            }
        }
    }

    impl PackedDna {
        // encode content with nucleotide at position, return encoded content
        fn encoder(content: u8, nuc: Nuc, pos: u32) -> u8 {
            match nuc {
                Nuc::C => content | (0b01_00_00_00) >> pos,
                Nuc::G => content | (0b10_00_00_00) >> pos,
                Nuc::T => content | (0b11_00_00_00) >> pos,
                _ => content | (0b00_00_00_00) >> pos,
            }
        }
        /// Get Nucleotide at idx (starting from 0)
        pub fn get(&self, idx: usize) -> Nuc {
            if idx >= self.size {
                panic!("Index out of bound!");
            }
            let idx: u32 = idx.try_into().unwrap(); // type casting
            let node_num: u32 = idx / 4; // which node does this nucleotide belongs to
            let remainder: u32 = idx % 4; // index of nucleotide in this node

            let node = self
                .list_of_nuc
                .iter()
                .nth(node_num.try_into().unwrap())
                .unwrap(); // reference of target node

            // byte representation of the nucleotide
            let res_in_byte = match remainder {
                0 => (node.content & 0b11_00_00_00) >> 6,
                1 => (node.content & 0b00_11_00_00) >> 4,
                2 => (node.content & 0b00_00_11_00) >> 2,
                3 => node.content & 0b00_00_00_11,
                _ => 0b0000_0000,
            };
            match res_in_byte {
                0b00_00_00_00 => Nuc::A,
                0b00_00_00_01 => Nuc::C,
                0b00_00_00_10 => Nuc::G,
                _ => Nuc::T,
            }
        }

        /// count number of each Nuc for given dna, return cnts as vector with 4 elems
        pub fn nuc_count_data(&self) -> Vec<u32> {
            let mut size: usize = self.size;
            let list: &Vec<Node> = &self.list_of_nuc;
            let mut a_cnt: u32 = 0;
            let mut c_cnt: u32 = 0;
            let mut g_cnt: u32 = 0;
            let mut t_cnt: u32 = 0;
            for i in list.iter() {
                let mut off_r: i32 = 6;
                while size > 0 && off_r >= 0 {
                    let nuc: u8 = (i.content & (0b11_00_00_00 >> (6 - off_r))) >> off_r;
                    match nuc {
                        0b00_00_00_00 => a_cnt += 1,
                        0b00_00_00_01 => c_cnt += 1,
                        0b00_00_00_10 => g_cnt += 1,
                        0b00_00_00_11 => t_cnt += 1,
                        _ => a_cnt += 99999,
                    };
                    off_r -= 2;
                    size -= 1;
                }
            }
            vec![a_cnt, c_cnt, g_cnt, t_cnt]
        }

        /// calls nuc_count_data and print out its results
        pub fn nuc_counter(&self) {
            let nucs: Vec<u32> = self.nuc_count_data();
            println!("A: {}", nucs.get(0).unwrap());
            println!("C: {}", nucs.get(1).unwrap());
            println!("G: {}", nucs.get(2).unwrap());
            println!("T: {}", nucs.get(3).unwrap());
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

/// Create nucleotide from char given; raise ParseNucError if invalid parameter is given
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

/// Create nucleotide from string slice given; raise ParseNucError if invalid parameter is given
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
    use crate::packed::PackedDna;

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
    #[should_panic]
    fn fropmstr_empty() {
        // test from_str with empty string, try to get nucleotide from empty dna
        let empty = PackedDna::from_str("").unwrap();
        empty.get(0);
    }

    #[test]
    fn fromiter_dna() {
        let acgt: PackedDna =
            PackedDna::from_iter(vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T].into_iter());
        assert_eq!(acgt.get(0), Nuc::A);
        assert_eq!(acgt.get(1), Nuc::C);
        assert_eq!(acgt.get(2), Nuc::G);
        assert_eq!(acgt.get(3), Nuc::T);
    }

    #[test]
    #[should_panic]
    fn fromiter_empty() {
        // test from_iter with empty string, try to get nucleotide from empty dna
        let acgt: PackedDna = PackedDna::from_iter(vec![].into_iter());
        assert_eq!(acgt.get(0), Nuc::A);
    }

    #[test]
    fn nuc_counter() {
        let dna: PackedDna = PackedDna::from_str("AAACCCCTTTTTGGGGGG").unwrap();
        let data: Vec<u32> = dna.nuc_count_data();
        assert!(*data.get(0).unwrap() == 3);
        assert!(*data.get(1).unwrap() == 4);
        assert!(*data.get(2).unwrap() == 6);
        assert!(*data.get(3).unwrap() == 5);
    }
}
