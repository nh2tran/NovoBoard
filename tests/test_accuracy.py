"""Tests for accuracy module."""

import pytest
from novoboard.accuracy import parse_raw_sequence
from novoboard import config


class TestParseRawSequence:
    """Test parse_raw_sequence function."""
    
    def test_simple_peptide(self):
        """Test parsing a simple peptide without modifications."""
        ok, peptide = parse_raw_sequence("PEPTIDE")
        assert ok is True
        assert peptide == ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    
    def test_carbamidomethylation(self):
        """Test parsing peptide with carbamidomethylation on cysteine."""
        ok, peptide = parse_raw_sequence("PEPC(+57.02)TIDE")
        assert ok is True
        assert peptide == ['P', 'E', 'P', 'C(Carbamidomethylation)', 'T', 'I', 'D', 'E']
    
    def test_oxidation(self):
        """Test parsing peptide with oxidation on methionine."""
        ok, peptide = parse_raw_sequence("PEPM(+15.99)TIDE")
        assert ok is True
        assert peptide == ['P', 'E', 'P', 'M(Oxidation)', 'T', 'I', 'D', 'E']
    
    def test_deamidation_n(self):
        """Test parsing peptide with deamidation on asparagine.
        
        Note: N(Deamidation) is commented out in the default vocab,
        so this returns False with the parsed peptide.
        """
        ok, peptide = parse_raw_sequence("PEPN(+0.98)TIDE")
        # Deamidation is parsed but not in vocab
        assert peptide == ['P', 'E', 'P', 'N(Deamidation)', 'T', 'I', 'D', 'E']
        assert ok is False  # Not in default vocab
    
    def test_deamidation_q(self):
        """Test parsing peptide with deamidation on glutamine.
        
        Note: Q(Deamidation) is commented out in the default vocab,
        so this returns False with the parsed peptide.
        """
        ok, peptide = parse_raw_sequence("PEPQ(+0.98)TIDE")
        # Deamidation is parsed but not in vocab
        assert peptide == ['P', 'E', 'P', 'Q(Deamidation)', 'T', 'I', 'D', 'E']
        assert ok is False  # Not in default vocab
    
    def test_phosphorylation_s(self):
        """Test parsing peptide with phosphorylation on serine.
        
        Note: S(Phosphorylation) is commented out in the default vocab,
        so this returns False with the parsed peptide.
        """
        ok, peptide = parse_raw_sequence("PEPS(+79.97)TIDE")
        # Phosphorylation is parsed but not in vocab
        assert peptide == ['P', 'E', 'P', 'S(Phosphorylation)', 'T', 'I', 'D', 'E']
        assert ok is False  # Not in default vocab
    
    def test_phosphorylation_t(self):
        """Test parsing peptide with phosphorylation on threonine.
        
        Note: T(Phosphorylation) is commented out in the default vocab,
        so this returns False with the parsed peptide.
        """
        ok, peptide = parse_raw_sequence("PEPT(+79.97)IDE")
        # Phosphorylation is parsed but not in vocab
        assert peptide == ['P', 'E', 'P', 'T(Phosphorylation)', 'I', 'D', 'E']
        assert ok is False  # Not in default vocab
    
    def test_phosphorylation_y(self):
        """Test parsing peptide with phosphorylation on tyrosine.
        
        Note: Y(Phosphorylation) is commented out in the default vocab,
        so this returns False with the parsed peptide.
        """
        ok, peptide = parse_raw_sequence("PEPY(+79.97)IDE")
        # Phosphorylation is parsed but not in vocab
        assert peptide == ['P', 'E', 'P', 'Y(Phosphorylation)', 'I', 'D', 'E']
        assert ok is False  # Not in default vocab
    
    def test_unknown_modification(self):
        """Test parsing peptide with unknown modification returns False."""
        ok, peptide = parse_raw_sequence("PEPX(+99.99)TIDE")
        assert ok is False
    
    def test_multiple_modifications(self):
        """Test parsing peptide with multiple modifications."""
        ok, peptide = parse_raw_sequence("C(+57.02)PEPM(+15.99)TIDE")
        assert ok is True
        assert peptide == ['C(Carbamidomethylation)', 'P', 'E', 'P', 'M(Oxidation)', 'T', 'I', 'D', 'E']
    
    def test_empty_sequence(self):
        """Test parsing empty sequence."""
        ok, peptide = parse_raw_sequence("")
        assert ok is True
        assert peptide == []


class TestConfig:
    """Test config module values."""
    
    def test_vocab_size(self):
        """Test vocabulary size is correct."""
        assert config.vocab_size == 24
    
    def test_vocab_contains_amino_acids(self):
        """Test vocabulary contains standard amino acids."""
        standard_aa = ['A', 'R', 'N', 'D', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        for aa in standard_aa:
            assert aa in config.vocab
    
    def test_mass_aa_values(self):
        """Test some amino acid mass values are correct."""
        assert abs(config.mass_AA['A'] - 71.03711) < 0.0001
        assert abs(config.mass_AA['G'] - 57.02146) < 0.0001
    
    def test_mz_max(self):
        """Test MZ_MAX value."""
        assert config.MZ_MAX == 8000.0

