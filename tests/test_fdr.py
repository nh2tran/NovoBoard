"""Tests for FDR module."""

import os
import tempfile
import pytest
import pandas as pd
import numpy as np
from novoboard.fdr import read_denovo, calculate_FDR


class TestReadDenovo:
    """Test read_denovo function."""
    
    @pytest.fixture
    def sample_denovo_csv(self):
        """Create a minimal sample de novo CSV file for testing."""
        data = {
            'Source File': ['test.mgf', 'test.mgf', 'test.mgf'],
            'Scan': [1, 2, 3],
            'Peptide': ['PEPTIDE', 'PROTEIN', 'SAMPLE'],
            'ALC (%)': [95.0, 85.0, 75.0],
            'local confidence (%)': ['95,94,93,92,91,90,89', '85,84,83,82,81,80,79', '75,74,73,72,71,70']
        }
        df = pd.DataFrame(data)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df.to_csv(f.name, index=False)
            temp_path = f.name
        yield temp_path
        os.unlink(temp_path)
    
    def test_read_denovo_creates_feature_id(self, sample_denovo_csv):
        """Test that read_denovo creates feature_id column."""
        result = read_denovo(sample_denovo_csv)
        
        assert 'feature_id' in result.columns
        assert len(result) == 3
    
    def test_feature_id_format(self, sample_denovo_csv):
        """Test that feature_id has correct format."""
        result = read_denovo(sample_denovo_csv)
        
        # Format should be: source_file||scan
        expected_ids = ['test.mgf||1', 'test.mgf||2', 'test.mgf||3']
        assert list(result['feature_id']) == expected_ids
    
    def test_selected_features_filter(self, sample_denovo_csv):
        """Test that selected_features parameter filters correctly."""
        selected = {'test.mgf||1', 'test.mgf||3'}
        result = read_denovo(sample_denovo_csv, selected_features=selected)
        
        assert len(result) == 2
        assert set(result['feature_id']) == selected


class TestCalculateFDR:
    """Test calculate_FDR function."""
    
    @pytest.fixture
    def target_and_decoy_csvs(self):
        """Create sample target and decoy CSV files for testing."""
        # Target PSMs with higher scores
        target_data = {
            'Source File': ['test.mgf'] * 5,
            'Scan': [1, 2, 3, 4, 5],
            'Peptide': ['PEPTIDE', 'PROTEIN', 'SAMPLE', 'TARGET', 'MATCH'],
            'ALC (%)': [95.0, 90.0, 85.0, 80.0, 75.0],
        }
        target_df = pd.DataFrame(target_data)
        
        # Decoy PSMs with lower scores
        decoy_data = {
            'Source File': ['decoy.mgf'] * 5,
            'Scan': [1, 2, 3, 4, 5],
            'Peptide': ['EDITPEP', 'NIETORP', 'ELPMAS', 'TEGRAT', 'HCTAM'],
            'ALC (%)': [70.0, 65.0, 60.0, 55.0, 50.0],
        }
        decoy_df = pd.DataFrame(decoy_data)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            target_df.to_csv(f.name, index=False)
            target_path = f.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            decoy_df.to_csv(f.name, index=False)
            decoy_path = f.name
        
        yield target_path, decoy_path
        os.unlink(target_path)
        os.unlink(decoy_path)
    
    def test_calculate_fdr_returns_dataframes(self, target_and_decoy_csvs):
        """Test that calculate_FDR returns expected outputs."""
        target_csv, decoy_csv = target_and_decoy_csvs
        
        dfs, dfs_fdr, score_list, count_list = calculate_FDR(
            target_csv, decoy_csv, 'ALC (%)', [0.01, 0.05, 0.10]
        )
        
        assert isinstance(dfs, pd.DataFrame)
        assert isinstance(dfs_fdr, pd.DataFrame)
        assert isinstance(score_list, list)
        assert isinstance(count_list, list)
    
    def test_calculate_fdr_has_is_target_column(self, target_and_decoy_csvs):
        """Test that output has is_target column."""
        target_csv, decoy_csv = target_and_decoy_csvs
        
        dfs, dfs_fdr, _, _ = calculate_FDR(
            target_csv, decoy_csv, 'ALC (%)', [0.01, 0.05]
        )
        
        assert 'is_target' in dfs.columns
        assert 'is_target' in dfs_fdr.columns
    
    def test_calculate_fdr_has_estimated_fdr(self, target_and_decoy_csvs):
        """Test that output has estimated_fdr column."""
        target_csv, decoy_csv = target_and_decoy_csvs
        
        _, dfs_fdr, _, _ = calculate_FDR(
            target_csv, decoy_csv, 'ALC (%)', [0.01, 0.05]
        )
        
        assert 'estimated_fdr' in dfs_fdr.columns
    
    def test_score_list_length_matches_fdr_list(self, target_and_decoy_csvs):
        """Test that score_list has same length as fdr_list."""
        target_csv, decoy_csv = target_and_decoy_csvs
        fdr_list = [0.01, 0.05, 0.10]
        
        _, _, score_list, count_list = calculate_FDR(
            target_csv, decoy_csv, 'ALC (%)', fdr_list
        )
        
        assert len(score_list) == len(fdr_list)
        assert len(count_list) == len(fdr_list)

