"""Command-line interface for NovoBoard."""

import argparse
import os
import sys

from novoboard.accuracy import WorkerTest
from novoboard.decoy import generate_decoy_mgf
from novoboard.fdr import validate_FDR
from novoboard.plotting import plot_fdr_validation


def download_data(data_dir):
    """Download data from Google Drive using gdown."""
    import gdown
    
    folder_url = "https://drive.google.com/drive/folders/1_6azR4-YjTUfRYdsXbFhZL9lFvjdrIDh"
    
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    
    print(f"Downloading data to {data_dir}...")
    gdown.download_folder(folder_url, output=data_dir, quiet=False)
    print("Download complete.")


def run_accuracy(data_dir, col_score, col_aa_score):
    """Run accuracy calculation (Cell 2 equivalent)."""
    folder = data_dir + '/'
    target_file = folder + 'pd_merged.csv.db.psms.csv'
    spectrum_file = folder + '2017-12-4_ABRF_200_DDA1.mgf'
    
    for x in range(10, 10+1):
        predicted_file = folder + 'PEAKS/Sample {0:s}.denovo.csv'.format(str(x))
        worker_test = WorkerTest(target_file, predicted_file, spectrum_file, col_score, col_aa_score)
        worker_test.test_accuracy()


def run_decoy_generation(data_dir, peak_sampling='random', sampling_rate=0.5):
    """Run decoy MGF generation (Cell 4 equivalent)."""
    folder = data_dir + '/'
    input_mgf_list = [
        '2017-12-4_ABRF_200_DDA1.mgf',
    ]
    input_mgf_list = [folder + x for x in input_mgf_list]
    
    generate_decoy_mgf(input_mgf_list, peak_sampling, sampling_rate)


def run_fdr_validation(data_dir, output_dir, col_score, col_aa_score):
    """Run FDR validation (Cell 8 equivalent)."""
    folder = data_dir + '/'
    db_csv = folder + 'pd_merged.csv.db.psms.csv'
    spectrum_file = folder + '2017-12-4_ABRF_200_DDA1.mgf'

    p_decoy = [x/1000. for x in range(0, 50, 1)]
    T_pct = 0.90

    samples = range(3, 7+1)
    target_csv = folder + 'PEAKS/Sample 10.denovo.csv'
    decoy_csv_list = [folder + 'PEAKS/Sample {0:s}.denovo.csv'.format(str(x)) for x in samples]
    engine_score = col_score
    
    results_list = [validate_FDR(target_csv, decoy_csv, engine_score, db_csv, spectrum_file, p_decoy, T_pct, col_score, col_aa_score) 
                    for decoy_csv in decoy_csv_list]

    output_path = os.path.join(output_dir, 'fig.decoy_fdr_valid_X_random_abrf_peaks.png')
    plot_fdr_validation(results_list, samples, output_path)


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description='NovoBoard - Framework for evaluating de novo peptide sequencing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  novoboard --data-dir data --output-dir fig
  novoboard --download --data-dir data --output-dir fig
        """
    )
    parser.add_argument(
        '--data-dir', 
        default='data',
        help='Path to data directory (default: data/)'
    )
    parser.add_argument(
        '--output-dir', 
        default='fig',
        help='Path for output files/plots (default: fig/)'
    )
    parser.add_argument(
        '--download',
        action='store_true',
        help='Download data using gdown before running'
    )
    parser.add_argument(
        '--skip-accuracy',
        action='store_true',
        help='Skip accuracy calculation step'
    )
    parser.add_argument(
        '--skip-decoy',
        action='store_true',
        help='Skip decoy MGF generation step'
    )
    parser.add_argument(
        '--skip-fdr',
        action='store_true',
        help='Skip FDR validation step'
    )
    
    args = parser.parse_args()
    
    # Column names for score values
    col_score = "ALC (%)"
    col_aa_score = "local confidence (%)"
    
    # Download data if requested
    if args.download:
        download_data(args.data_dir)
    
    # Check if data directory exists
    if not os.path.exists(args.data_dir):
        print(f"Error: Data directory '{args.data_dir}' does not exist.")
        print("Use --download to download the data first, or specify a valid --data-dir.")
        sys.exit(1)
    
    # Create output directory if needed
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Run all steps sequentially (like the notebook)
    if not args.skip_accuracy:
        print("\n" + "="*80)
        print("Step 1: Running accuracy calculation...")
        print("="*80)
        run_accuracy(args.data_dir, col_score, col_aa_score)
    
    if not args.skip_decoy:
        print("\n" + "="*80)
        print("Step 2: Generating decoy MGF...")
        print("="*80)
        run_decoy_generation(args.data_dir)
    
    if not args.skip_fdr:
        print("\n" + "="*80)
        print("Step 3: Running FDR validation...")
        print("="*80)
        run_fdr_validation(args.data_dir, args.output_dir, col_score, col_aa_score)
    
    print("\n" + "="*80)
    print("NovoBoard analysis complete!")
    print("="*80)


if __name__ == '__main__':
    main()

