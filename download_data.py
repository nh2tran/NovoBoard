#!/usr/bin/env python
"""Download NovoBoard example data from Google Drive.

This script downloads the ABRF dataset used for NovoBoard evaluation.
Data source: https://drive.google.com/drive/folders/1_6azR4-YjTUfRYdsXbFhZL9lFvjdrIDh

Usage:
    python download_data.py [--output-dir data]
    
Or via CLI:
    novoboard --download --data-dir data
"""

import argparse
import os
import sys

try:
    import gdown
except ImportError:
    print("Error: gdown is required. Install it with: uv add gdown")
    sys.exit(1)


# Google Drive folder containing the ABRF dataset
GDRIVE_FOLDER_URL = "https://drive.google.com/drive/folders/1_6azR4-YjTUfRYdsXbFhZL9lFvjdrIDh"


def download_data(output_dir: str = "data") -> None:
    """Download NovoBoard example data from Google Drive.
    
    Args:
        output_dir: Directory to save downloaded data (default: 'data')
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    print(f"Downloading NovoBoard data to '{output_dir}'...")
    print(f"Source: {GDRIVE_FOLDER_URL}")
    print()
    
    gdown.download_folder(GDRIVE_FOLDER_URL, output=output_dir, quiet=False)
    
    print()
    print("Download complete!")
    print(f"Data saved to: {output_dir}")


def main():
    """Main entry point for the download script."""
    parser = argparse.ArgumentParser(
        description="Download NovoBoard example data from Google Drive",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python download_data.py --output-dir data

This downloads the ABRF dataset used in the NovoBoard paper for
evaluating de novo peptide sequencing methods.
        """
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="data",
        help="Directory to save downloaded data (default: data)"
    )
    
    args = parser.parse_args()
    download_data(args.output_dir)


if __name__ == "__main__":
    main()

