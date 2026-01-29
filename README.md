# NovoBoard

A comprehensive framework for evaluating the false discovery rate and accuracy of de novo peptide sequencing.

## Features

- Calculate fragment ion, amino acid, and peptide accuracies
- Generate decoy spectra for FDR estimation
- Validate FDR estimation against known database matches
- Visualization of FDR validation results

## Installation

### Using uv (recommended)

```bash
# Clone the repository
git clone https://github.com/BioGeek/NovoBoard.git
cd NovoBoard

# Create virtual environment and install dependencies
uv sync

# Activate the virtual environment
source .venv/bin/activate
```

### Using pip

```bash
pip install -e .
```

## Data Download

The ABRF dataset can be downloaded from Google Drive:

```bash
# Using the provided script
python download_data.py --output-dir data

# Or via CLI
novoboard --download --data-dir data
```

Manual download link: https://drive.google.com/drive/folders/1_6azR4-YjTUfRYdsXbFhZL9lFvjdrIDh

## Usage

### Command Line Interface

Run all analysis steps sequentially (like the notebook):

```bash
novoboard --data-dir data --output-dir fig
```

Available options:

```
--data-dir      Path to data directory (default: data/)
--output-dir    Path for output files/plots (default: fig/)
--download      Download data using gdown before running
--skip-accuracy Skip accuracy calculation step
--skip-decoy    Skip decoy MGF generation step
--skip-fdr      Skip FDR validation step
```

### Jupyter Notebook

The original notebook `aa.fdr_github.ipynb` is still available for interactive analysis:

```bash
source .venv/bin/activate
jupyter notebook aa.fdr_github.ipynb
```

## Development

### Running Tests

```bash
source .venv/bin/activate
pytest tests/ -v
```

### Project Structure

```
novoboard/
├── src/novoboard/
│   ├── __init__.py
│   ├── cli.py           # Command-line interface
│   ├── config.py        # Vocabulary and mass definitions
│   ├── accuracy.py      # Accuracy calculation (WorkerTest)
│   ├── decoy.py         # Decoy MGF generation
│   ├── fdr.py           # FDR calculation and validation
│   └── plotting.py      # Visualization functions
├── tests/               # Unit tests
├── data/                # Data directory (not in git)
├── aa.fdr_github.ipynb  # Original Jupyter notebook
├── config.py            # Compatibility shim for notebook
├── download_data.py     # Data download script
└── pyproject.toml       # Project configuration
```

## Citation

If you use NovoBoard in your research, please cite:

> "NovoBoard: a comprehensive framework for evaluating the false discovery rate and accuracy of de novo peptide sequencing"
> https://doi.org/10.1101/2024.04.16.589668
