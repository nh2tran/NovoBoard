"""Plotting functions for FDR validation results."""

import os
from matplotlib import pyplot


def plot_fdr_validation(results_list, samples, output_path):
    """Plot FDR validation results.
    
    Args:
        results_list: List of validation result dictionaries from validate_FDR
        samples: Range of sample numbers used
        output_path: Path to save the output figure
    """
    # Create output directory if needed
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    fig, ax = pyplot.subplots(1, 2, figsize=(9,4))
    labels = ['decoy {0:s}0%'.format(str(10-x)) for x in samples]
    colors9 = ['b', 'c', 'g', 'k', 'm', 'r', 'y', 'orange', 'pink']
    colors = [colors9[x-1] for x in samples]
    ylim = 0
    for results, c, l in zip(results_list, colors, labels):
        ax[0].plot(results['estimated_fdr'], results['true_fdr_T'], color=c, label=l)
        ax[1].plot(results['cumsum'], results['estimated_fdr'], color=c, label=l)
    ax[0].plot([0, 0.05], [0, 0.05], color='black', linestyle='--', label='True FDR')
    ax[0].set_xlim(0, 0.05)
    # ax[0].set_ylim(0, 0.10)
    ax[0].set_xlabel('Estimated FDR'); ax[0].set_ylabel('True FDR')
    ax[0].legend()
    ax[1].plot(results_list[2]['cumsum'], results_list[2]['true_fdr_T'], color='black', label='True FDR', linestyle='--')
    ax[1].set_ylim(0, 0.05)
    ax[1].set_xlabel('Number of PSMs'); ax[1].set_ylabel('FDR')
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(output_path)
    print(f"Figure saved to {output_path}")
    pyplot.close(fig)

