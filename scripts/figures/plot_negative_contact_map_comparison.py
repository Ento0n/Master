#!/usr/bin/env python3
"""Render an exploratory, paired comparison of validation contact maps.

The inference output supplies one unpadded [chain_1_length, chain_2_length]
probability map per checkpoint and sample. Positive reference maps encode
unknown/noncontact/contact as -1/0/1. Unknown positive cells are masked in every
panel; synthetic negative pairs have no measured map and retain all predicted
cells. Selection uses the saved per-pair differences; split provenance is
displayed from the inference manifest so historical split changes stay visible.
"""

import argparse
import json
from pathlib import Path
import textwrap

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


# =============================================================================
# Deterministic example selection and input checks
# =============================================================================

def select_samples(candidates):
    """Choose distinct maximum, class-median-nearest, and random examples.

    The median is calculated on the entire supplied class cohort. Already
    selected examples are excluded from subsequent choices. Ties follow the
    validation index, and each class uses its own fixed seed (42).
    """
    selected = []
    for label in (1, 0):
        group = candidates.loc[candidates['label'].eq(label)].sort_values(
            'validation_index', kind='stable').copy()
        if group.empty:
            continue
        median = group['mean_absolute_difference'].median()
        maximum = group['mean_absolute_difference'].idxmax()
        row = group.loc[maximum].copy()
        row['selection_reason'] = 'Largest mean absolute difference in class cohort'
        selected.append(row)
        remaining = group.drop(index=maximum)
        if not remaining.empty:
            nearest = (remaining['mean_absolute_difference'] - median).abs().idxmin()
            row = remaining.loc[nearest].copy()
            row['selection_reason'] = 'Nearest class median difference among remaining pairs'
            selected.append(row)
            remaining = remaining.drop(index=nearest)
        if not remaining.empty:
            choice = np.random.default_rng(42).choice(remaining.index.to_numpy())
            row = remaining.loc[choice].copy()
            row['selection_reason'] = 'Random remaining pair (seed 42)'
            selected.append(row)
    return pd.DataFrame(selected).reset_index(drop=True)


def read_maps(input_dir, row):
    """Validate shapes, probability ranges, and the -1/0/1 reference convention."""
    path = input_dir / str(row['map_file'])
    with np.load(path, allow_pickle=False) as archive:
        disabled, enabled, reference = (
            archive[key] for key in ('disabled', 'enabled', 'reference'))
    expected = (int(row['length_1']), int(row['length_2']))
    if any(array.shape != expected for array in (disabled, enabled, reference)):
        raise ValueError(f'{path}: maps must all have unpadded shape {expected}')
    for name, array in (('disabled', disabled), ('enabled', enabled)):
        if not np.isfinite(array).all() or np.any((array < 0) | (array > 1)):
            raise ValueError(f'{path}: {name} probabilities must be finite and in [0, 1]')
    if not np.isin(reference, (-1, 0, 1)).all():
        raise ValueError(f'{path}: reference must use -1/0/1 labels')
    if int(row['label']) == 0 and np.any(reference != -1):
        raise ValueError(f'{path}: negative pairs must have an unknown reference')
    return disabled, enabled, reference


# =============================================================================
# Contact-map panels: absolute probability and signed changes use fixed scales
# =============================================================================

def identity_text(row):
    """Keep source identifiers visible even when synthetic pairs have no PDB."""
    pdb_id = str(row.get('pdb_id', '')).strip()
    assembly = str(row.get('assembly_number', '')).strip()
    if assembly:
        assembly = f'{float(assembly):g}'
    source = f'PDB {pdb_id}' if pdb_id else 'Synthetic pair'
    if pdb_id and assembly:
        source += f', assembly {assembly}'
    return (f"Validation {int(row['validation_index'])} | {source} | "
            f"{row.get('cluster_1', '')} / {row.get('cluster_2', '')} | "
            f"{int(row['length_1'])} × {int(row['length_2'])} residues")


def comparison_figure(input_dir, rows, title, manifest):
    """Create one or more paired rows without normalizing individual maps.

    Axes refer to 1-based positions in each model input sequence, not PDB residue
    numbering. Pixel interpolation is disabled. The arrays remain at their
    original lengths; unknown positive cells appear gray in all four panels.
    """
    count = len(rows)
    subtitle_lines = [
        'Exploratory validation examples; selection shown per row. '
        'These examples are not population estimates.']
    for field in ('comparison_description', 'cohort_description'):
        description = str(manifest.get(field, '')).strip()
        if description:
            subtitle_lines.extend(textwrap.wrap(description, width=155))
    # Reserve physical space for provenance even on a single-example figure;
    # long split-history notes must not collide with the first sample heading.
    extra_lines = max(0, len(subtitle_lines) - 2)
    height = 3.9 * count + 2.0 + 0.18 * extra_lines
    figure = plt.figure(figsize=(15, height), facecolor='white')
    grid = figure.add_gridspec(
        2 * count + 1, 4, height_ratios=[0.22, 1.0] * count + [0.07],
        left=0.075, right=0.97, bottom=max(0.10, 1.15 / height),
        top=1 - (1.05 + 0.15 * extra_lines) / height,
        wspace=0.26, hspace=0.68)
    figure.suptitle(title, fontsize=17, fontweight='bold', y=1 - 0.10 / height)
    figure.text(0.5, 1 - 0.47 / height, '\n'.join(subtitle_lines),
                ha='center', va='top', fontsize=9)

    probability_cmap = plt.get_cmap('viridis').copy()
    difference_cmap = plt.get_cmap('RdBu_r').copy()
    reference_cmap = ListedColormap(['white', 'black'])
    for cmap in (probability_cmap, difference_cmap, reference_cmap):
        cmap.set_bad('#bcbcbc')

    for position, (_, row) in enumerate(rows.iterrows()):
        disabled, enabled, reference = read_maps(input_dir, row)
        positive = int(row['label']) == 1
        # Structural coverage defines comparable cells for positives. Negative
        # pairs have no structural target, so all unpadded predictions are shown.
        unknown = reference < 0 if positive else np.zeros(reference.shape, dtype=bool)
        valid = ~unknown
        maps = [np.ma.masked_where(reference < 0, reference),
                np.ma.masked_where(unknown, disabled),
                np.ma.masked_where(unknown, enabled),
                np.ma.masked_where(unknown, enabled - disabled)]
        heading = figure.add_subplot(grid[2 * position, :])
        heading.axis('off')
        heading.text(0, 0.78, identity_text(row), va='center', fontsize=10,
                     fontweight='bold', transform=heading.transAxes)
        heading.text(0, -0.05, row['selection_reason'], va='center', fontsize=9,
                     color='#4c4c4c', transform=heading.transAxes)
        titles = ['Measured reference', 'Negative-map supervision OFF',
                  'Negative-map supervision ON', 'ON minus OFF']
        for column, array in enumerate(maps):
            axis = figure.add_subplot(grid[2 * position + 1, column])
            cmap = (reference_cmap if column == 0 else
                    difference_cmap if column == 3 else probability_cmap)
            norm = Normalize(-1, 1) if column == 3 else Normalize(0, 1)
            # Fill the same panel area even for very unequal chain lengths.
            # Explicit residue axes retain the true dimensions; rectangular
            # display pixels do not resample or normalize the underlying map.
            axis.imshow(array, cmap=cmap, norm=norm, origin='upper',
                        interpolation='none', aspect='auto',
                        extent=(0.5, array.shape[1] + 0.5, array.shape[0] + 0.5, 0.5))
            axis.set_xlabel('Chain 2 residue index', fontsize=8)
            axis.set_ylabel('Chain 1 residue index' if column == 0 else '', fontsize=8)
            axis.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
            axis.tick_params(labelsize=7, length=3)
            if column in (1, 2) and valid.any():
                values = disabled[valid] if column == 1 else enabled[valid]
                stats = (f'Mean p = {values.mean():.4f}; '
                         f'p ≥ 0.5: {(values >= 0.5).mean():.1%}')
            elif column == 3 and valid.any():
                difference = (enabled - disabled)[valid]
                stats = (f'Mean Δp = {difference.mean():+.4f}; '
                         f'mean |Δp| = {np.abs(difference).mean():.4f}')
            elif column == 0 and positive:
                stats = f'0 noncontact; 1 contact; {valid.mean():.1%} measured'
            else:
                stats = ''
            axis.set_title(titles[column] + '\n' + stats, fontsize=9, pad=7)
            if column == 0 and not positive:
                axis.text(0.5, 0.5,
                          'No measured contact map\nNon-interacting pair (label 0)',
                          ha='center', va='center', transform=axis.transAxes,
                          fontsize=9, linespacing=1.6,
                          bbox=dict(facecolor='white', edgecolor='none', alpha=0.9))
            for spine in axis.spines.values():
                spine.set_color('#777777')
                spine.set_linewidth(0.5)

    # One shared color scale makes absolute changes comparable across samples
    # and checkpoints. The difference limits are fixed at the full [-1, 1].
    for columns, cmap, norm, label, ticks in (
            (slice(1, 3), probability_cmap, Normalize(0, 1),
             'Predicted contact probability', [0, 0.25, 0.5, 0.75, 1]),
            (3, difference_cmap, Normalize(-1, 1),
             'Contact probability change (ON − OFF)', [-1, -0.5, 0, 0.5, 1])):
        color_axis = figure.add_subplot(grid[-1, columns])
        colorbar = figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                                  cax=color_axis, orientation='horizontal', ticks=ticks)
        colorbar.set_label(label, fontsize=9)
        colorbar.ax.tick_params(labelsize=8)
    figure.text(0.075, 0.025,
                'Gray = unknown structural coverage. Positive-map statistics use measured cells only; '
                'negative-map statistics use every unpadded residue pair.\n'
                'Residue indices follow the model input sequences. The 0.5 cutoff is descriptive; '
                'it is not a calibrated interaction decision threshold.',
                ha='left', va='bottom', fontsize=8, color='#444444')
    return figure


def cohort_figure(candidates, manifest):
    """Show each candidate pair's OFF/ON mean on the same [0, 1] axes.

    Points are paired checkpoint evaluations, not independent observations from
    separate model cohorts. The diagonal marks equal mean contact probability.
    This exploratory cohort does not establish validation-wide performance.
    """
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 5.5))
    for axis, label, color in zip(axes, (1, 0), ('#277da1', '#c57035')):
        group = candidates.loc[candidates['label'].eq(label)]
        axis.plot([0, 1], [0, 1], color='#858585', linestyle='--', lw=1, zorder=1)
        axis.scatter(group['mean_disabled'], group['mean_enabled'],
                     s=28, c=color, alpha=0.75, linewidths=0.3,
                     edgecolors='white', zorder=2)
        axis.set(xlim=(0, 1), ylim=(0, 1), aspect='equal',
                 xlabel='Mean contact probability — OFF',
                 ylabel='Mean contact probability — ON',
                 title=f"{'Interacting' if label else 'Non-interacting'} pairs (label {label}, n = {len(group)})")
        axis.grid(color='#dddddd', linewidth=0.5, alpha=0.7)
        axis.set_axisbelow(True)
        if not group.empty:
            delta = group['mean_enabled'] - group['mean_disabled']
            axis.text(0.03, 0.97, f'Median paired change: {delta.median():+.4f}',
                      transform=axis.transAxes, fontsize=9, va='top')
    figure.suptitle('Exploratory validation cohort: paired contact-map means',
                     fontsize=14, fontweight='bold')
    description = str(manifest.get('cohort_description', '')).strip()
    footer = ('Each point is one pair evaluated with both checkpoints; dashed line: ON = OFF.\n'
              'Positive means use measured cells; negative means use all unpadded cells. '
              'This cohort is not a population estimate.')
    if description:
        footer += '\n' + '\n'.join(textwrap.wrap(description, width=130))
    figure.text(0.07, 0.025, footer, ha='left', va='bottom', fontsize=8)
    figure.subplots_adjust(left=0.075, right=0.985, bottom=0.24, top=0.85, wspace=0.28)
    return figure


# =============================================================================
# Standalone rendering entry point
# =============================================================================

def render(input_dir):
    """Write selection records, class/individual figures, and a multipage PDF."""
    input_dir = Path(input_dir).resolve()
    manifest_path = input_dir / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    candidates = pd.read_csv(input_dir / 'candidates.tsv', sep='\t', keep_default_na=False,
                             dtype={'cluster_1': str, 'cluster_2': str,
                                    'pdb_id': str, 'assembly_number': str})
    required = {'validation_index', 'label', 'length_1', 'length_2', 'mean_disabled',
                'mean_enabled', 'mean_absolute_difference', 'map_file'}
    if not required.issubset(candidates.columns):
        raise ValueError(f'Missing candidate columns: {sorted(required - set(candidates.columns))}')
    if candidates.empty or candidates['validation_index'].duplicated().any():
        raise ValueError('Candidates must contain at least one pair and unique validation indices')
    if not candidates['label'].isin((0, 1)).all():
        raise ValueError('Candidate labels must be 0 or 1')
    for name in ('mean_disabled', 'mean_enabled', 'mean_absolute_difference'):
        if not np.isfinite(candidates[name]).all() or not candidates[name].between(0, 1).all():
            raise ValueError(f'Candidate {name} values must be finite and in [0, 1]')

    figure_dir = input_dir / 'figures'
    figure_dir.mkdir(parents=True, exist_ok=True)
    selected = select_samples(candidates)
    selected_path = input_dir / 'selected_samples.tsv'
    selected.to_csv(selected_path, sep='\t', index=False)
    outputs = [selected_path]
    combined_path = input_dir / 'contact_map_comparison.pdf'
    # Rendering is local and uses only saved arrays; no checkpoint or training
    # state is loaded. Text remains selectable in the manuscript PDF/SVG files.
    with plt.rc_context({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'svg.fonttype': 'none'}):
        with PdfPages(combined_path) as combined:
            for label, name, title in (
                    (1, 'positives_comparison', 'Interacting validation pairs: negative-map supervision OFF versus ON'),
                    (0, 'negatives_comparison', 'Non-interacting validation pairs: negative-map supervision OFF versus ON')):
                group = selected.loc[selected['label'].eq(label)]
                if group.empty:
                    continue
                figure = comparison_figure(input_dir, group, title, manifest)
                for suffix in ('png', 'pdf', 'svg'):
                    path = input_dir / f'{name}.{suffix}'
                    figure.savefig(path, dpi=180, facecolor='white')
                    outputs.append(path)
                combined.savefig(figure)
                plt.close(figure)
                for index, row in group.iterrows():
                    figure = comparison_figure(input_dir, selected.loc[[index]], title, manifest)
                    for suffix in ('png', 'pdf', 'svg'):
                        path = figure_dir / f"val_{int(row['validation_index']):05d}.{suffix}"
                        figure.savefig(path, dpi=180, facecolor='white')
                        outputs.append(path)
                    plt.close(figure)
            figure = cohort_figure(candidates, manifest)
            for suffix in ('png', 'pdf'):
                path = input_dir / f'cohort_summary.{suffix}'
                figure.savefig(path, dpi=180, facecolor='white')
                outputs.append(path)
            combined.savefig(figure)
            plt.close(figure)
    outputs.append(combined_path)
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', required=True, type=Path,
                        help='Inference output containing candidates.tsv, manifest.json and maps/')
    args = parser.parse_args()
    for path in render(args.input_dir):
        print(path)


if __name__ == '__main__':
    main()
