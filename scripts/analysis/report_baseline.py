#!/usr/bin/env python3
"""Write the ONE-HOT baseline's shared-test results as per-seed and aggregate CSVs.

Saved interaction probabilities are checked against their original exports, then
accuracy and AUPRC are recomputed on the same 950 pairs for both dataset versions.
Each seed keeps its validation-selected cutoff and the accuracy at that cutoff.
No model is loaded or threshold fitted. A separate aggregate CSV contains one
row per dataset, with test means and sample standard deviations across seeds.
All metric values remain fractions in [0, 1], rather than percentages.
"""

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

if __package__:
    from .report_seed_thresholds import average_precision, digest, metrics, read_tsv, require_equal
else:
    from report_seed_thresholds import average_precision, digest, metrics, read_tsv, require_equal


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ANALYSIS = ROOT / 'jobs/training_pipeline/results/fully_connected_thresholds/20260916'
DATASETS = {'strict': 'strict_cd_hit', 'keep_homomers': 'keep_homomers_cd_hit'}
SEEDS = (42, 123, 456)
FIELDS = ['dataset', 'seed', 'threshold', 'validation_accuracy', 'test_auprc',
          'test_overall_accuracy_at_0_5', 'test_overall_accuracy_selected',
          'test_homo_accuracy_at_0_5', 'test_homo_accuracy_selected',
          'test_hetero_accuracy_at_0_5', 'test_hetero_accuracy_selected']


# =============================================================================
# Verify saved decisions and prediction provenance before selecting shared pairs
# =============================================================================

def verify_run(directory, dataset, seed):
    """Return preserved test rows and the matching frozen validation summary."""
    summary = json.loads((directory / 'summary.json').read_text())
    frozen = json.loads((directory / 'threshold.json').read_text())
    if 'test' in frozen or 'cohorts' in frozen:
        raise ValueError(f'{directory}: frozen cutoff contains test results')
    for key, value in frozen.items():
        require_equal(summary[key], value, f'{directory.name}: frozen {key}')
    for key, value in (('dataset', dataset), ('seed', seed), ('selection_split', 'val'),
                       ('selection_metric', 'unweighted_accuracy'),
                       ('decision_rule', 'probability >= threshold')):
        require_equal(frozen[key], value, f'{directory.name}: {key}')
    if not math.isfinite(frozen['threshold']):
        raise ValueError(f'{directory}: nonfinite cutoff')

    # The original validation export and inference record anchor the cutoff to
    # its checkpoint. These checks do not rerun inference or optimize a cutoff.
    for name, key in (('validation_predictions.tsv', 'validation_sha256'),
                      ('inference_provenance.json', 'inference_provenance_sha256')):
        require_equal(digest(directory / name), frozen[key], f'{directory.name}: {name} hash')
    validation = read_tsv(directory / 'validation_predictions.tsv')
    require_equal(metrics(validation, frozen['threshold'])['accuracy'],
                  frozen['validation_at_selected']['accuracy'], f'{directory.name}: validation accuracy')
    require_equal(digest(directory.parent / 'threshold_selection_source.py'),
                  frozen['selection_implementation_sha256'], 'frozen selector hash')
    provenance = json.loads((directory / 'inference_provenance.json').read_text())
    require_equal(provenance['seed'], seed, f'{directory.name}: inference seed')
    require_equal(provenance['protocol'], dataset, f'{directory.name}: inference dataset')
    require_equal(provenance['metric_verification']['passed'], True, 'validation verification')
    probe = json.loads((directory / 'test_score_verification.json').read_text())
    require_equal(probe['passed'], True, f'{directory.name}: saved inference probe')
    require_equal(probe['checkpoint_sha256'], provenance['checkpoint_sha256'], 'probe checkpoint')
    require_equal(probe['frozen_threshold_sha256'], digest(directory / 'threshold.json'), 'probe cutoff')

    original_path = Path(summary['test_predictions'])
    require_equal(digest(original_path), summary['original_test_sha256'], 'original test hash')
    require_equal(probe['test_export_sha256'], summary['original_test_sha256'], 'probe test hash')
    original = read_tsv(original_path)
    rows = read_tsv(directory / 'test_predictions.tsv')
    require_equal(len(rows), len(original), f'{directory.name}: preserved row count')
    require_equal(len(rows), 950 if dataset == 'strict_cd_hit' else 4054, 'native test size')
    require_equal(len({row['entity_pair'] for row in rows}), len(rows), 'unique test pairs')
    for row, saved in zip(rows, original):
        # Preserve every original identifier, label, and probability exactly;
        # the analysis adds only the decision at the already frozen cutoff.
        for key, value in saved.items():
            require_equal(row[key], value, f'{directory.name}: preserved {key}')
        probability, label = float(row['predicted_probability']), int(row['label'])
        if label not in (0, 1) or not 0 <= probability <= 1:
            raise ValueError(f'{directory.name}: invalid binary label or probability')
        require_equal(row['split'], 'test', f'{directory.name}: prediction split')
        require_equal(int(row['predicted_label']), int(probability >= 0.5), 'original decision')
        require_equal(int(row['predicted_label_selected_threshold']),
                      int(probability >= frozen['threshold']), 'frozen decision')
    return rows, summary


# =============================================================================
# Shared-cohort metrics and seed aggregation; export one CSV for each row format
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root', type=Path, default=DEFAULT_ANALYSIS,
                        help='Completed six-run baseline threshold analysis directory')
    parser.add_argument('--output-stem', '--output-prefix', type=Path,
                        default=ROOT / 'results/baseline_thresholds',
                        help='Output prefix for the per-seed .csv and _aggregate.csv files')
    args = parser.parse_args()
    analysis_root = args.analysis_root.resolve()

    # Strict seed 42 defines the common cohort; keep-homomers has additional
    # native test pairs that must not enter this direct dataset comparison.
    reference_summary = json.loads((analysis_root / 'strict_cd_hit_seed_42/summary.json').read_text())
    reference = read_tsv(Path(reference_summary['test_predictions']))
    reference_by_pair = {row['entity_pair']: row for row in reference}
    require_equal(len(reference), 950, 'shared reference size')
    require_equal(len(reference_by_pair), 950, 'unique shared reference pairs')

    records = []
    for name, dataset in DATASETS.items():
        for seed in SEEDS:
            rows, summary = verify_run(analysis_root / f'{dataset}_seed_{seed}', dataset, seed)
            by_pair = {row['entity_pair']: row for row in rows}
            shared = [by_pair[pair] for pair in reference_by_pair]
            for row in shared:
                for key in ('label', 'dimer_type'):
                    require_equal(row[key], reference_by_pair[row['entity_pair']][key], f'shared {key}')
            cutoff = summary['threshold']
            record = dict(dataset=name, seed=seed, threshold=cutoff,
                          validation_accuracy=summary['validation_at_selected']['accuracy'])
            for group, size in (('all', 950), ('homo', 606), ('hetero', 344)):
                subset = shared if group == 'all' else [row for row in shared if row['dimer_type'] == group]
                require_equal(len(subset), size, f'{dataset}/{seed}/{group}: shared size')
                require_equal(sum(int(row['label']) for row in subset), size // 2, 'balanced shared labels')
                saved = summary['cohorts']['shared_test'][group]
                # AP ranks all pair probabilities and is independent of either
                # cutoff. Accuracy is unweighted, including for keep-homomers.
                ap = average_precision(subset)
                require_equal(ap, saved['interaction_ap'], f'{dataset}/{seed}/{group}: AP')
                if group == 'all':
                    record['test_auprc'] = ap
                for rule, suffix, threshold in (('fixed_0.5', 'at_0_5', 0.5),
                                                ('validation_selected', 'selected', cutoff)):
                    computed = metrics(subset, threshold)
                    for key, value in computed.items():
                        require_equal(value, saved[rule][key], f'{dataset}/{seed}/{group}/{rule}/{key}')
                    column_group = 'overall' if group == 'all' else group
                    record[f'test_{column_group}_accuracy_{suffix}'] = computed['accuracy']
            records.append(record)

    # Average test metrics computed separately for each trained seed. The
    # cutoff and its validation score remain per-seed selection information.
    aggregates = []
    for dataset in DATASETS:
        selected = [row for row in records if row['dataset'] == dataset]
        aggregate = {'dataset': dataset, 'n_seeds': len(selected)}
        for key in FIELDS:
            if key.startswith('test_'):
                values = [row[key] for row in selected]
                aggregate[f'{key}_mean'] = statistics.fmean(values)
                aggregate[f'{key}_sample_sd'] = statistics.stdev(values)
        aggregates.append(aggregate)

    prefix = args.output_stem.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    output = prefix.with_suffix('.csv')
    aggregate_output = prefix.with_name(prefix.name + '_aggregate.csv')
    for path, rows in ((output, records), (aggregate_output, aggregates)):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(f'Verified all six baseline analyses; wrote {output} and {aggregate_output}')


if __name__ == '__main__':
    main()
