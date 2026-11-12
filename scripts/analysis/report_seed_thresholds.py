#!/usr/bin/env python3
"""Verify embedding comparisons and export per-seed and aggregate test CSVs.

Only stdlib CSV/JSON readers are used: no model is loaded and no cutoff is fitted
on test labels. Accuracy, confusion counts and tie-aware average precision are
recomputed from saved predictions before any report is written. Aggregates use
three training seeds evaluated on the same held-out pairs. The per-seed CSV
combines embeddings and contains run identifiers, each fitted cutoff, its
validation accuracy, test AUPRC and test accuracies at that cutoff. A separate
aggregate CSV contains one row per embedding/model with test means and sample
standard deviations. All metrics are fractions on a 0-1 scale. Other validation
metrics and provenance remain part of verification, without report columns.
"""

import argparse
import csv
import hashlib
from itertools import groupby
import json
import math
from pathlib import Path
import statistics


MODELS = {'dscript': 'dscript', 'qp_conv_ic': 'query_patch'}
SEEDS = (42, 123, 456)
IDENTITY = ('split', 'new_cluster_pair', 'cluster_1', 'cluster_2',
            'pdb_id', 'assembly_number', 'dimer_type', 'label')


# =============================================================================
# Independent checks against saved prediction tables and provenance
# =============================================================================

def read_tsv(path):
    """Keep identifiers as strings; empty synthetic-negative IDs stay empty."""
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def digest(path):
    """Hash in bounded chunks so the same helper handles tables and source."""
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def metrics(rows, cutoff):
    """Return binary metrics with positive := probability >= cutoff, label 0/1."""
    counts = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}
    for row in rows:
        positive = float(row['predicted_probability']) >= cutoff
        actual = int(row['label']) == 1
        counts[('t' if positive == actual else 'f') + ('p' if positive else 'n')] += 1
    tp, tn, fp, fn = (counts[key] for key in ('tp', 'tn', 'fp', 'fn'))
    return dict(n=len(rows), accuracy=(tp + tn) / len(rows),
                precision=tp / (tp + fp) if tp + fp else 0.0,
                recall=tp / (tp + fn) if tp + fn else 0.0,
                f1=2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0,
                predicted_positive=tp + fp, **counts)


def average_precision(rows):
    """Integrate precision at each recall increment; equal scores enter together."""
    ordered = sorted(rows, key=lambda row: -float(row['predicted_probability']))
    positives = sum(int(row['label']) for row in rows)
    tp = seen = 0
    terms = []
    for _, tied in groupby(ordered, key=lambda row: float(row['predicted_probability'])):
        tied = list(tied)
        added = sum(int(row['label']) for row in tied)
        tp += added
        seen += len(tied)
        terms.append((added / positives) * (tp / seen))
    return math.fsum(terms)


def require_equal(actual, expected, context):
    """Allow roundoff only for derived floating metrics; counts remain exact."""
    same = (math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12)
            if isinstance(actual, float) else actual == expected)
    if not same:
        raise ValueError(f'{context}: actual {actual!r} != recorded {expected!r}')


def verify_run(directory, summary):
    """Check both splits, subgroups, original test scores and the 32-pair probe."""
    cutoff = summary['threshold']
    threshold = json.loads((directory / 'threshold.json').read_text())
    # The validation-only record freezes the rule and the evaluated checkpoint.
    # It must not contain test results or disagree with the later summary.
    if 'test' in threshold:
        raise ValueError(f'{directory.name}: frozen threshold record contains test results')
    for key in ('threshold', 'selection_split', 'selection_metric', 'tie_break',
                'model', 'seed', 'checkpoint', 'checkpoint_sha256', 'source_root',
                'source_sha256', 'input_sha256', 'submission_record', 'submission_record_sha256'):
        require_equal(threshold[key], summary[key], f'{directory.name}: frozen {key}')
    if summary['selection_split'] != 'val' or summary['selection_metric'] != 'accuracy':
        raise ValueError(f'{directory}: threshold must be fitted on validation accuracy')
    require_equal(summary['tie_break'], 'closest to 0.5, then smaller threshold',
                  f'{directory.name}: threshold tie break')
    require_equal(summary['model_hyperparameters']['interaction_module_type'], 'normal',
                  f'{directory.name}: interaction head')
    if not math.isfinite(cutoff):
        raise ValueError(f'{directory}: nonfinite threshold')
    require_equal(digest(Path(summary['checkpoint'])), summary['checkpoint_sha256'],
                  f'{directory.name}: checkpoint hash')
    tables = {}
    for split, expected_n in (('validation', 1480), ('test', 950)):
        path = directory / f'{split}_predictions.tsv'
        rows = read_tsv(path)
        require_equal(len(rows), expected_n, f'{directory.name}/{split} cohort size')
        require_equal(digest(path), summary['accuracy_reporting']['prediction_sha256'][path.name], str(path))
        for row in rows:
            label = int(row['label'])
            probability = float(row['predicted_probability'])
            if label not in (0, 1) or not (0 <= probability <= 1):
                raise ValueError(f'{path}: invalid binary label or probability')
            require_equal(row['split'], 'val' if split == 'validation' else 'test', str(path))
            require_equal(int(row['predicted_label']), int(probability >= 0.5), f'{path}: original label')
            require_equal(int(row['predicted_label_selected_threshold']), int(probability >= cutoff),
                          f'{path}: selected label')
        require_equal(set(row['dimer_type'] for row in rows), {'homo', 'hetero'}, f'{path}: groups')
        expected_sizes = {'overall': expected_n, 'homo': 924 if split == 'validation' else 606,
                          'hetero': 556 if split == 'validation' else 344}
        for group in ('overall', 'homo', 'hetero'):
            subset = rows if group == 'overall' else [row for row in rows if row['dimer_type'] == group]
            require_equal(len(subset), expected_sizes[group], f'{path}/{group}: group size')
            require_equal(sum(int(row['label']) for row in subset), len(subset) // 2,
                          f'{path}/{group}: balanced labels')
            recorded = summary[split] if group == 'overall' else summary[split]['by_dimer_type'][group]
            for decision, value in (('at_0.5', 0.5), ('at_selected_threshold', cutoff)):
                for key, actual in metrics(subset, value).items():
                    require_equal(actual, recorded[decision][key], f'{directory.name}/{split}/{group}/{decision}/{key}')
        require_equal(average_precision(rows), summary[split]['interaction_auprc'], f'{directory.name}/{split}/AP')
        tables[split] = rows
    for key, expected in threshold['validation'].items():
        require_equal(summary['validation'][key], expected, f'{directory.name}: frozen validation {key}')

    # Test output must preserve the original probabilities and identifiers. The
    # verification probe may differ within its recorded numeric tolerance only.
    original_path = Path(summary['test_predictions'])
    require_equal(digest(original_path), summary['test_predictions_sha256'], 'original test file hash')
    original = read_tsv(original_path)
    require_equal(len(original), len(tables['test']), 'original test row count')
    for original_row, row in zip(original, tables['test']):
        for key in IDENTITY + ('predicted_probability', 'predicted_label'):
            require_equal(row[key], original_row[key], f'{directory.name}: preserved {key}')
    probe = read_tsv(directory / 'test_score_verification.tsv')
    verification = summary['test_score_verification']
    require_equal(len(probe), 32, f'{directory.name}: probe length')
    require_equal(verification['n'], 32, f'{directory.name}: recorded probe length')
    errors = []
    disagreements = {'at_0.5': 0, 'at_selected_threshold': 0}
    for check, row in zip(probe, tables['test']):
        for key in IDENTITY:
            require_equal(check[key], row[key], f'{directory.name}: probe identity {key}')
        actual = float(check['predicted_probability'])
        saved = float(row['predicted_probability'])
        require_equal(float(check['original_predicted_probability']), saved, 'probe original probability')
        error = abs(actual - saved)
        if error > verification['absolute_tolerance'] + verification['relative_tolerance'] * abs(saved):
            raise ValueError(f'{directory.name}: failed reproduction tolerance')
        errors.append(error)
        for name, value in (('at_0.5', 0.5), ('at_selected_threshold', cutoff)):
            disagreements[name] += (actual >= value) != (saved >= value)
    require_equal(max(errors), verification['max_absolute_probability_error'], 'probe maximum error')
    require_equal(statistics.mean(errors), verification['mean_absolute_probability_error'], 'probe mean error')
    for name, count in disagreements.items():
        require_equal(count, verification[f'classification_disagreements_{name}'], f'probe {name}')
    return tables


# =============================================================================
# Comparison verification: each embedding has its own frozen inputs and source
# =============================================================================

def load_verified_comparison(analysis_root):
    """Verify all six runs for one embedding and return rows plus split identities."""
    if analysis_root.parent.name != 'threshold_selection':
        raise ValueError('--analysis-root must be <training>/threshold_selection/<analysis-id>')
    training_root = analysis_root.parents[1]
    embedding = training_root.parent.name
    if embedding not in ('esm2', 'prostt5'):
        raise ValueError('--analysis-root must belong to an esm2 or prostt5 comparison')

    directories = [analysis_root / f'{model}_normal_seed_{seed}' for model in MODELS for seed in SEEDS]
    missing = [str(path / 'summary.json') for path in directories
               if not (path / 'summary.json').exists() or not (path / 'accuracy_by_dimer_type.tsv').exists()]
    if missing:
        raise SystemExit('Awaiting all six completed analyses; no CSV written:\n' + '\n'.join(missing))

    records, summaries, cohorts = [], [], {}
    for directory in directories:
        summary = json.loads((directory / 'summary.json').read_text())
        tables = verify_run(directory, summary)
        for split, rows in tables.items():
            identities = [tuple(row[key] for key in IDENTITY) for row in rows]
            if split in cohorts:
                require_equal(identities, cohorts[split], f'{directory.name}: identical {split} cohort')
            cohorts[split] = identities
        model = 'dscript' if directory.name.startswith('dscript') else 'qp_conv_ic'
        seed = int(directory.name.rsplit('_', 1)[1])
        require_equal(summary['seed'], seed, f'{directory.name}: seed')
        require_equal(summary['model'], MODELS[model], f'{directory.name}: model')
        checkpoint_embeddings = [name for name in ('esm2', 'prostt5')
                                 if name in Path(summary['checkpoint']).parts]
        require_equal(checkpoint_embeddings, [embedding], f'{directory.name}: embedding')
        # The thesis uses the validation-fitted cutoff for all three test
        # groups. Keep its validation accuracy beside the cutoff to show the
        # selection score; AP uses the probability ranking and needs no cutoff.
        record = dict(embedding=embedding, model=MODELS[model], seed=seed,
                      threshold=summary['threshold'],
                      validation_accuracy=summary['validation']['at_selected_threshold']['accuracy'],
                      test_auprc=summary['test']['interaction_auprc'])
        for group in ('overall', 'homo', 'hetero'):
            block = summary['test'] if group == 'overall' else summary['test']['by_dimer_type'][group]
            record[f'test_{group}_accuracy'] = block['at_selected_threshold']['accuracy']
        records.append(record)
        summaries.append(summary)

    # All seeds must evaluate the same source and inputs. When the evaluation
    # received a submission record, also verify that record's own hash and its
    # source/input entries; older analyses can legitimately have no such record.
    source_hashes = summaries[0]['source_sha256']
    input_hashes = summaries[0]['input_sha256']
    for summary in summaries:
        require_equal(summary['source_sha256'], source_hashes, 'source hashes across runs')
        require_equal(summary['input_sha256'], input_hashes, 'input hashes across runs')
        require_equal(summary['source_root'], summaries[0]['source_root'], 'source root across runs')
        if summary['submission_record']:
            path = Path(summary['submission_record'])
            require_equal(digest(path), summary['submission_record_sha256'], f'submission record: {path}')
            submission = json.loads(path.read_text())
            for key, expected_hashes in (('source_sha256', source_hashes), ('input_sha256', input_hashes)):
                for filename, expected in expected_hashes.items():
                    require_equal(submission[key][filename], expected, f'submission {key}: {filename}')
        else:
            require_equal(summary['submission_record_sha256'], None, 'absent submission record hash')
    source_root = Path(summaries[0]['source_root'])
    for relative, expected in source_hashes.items():
        require_equal(digest(source_root / relative), expected, f'final source hash: {relative}')
    for filename, expected in input_hashes.items():
        require_equal(digest(Path(filename)), expected, f'final input hash: {filename}')
    return records, cohorts


# =============================================================================
# Reporting: one row per seed, with test means and sample SD in a separate CSV
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root', type=Path, action='append', required=True,
                        help='Completed <training>/threshold_selection/<analysis-id> directory; repeat to combine embeddings')
    parser.add_argument('--output-prefix', '--output-stem', dest='output_prefix', type=Path, required=True,
                        help='Output prefix for the per-seed .csv and _aggregate.csv files')
    args = parser.parse_args()
    records, reference_cohorts = [], None
    for root in args.analysis_root:
        rows, cohorts = load_verified_comparison(root.resolve(strict=True))
        # Embeddings change features, but all final runs must evaluate the
        # same ordered validation/test pairs before their results are combined.
        if reference_cohorts is not None:
            require_equal(cohorts, reference_cohorts, f'{root}: cohorts across embeddings')
        reference_cohorts = cohorts
        records.extend(rows)

    # Each input root has all six model/seed combinations. Reject duplicate
    # embeddings so repeated inputs cannot inflate the aggregate seed counts.
    combinations = [(row['embedding'], row['model'], row['seed']) for row in records]
    if len(combinations) != len(set(combinations)):
        raise ValueError('Duplicate embedding/model/seed combinations in the requested analyses')
    records.sort(key=lambda row: (row['embedding'], row['model'], row['seed']))

    # Match the contact report: one aggregate row per configuration, with
    # paired mean/sample-SD columns. Thresholds and selection scores stay in
    # the per-seed file because each trained model has its own fitted cutoff.
    aggregate_keys = [key for key in records[0] if key.startswith('test_')]
    aggregates = []
    for embedding in sorted({row['embedding'] for row in records}):
        for model in MODELS.values():
            group = [row for row in records if row['embedding'] == embedding and row['model'] == model]
            aggregate = dict(embedding=embedding, model=model, n_seeds=len(group))
            for key in aggregate_keys:
                values = [row[key] for row in group]
                aggregate[f'{key}_mean'] = statistics.fmean(values)
                aggregate[f'{key}_sample_sd'] = statistics.stdev(values)
            aggregates.append(aggregate)

    # All inputs and results are verified before either CSV report is written.
    prefix = args.output_prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = prefix.with_suffix('.csv')
    aggregate_path = prefix.with_name(prefix.name + '_aggregate.csv')
    for path, rows in ((csv_path, records), (aggregate_path, aggregates)):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(f'Verified {len(records)} analyses; wrote {csv_path} and {aggregate_path}')


if __name__ == '__main__':
    main()
