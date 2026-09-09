# Query-Patch Model History and Sweep Index

Last reviewed: 2026-09-09

This document records how the query-patch model evolved and maps each experimental
stage to the `max_epochs` value used to identify its artifacts under
[`jobs/training_pipeline/results`](../jobs/training_pipeline/results). It covers
architecture, training objectives, evaluation changes, and the sweeps used to test
them.

## Selecting the Final Model Variants

The [thesis](https://sharelatex.tum.de/read/xywkwfhfxsfc) variants are selectable
through the training CLI. All four support either `--int-mod-type normal` or
`--int-mod-type max` for interaction prediction.

| Thesis variant | CLI settings |
|---|---|
| D-SCRIPT | `--model dscript` |
| QP-Additive | `--model query_patch --query-evidence-mode additive --compatibility-head bilinear` |
| QP-Bilin-IC | `--model query_patch --query-evidence-mode gated --compatibility-head bilinear` |
| QP-Conv-IC | `--model query_patch --query-evidence-mode gated --compatibility-head dscript` |

For example, run QP-Additive with balanced contact BCE and max interaction pooling:

```bash
python -m scripts.training_pipeline.execute_pipeline \
    --model query_patch \
    --query-evidence-mode additive \
    --compatibility-head bilinear \
    --compatibility-scale 2 \
    --loss-type balanced_bce \
    --int-mod-type max \
    --max-epochs 100
```

`--query-evidence-mode` defaults to `gated`, preserving the current calculation.
`additive` restores the fixed, normalized formulation from stage 3 below:
`contact_logits = query_evidence + scale * compatibility_evidence + contact_bias`.
It uses the main interaction/contact objectives only: both auxiliary interface
and standalone compatibility losses are disabled, and gate-bias/interface-loss
settings have no effect. The bilinear scale stays fixed; scale `0` gives the
query-only ablation. Both query and compatibility branches still receive gradients
through the final contact logits.

The mode can also be set in a W&B sweep using `query-evidence-mode` as the
parameter name. Checkpoints save the mode, and automatically generated result
directories include `qev=additive` or `qev=gated`. Names exceeding the filesystem's
255-byte limit are shortened with a hash of the full name so long sweep values
still identify distinct runs. Existing gated checkpoints without this setting
still load as gated. When explicitly loading an older
additive Lightning checkpoint, pass `query_evidence_mode="additive"` to
`QueryPatchLightningModule.load_from_checkpoint(...)`.

Additive combination also supports `--compatibility-head dscript` as an additional
ablation, but that combination is not the thesis's historical QP-Additive model.
The convolutional head does not use `--compatibility-rank` or
`--compatibility-scale`.

## Current Model at a Glance

The query-patch model reuses the D-SCRIPT residue projection and interaction-pooling
heads, but replaces the D-SCRIPT contact-map generator with learned query patches.
For one protein pair, the default gated design is:

1. Project both proteins' residue embeddings with the shared D-SCRIPT projection.
2. Let learned queries attend jointly to the two projected proteins.
3. Decode one soft residue mask on each protein and a positive strength per query.
4. Combine matched masks into residue-pair query evidence.
5. Convert that evidence into a soft gate over residue pairs.
6. Predict residue-pair compatibility with a signed low-rank bilinear head.
7. Compute the final contact probability as
   `sigmoid(query_gate_logits) * sigmoid(compatibility_contact_logits)`.
8. Pool the resulting contact map with the inherited D-SCRIPT interaction head.

The product in step 7 makes query evidence a prerequisite: the compatibility branch
cannot predict a high contact probability where the query gate is closed. The main
implementation is in [`scripts/models/query_patch.py`](../scripts/models/query_patch.py),
and its training wrapper is in
[`scripts/training_pipeline/execute_pipeline.py`](../scripts/training_pipeline/execute_pipeline.py).

## Chronological Change History

### 1. Initial learned-query contact generator

**Commit:** `4fbe2c8` — 2026-08-04

The initial model introduced a Transformer decoder with `K` learned queries. Each
query produced:

- a soft residue mask for protein 1;
- a soft residue mask for protein 2; and
- a non-negative patch strength.

The masks were combined efficiently as a sum of rank-one patches:
`A^T diag(strength) B / sqrt(K) + contact_bias`. This generated a dense contact map
without materializing a `[batch, queries, length_1, length_2]` tensor. The existing
D-SCRIPT projection and interaction heads remained in use.

This stage was tested by the **max-97 query-patch sweep**.

### 2. Additive signed compatibility branch

**Commit:** `2c615aa` — 2026-08-05

Query patches could identify plausible interface regions, but could not explicitly
decide whether residues in two such regions were compatible. A signed low-rank
bilinear branch was therefore added:

- residues were projected into a configurable compatibility rank;
- an alternating positive/negative trainable metric allowed matching and opposing
  residue features; and
- a trainable scale controlled the compatibility contribution.

The contact logits became
`query_evidence + scale * compatibility_evidence + contact_bias`.

This stage was tested by the **max-96 query-patch sweep**.

### 3. Fixed and normalized compatibility

**Commit:** `8ad3bcb` — 2026-08-11

The compatibility branch was stabilized by:

- L2-normalizing the projected compatibility features;
- bounding the learned signed metric with `tanh`; and
- replacing the trainable compatibility scale with a fixed configuration value.

Fixing the scale prevented training from silently disabling the compatibility
branch. The sweep still included scale `0` as an explicit query-only ablation,
alongside scales `1`, `2`, `5`, `10`, and `20`.

This stage was tested by the broad **max-95 query-patch sweep**.

### 4. Class-balanced contact objectives

**Commit:** `4318c37` — 2026-08-13

Because contacts are extremely sparse relative to known non-contacts, two new
contact-map objectives were introduced:

- `balanced_bce`, which averages positive and negative BCE independently within
  each supervised map; and
- `balanced_bce_dice`, which adds a soft Dice term.

The positive-class weight is controlled by `omega`. These losses were tested on
both query patch and D-SCRIPT so that loss effects could be separated from
architecture effects.

They were tested first by the **max-94 fixed grid**, then by the **max-93 Bayesian
sweep over `omega`**.

### 5. Query evidence changed from an additive term to a gate

**Commit:** `f49e619` — 2026-08-19

The additive formulation allowed the compatibility branch to dominate while the
query decoder was ignored. Query evidence was therefore converted into a learned
gate with a positive scale and trainable bias. The final probability became:

```text
p(contact) = sigmoid(query_gate_logits)
           * sigmoid(compatibility_contact_logits)
```

The first version of this change directly supervised pairwise query-gate logits
against known contact-map cells with an auxiliary balanced BCE.

This direct-gate objective was tested by the **max-92 auxiliary-loss sweep**.

### 6. Auxiliary supervision changed to interfaces and compatibility

**Commit:** `2a11075` — 2026-08-19

Directly comparing the broad query gate with sparse residue-pair contacts did not
match the intended role of the queries. The auxiliary objective was redesigned:

- a residue is an interface target if any known partner cell is a contact;
- max-pooling query-mask logits over queries produces one interface logit per
  residue and chain;
- chain-wise interface masks receive their own balanced BCE; and
- the compatibility branch receives a separate contact-map BCE.

The final gated contact logits continue to receive the main contact loss. The
legacy direct-gate metrics and their current replacements are documented in
[`docs/training_metrics.md`](training_metrics.md#legacy-query-gate-metrics-in-older-runs).

This redesigned objective was tested by the **max-91 auxiliary-loss sweep**.

### 7. Contact and interaction AUPRC metrics

**Commit:** `b26a7e9` — 2026-08-20

Interaction AUPRC plus macro and micro contact AUPRC were added. These metrics are
more informative than thresholded accuracy for sparse contact maps:

- macro contact AUPRC gives each supervised complex equal influence; and
- micro contact AUPRC pools all known residue-pair cells.

### 8. Separate interface and contact class balance

**Commit:** `55791a2` — 2026-08-31

Interface residues are much denser than residue-pair contacts, so the single
`omega` value was split into:

- `interface_omega`, defaulting to `0.2`, for chain-wise interface BCE; and
- `omega` for final and compatibility contact-map BCE.

This commit also changed query-patch checkpointing and sweep selection to maximize
macro validation contact AUPRC, and changed the pipeline's default query count from
32 to 4. The max-91 grid fixes `interface_omega=0.2` while varying both interface
loss weight and contact `omega`.

## Sweep and Result-Directory Index

The configured epoch cap is intentionally used as an experiment-generation marker.
It is not the number of epochs that every run completed: early stopping usually
ended training earlier.

| `max_epochs` | Run date | W&B sweep | SLURM job | Method and change tested | Result-directory pattern | Artifacts observed |
|---:|---|---|---:|---|---|---:|
| 97 | 2026-08-03/04 | `3i1q6nkt` | `6432851` | Bayesian initial query-decoder search without compatibility. Varied interaction head, query count, heads, layers, dropout, and contact bias. | `querypatch-*_97` | 28 runs, 20 unique directories |
| 96 | 2026-08-04/05 | `9833yis3` | `6433993` | Bayesian additive-compatibility search. Added compatibility rank and trainable scale initialization. | `querypatch-*_96` | 60 directories |
| 95 | 2026-08-06 onward | `uhjcq9sp` | `6453061` | Broad Bayesian fixed-compatibility search over five loss types, query architecture, dropout, bias, rank, scale, and interaction head. | `95_model=query_patch_*` | 230 directories; 226 completed, 4 partial |
| 94 | 2026-08-12 | `v6bjofvo` | `6504538` | Fixed grid: query patch vs. D-SCRIPT, balanced BCE vs. balanced BCE + Dice, max vs. normal interaction head, `omega=0.5`. | `94_model=query_patch_*` | 4 query-patch / 8 total |
| 93 | 2026-08-13/14 | `k7vtr11y` | `6507382` | Bayesian continuation of max 94, now varying `omega` uniformly from `0.2` to `0.8`. | `93_model=query_patch_*` | 12 query-patch / 28 total |
| 92 | 2026-08-17 | `7pu2koy9` | `6529755` | Grid over direct pairwise gate-loss weight `[0, 0.05, 0.2, 0.5]` and gate bias `[-2, 0]`. | `92_model=query_patch_*qgatew=*` | 8 directories |
| 91 | 2026-08-20 | `11zro799` | `6548787` | `5 x 5` grid over interface-loss weight `[0, 0.05, 0.2, 0.5, 1]` and contact `omega` `[0.02, 0.05, 0.1, 0.2, 0.5]`, with `interface_omega=0.2`. | `91_model=query_patch_*qifw=*ifomega=0.2*` | 25 directories |

Representative logs containing the authoritative W&B configurations are:

- [max 97](../jobs/training_pipeline/logs/query-patch-sweep-6432851_0.err)
- [max 96](../jobs/training_pipeline/logs/query-patch-sweep-6433993_0.err)
- [max 95](../jobs/training_pipeline/logs/query-patch-sweep-6453061_0.err)
- [max 94](../jobs/training_pipeline/logs/balanced-bce-sweep-6504538_0.err)
- [max 93](../jobs/training_pipeline/logs/balanced-bce-sweep-6507382_0.err)
- [max 92](../jobs/training_pipeline/logs/aux-loss-sweep-6529755_0.err)
- [max 91](../jobs/training_pipeline/logs/aux-loss-sweep-6548787_0.err)

The currently tracked sweep definitions are:

- [`query_patch_sweep.yaml`](../jobs/training_pipeline/query_patch_sweep.yaml),
  currently retaining marker 95;
- [`balanced_bce_sweep.yaml`](../jobs/training_pipeline/balanced_bce_sweep.yaml),
  currently retaining marker 93; and
- [`aux_loss_sweep.yaml`](../jobs/training_pipeline/aux_loss_sweep.yaml),
  currently using marker 91.

## How to Locate a Generation Under `results`

The two oldest query-patch generations used the epoch marker as a suffix:

```bash
jobs/training_pipeline/results/querypatch-*_97
jobs/training_pipeline/results/querypatch-*_96
```

From max 95 onward, the output naming convention puts `max_epochs` first:

```bash
jobs/training_pipeline/results/95_model=query_patch_*
jobs/training_pipeline/results/94_model=query_patch_*
jobs/training_pipeline/results/93_model=query_patch_*
jobs/training_pipeline/results/92_model=query_patch_*
jobs/training_pipeline/results/91_model=query_patch_*
```

The naming implementation explicitly keeps `max_epochs` first for grouping; see
[`execute_pipeline.py`](../scripts/training_pipeline/execute_pipeline.py#L2286).

Lightning's progress display uses zero-based epoch indices. For example,
`max_epochs=95` appears as `Epoch .../94`. Use the W&B configuration or result
directory marker, not the progress-bar denominator, when identifying a sweep.

## Interpretation Caveats

- Several experiments were run from working-tree changes and committed only after
  the runs. Run logs are therefore the primary evidence for the exact historical
  configuration, while commit dates show when the code was recorded in Git.
- The exact max-94 grid YAML was not committed. Its eight logs and result
  directories are the authoritative record.
- The old max-97 output name omitted continuous parameters such as dropout and
  contact bias. Consequently, 28 runs produced only 20 unique result directories.
- The max-95 SLURM script allowed 280 array entries, but only 230 run logs/result
  directories are present; four of those directories are partial.
- Historical max-97 through max-92 runs selected and early-stopped on validation
  loss. The max-91 runs selected on macro validation contact AUPRC. The current
  YAML files were later updated to the AUPRC objective and therefore do not
  reproduce every historical sweep exactly.

## Fixed Runs That Were Not Sweeps

- [`esm2_query_patch.sh`](../jobs/training_pipeline/esm2_query_patch.sh) uses
  `max_epochs=89`. Its initial runs used explicit result names such as
  `test_query_patch`, so 89 is not necessarily visible in their directory names.
- [`test_query_patch.sh`](../jobs/training_pipeline/test_query_patch.sh) currently
  uses `max_epochs=80`. Earlier smoke tests also used two epochs.

## Working-Tree Change Not Yet Covered by a Sweep

As observed on 2026-08-31, the working tree contains an uncommitted alternative
compatibility head:

- `bilinear` retains the signed low-rank compatibility head used by the historical
  sweeps; and
- `dscript` applies the local pairwise D-SCRIPT convolutional contact module before
  the same query gate.

The option is exposed as `--compatibility-head {bilinear,dscript}`. No tracked
sweep YAML currently varies it, and no result generation can yet be attributed to
this change. This section should be updated with a commit and a new sweep marker
once the alternative is tested.

## Updating This Overview

For each future query-patch generation, record:

1. the commit hash and architectural or objective change;
2. the W&B sweep ID and SLURM array job ID;
3. the configured `max_epochs` marker;
4. the exact varied and fixed parameters;
5. the result-directory pattern and number of complete/partial artifacts; and
6. the checkpoint/sweep-selection metric used at run time.
