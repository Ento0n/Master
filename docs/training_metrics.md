# Training metrics reference

This document explains every project-defined metric emitted by
[`scripts/training_pipeline/execute_pipeline.py`](../scripts/training_pipeline/execute_pipeline.py).
It covers the current fully
connected, D-SCRIPT, and query-patch training paths, followed by a short appendix
for metric names that still occur in older W&B runs.

The reference distinguishes four kinds of logged values:

- **quality metrics**, such as accuracy and AUPRC;
- **losses**, which are training objectives and are meaningful only together with
  their loss configuration;
- **model diagnostics**, such as query-gate scale and bias;
- **runtime metadata**, such as epoch, step, and contact-metric computation time.

## Quick conventions

- `{stage}` in a metric name expands to `train`, `val`, or `test`. For example,
  `{stage}_accuracy` represents `train_accuracy`, `val_accuracy`, and
  `test_accuracy`.
- All metrics produced through Lightning's `self.log` are epoch-level values:
  `on_step=False` and `on_epoch=True`. Training and validation values are emitted
  once per epoch; test values are emitted for the final test evaluation.
- Pair labels use `1 = interacting` and `0 = non-interacting`. Pair predictions
  always use a fixed interaction-probability threshold of `0.5`.
- Contact-map labels use `1 = contact`, `0 = known non-contact`, and
  `-1 = unknown or padding`. Unknown cells never enter contact losses or metrics.
- Thresholded contact metrics use `--contact-threshold`, whose default is `0.5`.
  Contact AUPRC does not use this threshold.
- **Average precision (AP)** is reported under the familiar `auprc` name. It is
  the step-wise precision-recall integral
  `AP = sum(delta_recall * precision)`, not a single-threshold measurement.
- Unless noted otherwise, losses are better when **lower**, while accuracy,
  precision, recall, and AUPRC are better when **higher**.

## Complete current metric inventory

Metric keys are conditional on the selected model and on the available data. A
missing validation/test split omits that stage. The homo/hetero metrics are
omitted if that subgroup is absent, and the positive/negative interface means are
omitted when their corresponding class is absent.

| Metric name or pattern | Fully connected | D-SCRIPT | Query patch | Stages |
|---|:---:|:---:|:---:|---|
| `{stage}_loss` | yes | yes | yes | train, val, test |
| `{stage}_accuracy` | yes | yes | yes | train, val, test |
| `{stage}_interaction_loss` | no | yes | yes | train, val, test |
| `val_interaction_auprc`, `test_interaction_auprc` | yes | yes | yes | val, test |
| `test_accuracy_homo`, `test_accuracy_hetero` | yes | yes | yes | test |
| `{stage}_contact_loss` | no | yes | yes | train, val, test |
| `{stage}_contact_accuracy` | no | yes | yes | train, val, test |
| `{stage}_contact_precision` | no | yes | yes | train, val, test |
| `{stage}_contact_recall` | no | yes | yes | train, val, test |
| `val_contact_auprc_macro`, `test_contact_auprc_macro` | no | yes | yes | val, test |
| `val_contact_auprc_micro`, `test_contact_auprc_micro` | no | yes | yes | val, test |
| `{stage}_contact_metrics_time_s` | no | yes | yes | train, val, test |
| `{stage}_optimization_loss` | no | no | yes | train, val, test |
| `{stage}_contact_objective` | no | no | yes | train, val, test |
| `{stage}_interface_mask_loss` | no | no | yes | train, val, test |
| `{stage}_interface_mask_weighted_loss` | no | no | yes | train, val, test |
| `{stage}_compatibility_contact_loss` | no | no | yes | train, val, test |
| `{stage}_interface_mask_mean` | no | no | yes | train, val, test |
| `{stage}_interface_mask_positive_mean` | no | no | yes | train, val, test |
| `{stage}_interface_mask_negative_mean` | no | no | yes | train, val, test |
| `{stage}_query_gate_scale` | no | no | yes | train, val, test |
| `{stage}_query_gate_bias` | no | no | yes | train, val, test |
| `best_val_loss` | yes | no | no | once after fitting; W&B only |
| `best_val_contact_auprc_macro` | no | yes | yes | once after fitting; W&B only |

When all conditional keys are present, the model families produce the following
number of project-defined metric keys, including the model family's `best_*`
selection metric:

| Model | Number of keys |
|---|---:|
| Fully connected | 11 |
| D-SCRIPT | 33 |
| Query patch | 63 |

There is deliberately no `train_interaction_auprc` or
`train_contact_auprc_*`. The AUPRC metrics are reserved for validation and test
evaluation.

## Interaction metrics and losses

Let `z_i` be a protein-pair interaction logit, `p_i = sigmoid(z_i)` its
probability, `y_i` its binary label, and `w_i` its optional loss weight.

### `{stage}_accuracy`

The fraction of protein pairs classified correctly:

```text
predicted_i = 1 if p_i >= 0.5 else 0
accuracy = number of correct predictions / number of pairs
```

Accuracy is unweighted even if interaction BCE uses dimer-type/label weights.
Lightning weights batch accuracies by batch size, so the epoch value is the exact
pair-level accuracy over the stage. Range: `[0, 1]`; higher is better.

Accuracy can conceal poor minority-class performance. Use
`val_interaction_auprc` or `test_interaction_auprc` when the positive/negative
classes are imbalanced.

### `test_accuracy_homo` and `test_accuracy_hetero`

The same `0.5`-threshold interaction accuracy, restricted respectively to test
homodimers and heterodimers. Rows whose dimer type is unknown are not included in
either subgroup. The epoch aggregation is weighted by the number of subgroup
members in each batch, so each result is the exact accuracy over that subgroup.

The metric is absent when its subgroup does not occur in the test split. Range:
`[0, 1]`; higher is better.

### `val_interaction_auprc` and `test_interaction_auprc`

Exact split-wide binary average precision over all protein-pair probabilities and
labels in the validation or test split. The complete split is ranked by predicted
probability, so the result is independent of the `0.5` classification threshold.
Loss weights are not applied: this metric describes the natural evaluation
population.

Range: `[0, 1]`; higher is better. A value should be interpreted relative to the
positive prevalence in the corresponding split.

### `{stage}_interaction_loss`

Logged by D-SCRIPT and query patch. It is the weighted mean pair-level binary
cross-entropy:

```text
interaction_loss = mean_i(w_i * BCEWithLogits(z_i, y_i))
```

In `basic` loss mode, every `w_i` is `1`. In `dimer_type` mode, the training-set
weight for a `(dimer_type, label)` cell is

```text
w_(d,y) = N_train / (C * n_(d,y))
```

where `C` is the number of observed training cells and `n_(d,y)` is that cell's
training count. The training-derived lookup is also used for validation and test
losses. Lower is better, but weighted losses should only be compared when the
weighting setup matches.

The fully connected model does not log a separate `interaction_loss`; its
`{stage}_loss` is this same interaction objective.

### `{stage}_loss`

This name provides a shared differentiable task objective.

For the fully connected model:

```text
stage_loss = interaction_loss
```

For D-SCRIPT and query patch, with
`lambda = --interaction-loss-lambda`:

```text
if the batch has known contact cells:
    stage_loss = lambda * interaction_loss
               + (1 - lambda) * contact_loss
else:
    stage_loss = interaction_loss
```

Query-patch auxiliary losses are deliberately excluded from `{stage}_loss`, even
though they participate in optimization. Ordinary batch values are aggregated
with the full protein-pair batch size. Lower is better.

For D-SCRIPT and query patch, `val_loss` drives the learning-rate scheduler.
Checkpointing and early stopping instead maximize
`val_contact_auprc_macro`. The fully connected model continues to minimize
`val_loss` for checkpointing and early stopping.

### `best_val_loss`

For the fully connected model, this is the lowest checkpoint-monitored
`val_loss` observed during fitting. It is logged once after training at the
final global step.

This is the selected checkpoint's validation score, not necessarily the last
epoch's `val_loss`. Lower is better. Unlike all `self.log` metrics,
`best_val_loss` is sent directly to W&B and is not written to the CSV metrics
file.

### `best_val_contact_auprc_macro`

For D-SCRIPT and query patch, this is the highest
`val_contact_auprc_macro` observed during fitting. It is used for checkpoint
selection, early stopping, and W&B contact-model sweeps.

This value is logged once after training at the final global step. Higher is
better. It is sent directly to W&B and is not written to the CSV metrics file.

## Contact-map population and masking

The following distinctions apply to every contact metric:

- A **known cell** has target `0` or `1`. A target of `-1` denotes an unresolved
  residue pair or batch padding and is always excluded.
- A **true structural map** was loaded from the contact-map data. Contact AUPRC
  evaluates only these maps.
- With `--negative-contact-maps`, an example without a structural map receives an
  optional synthetic all-zero map over its real residue pairs. Synthetic maps
  enter contact losses, thresholded contact metrics, and query-patch auxiliary
  diagnostics, but they are explicitly excluded from both contact AUPRC metrics.

## `{stage}_contact_loss`

The logged name remains the same for every `--loss-type`, but its meaning changes
substantially. Lower is better. Compare values only between runs using the same
loss type, `omega`, negative-map setting, and other relevant loss parameters.

### `--loss-type contact`

Plain BCE is averaged over known cells within each supervised map, then those map
losses are averaged equally:

```text
map_loss_m = mean_known_cells(BCEWithLogits(contact_logit, target))
contact_loss = mean_supervised_maps(map_loss_m)
```

### `--loss-type balanced_bce`

For each supervised map, positive and negative cells are averaged separately:

```text
map_loss_m = [omega * mean(BCE on positive cells)
              + (1 - omega) * mean(BCE on negative cells)]
             / sum(weights of classes that are present)
contact_loss = mean_supervised_maps(map_loss_m)
```

If a map contains only one class, the active class weight is renormalized so the
map's contribution does not shrink merely because the other class is absent.

### `--loss-type balanced_bce_dice`

The balanced BCE above plus a soft Dice loss. For a map with contact probability
`p`, target `y`, and `epsilon = 1e-6`:

```text
dice = (2 * sum(p * y) + epsilon) / (sum(p) + sum(y) + epsilon)
dice_loss = 1 - dice
contact_loss = balanced_bce + mean(dice_loss)
```

Dice loss is averaged only over maps containing at least one true contact.
All-negative maps still enter the balanced-BCE term.

### `--loss-type sparsity_11`, `sparsity_12`, `sparsity_21`, or `sparsity_22`

For each supervised map:

```text
predicted_sparsity = sum of contact probabilities / number of known cells
true_sparsity      = number of positive cells / number of known cells
soft_recall        = sum(p * y) / number of positive cells
recall_error       = 1 - soft_recall
```

For a map without positive cells, `recall_error` is defined as zero. The first
digit in the loss name selects the contact-density error, and the second selects
the recall-error term:

| Loss type | Per-map density term | Per-map recall term |
|---|---|---|
| `sparsity_11` | `abs(predicted_sparsity - true_sparsity)` | `recall_error` |
| `sparsity_12` | `abs(predicted_sparsity - true_sparsity)` | `recall_error^2` |
| `sparsity_21` | `(predicted_sparsity - true_sparsity)^2` | `recall_error` |
| `sparsity_22` | `(predicted_sparsity - true_sparsity)^2` | `recall_error^2` |

The two terms are added per map and supervised maps are averaged equally.

### Epoch-aggregation caveat for `contact_loss`

Within a batch, `contact_loss` is a mean over supervised maps. When Lightning
constructs the epoch value, however, the batch result is weighted by the full
protein-pair batch size, including examples without contact supervision. A batch
without any known contact cells logs contact loss zero. The resulting epoch value
is therefore a pair-batch-weighted average of batch map means, not a strict global
mean over supervised maps.

## Thresholded contact metrics

For these metrics, a known cell is predicted as a contact when
`contact_probability >= --contact-threshold`. Within each batch, all known cells
from all maps are pooled into one confusion matrix.

### `{stage}_contact_accuracy`

```text
(TP + TN) / (TP + TN + FP + FN)
```

Range: `[0, 1]`; higher is better. Contact maps are usually very sparse, so the
large number of true negatives can make this value look strong even when actual
contacts are missed.

### `{stage}_contact_precision`

```text
TP / (TP + FP)
```

Answers: “Of the residue pairs predicted as contacts, how many are true
contacts?” If the model predicts no positive cells, precision is defined as zero.
Range: `[0, 1]`; higher is better.

### `{stage}_contact_recall`

```text
TP / (TP + FN)
```

Answers: “Of the true residue contacts, how many did the model find?” If the
evaluated batch has no actual positive cells, recall is defined as zero. Range:
`[0, 1]`; higher is better.

### Epoch-aggregation caveat for threshold metrics

Contact accuracy, precision, and recall are cell-micro values **within each
batch**, but Lightning then averages the batch values using the number of protein
pairs as the weight. The epoch values are therefore not a single globally pooled
contact-cell confusion matrix. Precision and recall are likewise averages of
batch precision/recall rather than ratios of epoch-wide TP/FP/FN counts.

This behavior makes the threshold metrics useful as training diagnostics, but
the split-wide contact AUPRC metrics below are safer primary evaluation summaries.

## Contact AUPRC

### `val_contact_auprc_macro` and `test_contact_auprc_macro`

Exact AP is calculated independently for each true structural contact map after
removing its `-1` cells. The per-map AP values are then averaged equally, so a
small complex and a large complex have the same influence.

- A true structural map with known cells but no positive contacts contributes
  `AP = 0` and remains in the macro denominator.
- A map without any known cell is skipped because it cannot be evaluated.
- Synthetic negative maps are excluded.

Range: `[0, 1]`; higher is better. This is usually the most direct answer to:
“How well does the model rank contacts for a typical structural complex?”

### `val_contact_auprc_micro` and `test_contact_auprc_micro`

All known cells from all true structural maps are pooled into one precision-recall
curve. Larger or more completely resolved maps contribute more cells and therefore
have more influence. Synthetic negative maps are excluded.

To keep memory bounded, this metric places probabilities into 4,096 equal-width
bins before computing AP. Its probability resolution is approximately
`1 / 4096 = 0.000244`; predictions in the same bin are treated as tied. If the
evaluation population has no positive cells, the result is zero.

Range: `[0, 1]`; higher is better. This answers: “How well does the model rank
contacts across all evaluated residue pairs combined?”

Both custom contact AUPRC metrics keep additive distributed states, so their
results do not depend on batch boundaries.

## Query-patch-only losses and diagnostics

The final query-patch contact probability is the product of two branch
probabilities:

```text
final_contact_probability
    = sigmoid(query_gate_logit) * sigmoid(compatibility_contact_logit)
```

The query branch identifies possible interface regions; the compatibility branch
decides whether residue pairs in those regions are compatible. This is why both
branches receive separate diagnostics.

### `{stage}_optimization_loss`

The quantity returned by the shared step and actually optimized during training:

```text
auxiliary_loss = interface_mask_weighted_loss
               + compatibility_contact_loss

if the batch has known contact cells:
    optimization_loss = stage_loss + (1 - lambda) * auxiliary_loss
else:
    optimization_loss = stage_loss = interaction_loss
```

It is also computed and logged for validation/test for comparison, although no
optimization occurs there. Lower is better.

### `{stage}_contact_objective`

The complete contact-side objective before the outer `(1 - lambda)` task weight:

```text
contact_objective = contact_loss
                  + interface_mask_weighted_loss
                  + compatibility_contact_loss
```

This key is emitted only for batches with contact supervision. Its epoch mean is
weighted by the number of supervised maps. Lower is better.

### Interface target convention

An interface target is defined separately for each protein chain:

- a residue has a **known interface label** if any partner cell in its contact-map
  row/column is known;
- it is a **positive interface residue** if any such partner cell is a contact;
- otherwise it is a known negative/non-interface residue.

### `{stage}_interface_mask_loss`

Class-balanced BCE between the predicted chain-wise interface logits and the
interface targets. The balanced loss is computed for each chain using
`--interface-omega` (default `0.2`) as the positive-class weight and
`1 - interface_omega` as the negative-class weight, then the two chain losses
are averaged. This is independent of `--omega`, which controls the final and
compatibility contact-map BCEs. The epoch mean is weighted by the number of
supervised maps. Lower is better.

### `{stage}_interface_mask_weighted_loss`

```text
interface_mask_weighted_loss
    = query_interface_loss_weight * interface_mask_loss
```

This is the interface term that enters `auxiliary_loss`. The command-line option
is `--query-interface-loss-weight`; `--query-gate-loss-weight` remains an alias
for existing commands. Lower is better, but values should only be compared when
the configured weight matches.

### `{stage}_compatibility_contact_loss`

Class-balanced, per-map BCE between the standalone compatibility-branch logits
and known residue-pair contact targets. It uses the same `omega` behavior as
`balanced_bce`. The epoch mean is weighted by the number of supervised maps.
Lower is better.

### `{stage}_interface_mask_mean`

Mean predicted interface probability over every residue with a known interface
label, pooling both chains. The epoch mean is weighted by the number of such
residues. This is a calibration/activity diagnostic, not a quality metric: there
is no universally better direction.

### `{stage}_interface_mask_positive_mean`

Mean predicted interface probability on true interface residues. The epoch mean
is weighted by the number of positive interface residues. Higher is generally
better, ideally approaching `1`, but it should be interpreted together with the
negative mean. The metric is absent when no positive interface residue is logged.

### `{stage}_interface_mask_negative_mean`

Mean predicted interface probability on known non-interface residues. The epoch
mean is weighted by the number of negative interface residues. Lower is generally
better, ideally approaching `0`. The metric is absent when no negative interface
residue is logged.

Useful separation means
`interface_mask_positive_mean > interface_mask_negative_mean`; a high overall
`interface_mask_mean` alone does not establish good discrimination.

### `{stage}_query_gate_scale`

The learned positive multiplier in

```text
query_gate_logit = query_gate_scale * nonnegative_query_evidence
                 + query_gate_bias
query_gate_scale = softplus(raw_scale) > 0
```

A larger scale makes gate probability respond more strongly to query evidence.
It is a learned parameter diagnostic, not a quality score, so neither direction
is inherently better. Because the parameter changes during training, an epoch
value is a supervised-map-weighted average of the values seen across batches,
not necessarily the exact end-of-epoch parameter value.

### `{stage}_query_gate_bias`

The learned gate logit when query evidence is zero. `sigmoid(query_gate_bias)` is
therefore the baseline gate probability at zero evidence. A more positive bias
makes the gate more permissive; a more negative bias makes it more restrictive.
There is no universally better direction. As with gate scale, the epoch value is
a supervised-map-weighted average across logged batches.

## Performance diagnostic

### `{stage}_contact_metrics_time_s`

Mean wall-clock seconds per batch spent on:

- main contact-loss and threshold-metric computation; and
- query-patch auxiliary contact/interface computation, when applicable.

GPU synchronization surrounds the timed region so queued CUDA work is included.
The model forward pass and interaction-loss calculation are outside the timer.
Lightning uses `batch_size=1`, so the epoch value is an unweighted mean over
batches. Lower means faster for a comparable batch/map-size and hardware setup,
but this is not a model-quality metric.

## Logger bookkeeping and output locations

All ordinary metrics are sent to both configured loggers:

- Lightning CSV under
  `<output-dir>/<output-subdir>/metrics/version_<n>/metrics.csv`;
- W&B project `master`.

The CSV is sparse: a row can contain a value for only some metric columns. Its
additional columns are:

| Field | Meaning |
|---|---|
| `epoch` | Zero-based Lightning epoch associated with the log row. |
| `step` | Optimizer/global step associated with the log row. |

W&B can additionally show `epoch`, `trainer/global_step`, `_step`, `_timestamp`,
`_runtime`, `_wandb`, and system-monitoring fields. These are logger/run metadata,
not project-defined model metrics. Saved hyperparameters are configuration values,
not metrics.

`test_predictions.tsv` contains per-example labels, predicted labels, and
probabilities; it is an output table rather than an aggregate metric.

## Aggregation summary

| Metric family | Epoch/split aggregation |
|---|---|
| Interaction accuracy | Exact pair mean via protein-pair batch-size weighting. |
| Interaction loss | Protein-pair batch-size-weighted mean of batch losses. |
| Interaction AUPRC | Exact state over the complete validation/test split. |
| Contact loss | Pair-batch-weighted mean of per-batch supervised-map means. |
| Contact accuracy/precision/recall | Known-cell pooling inside each batch, followed by protein-pair batch-size weighting. Not a global cell confusion matrix. |
| Contact macro AUPRC | Exact AP per true structural map, then equal-map mean over the split. |
| Contact micro AUPRC | 4,096-bin approximate AP over all known cells from true structural maps. |
| Query auxiliary losses | Supervised-map-count-weighted epoch means. |
| Interface probability means | Known/positive/negative-residue-count-weighted epoch means. |
| Contact metric time | Unweighted mean seconds per batch. |

## Which metrics should be primary?

- For protein-pair interaction quality, use `val_interaction_auprc` during model
  development and `test_interaction_auprc` for final reporting. Keep accuracy as
  an intuitive thresholded companion.
- For contact prediction, report both macro and micro contact AUPRC. Macro gives
  every structural complex equal influence; micro describes the pooled
  residue-pair population.
- Use contact precision and recall to understand behavior at the configured
  threshold, but remember their epoch aggregation and sensitivity to that
  threshold.
- For D-SCRIPT and query patch, use `best_val_contact_auprc_macro` for checkpoint
  and sweep selection. Do not compare `contact_loss` across different loss types.
- For query patch, use `optimization_loss` to understand the full optimized
  objective and the interface/compatibility metrics to diagnose which branch is
  improving or collapsing.

## Legacy query-gate metrics in older runs

Older W&B/CSV runs may contain the following names. They came from a previous
query-patch objective that directly supervised the residue-pair query gate against
contact-map cells. The current implementation instead supervises chain-wise
interface masks and the compatibility branch.

As everywhere else in this document, `{stage}` expands to `train`, `val`, and
`test`.

| Legacy metric | Meaning |
|---|---|
| `{stage}_query_gate_loss` | Class-balanced per-map BCE between query-gate logits and known residue-pair contact targets. Lower is better. |
| `{stage}_query_gate_weighted_loss` | `query_gate_loss_weight * query_gate_loss`; the former auxiliary term. Lower is better for the same weight. |
| `{stage}_query_gate_mean` | Mean gate probability over all known contact-map cells. Diagnostic only. |
| `{stage}_query_gate_positive_mean` | Mean gate probability over true contact cells. Higher was generally desirable. |
| `{stage}_query_gate_negative_mean` | Mean gate probability over known non-contact cells. Lower was generally desirable. |
| `{stage}_query_gate_open_fraction` | Fraction of known cells with gate probability at least `0.9`. Diagnostic for an always-open gate. |
| `{stage}_query_gate_closed_fraction` | Fraction of known cells with gate probability at most `0.1`. Diagnostic for an always-closed gate. |

The still-current `{stage}_query_gate_scale` and `{stage}_query_gate_bias` retain
the meanings documented above. Current replacements for the removed quality/loss
keys are `{stage}_interface_mask_*` and
`{stage}_compatibility_contact_loss`.

## Implementation pointers

- Metric logging and loss relationships:
  [`execute_pipeline.py`](../scripts/training_pipeline/execute_pipeline.py)
- Contact AUPRC implementations:
  [`metrics.py`](../scripts/training_pipeline/metrics.py)
- Reusable balanced-BCE implementation:
  [`losses.py`](../scripts/training_pipeline/losses.py)
- Query-gate and compatibility equations:
  [`query_patch.py`](../scripts/models/query_patch.py)
- Metric and plot behavior tests:
  [`test_training_metrics.py`](../tests/test_training_metrics.py) and
  [`test_training_plots.py`](../tests/test_training_plots.py)
