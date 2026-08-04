# Model Implementation Suggestions

This document records model and pipeline ideas that have been discussed but are
not approved for implementation. A suggestion remains documentation only until
its status is explicitly changed to `Approved` by the project owner.

## Implementation Control

Allowed statuses are:

- `Proposed`: recorded for review; do not implement.
- `Approved`: implementation scope has been explicitly accepted.
- `Implemented`: code and tests have been completed and linked here.
- `Rejected`: retained as decision history but should not be implemented.

When approving a suggestion, also record the intended files, compatibility
requirements, experiment configuration, and required tests. Implementations
must not silently include other proposals from this document.

## Baseline Invariants

- Preserve the existing public API and behavior of `scripts/models/dscript.py`.
- Preserve D-SCRIPT's sequence-padding and structural-annotation mask semantics.
- Put experimental architectures in separate model modules.
- Reuse baseline behavior through composition or model-specific overrides rather
  than changing the baseline to accommodate an experiment.
- Require a focused D-SCRIPT regression check for any approved change that must
  touch shared training code.
- Record optional improvements here before implementing them.

## Current Implemented Scope

The `query_patch` experiment currently consists of:

- a separate learned-query contact module in `scripts/models/query_patch.py`;
- query-specific contact-logit generation while reusing the established
  D-SCRIPT projection and interaction heads;
- model selection and query hyperparameters in the training pipeline;
- the `jobs/training_pipeline/esm2_query_patch.sh` Slurm entry point; and
- focused tests for the new model's shapes, gradients, symmetry, masking, and
  inherited interaction-head contract.

The D-SCRIPT model implementation remains unchanged.

## Proposal Register

| ID | Proposal | Status | Baseline impact |
| --- | --- | --- | --- |
| BASE-01 | Mask padded pair features inside the D-SCRIPT contact CNN | Proposed | Changes D-SCRIPT outputs and BatchNorm statistics |
| PIPE-01 | Pool interaction evidence over all sequence residues instead of structurally known cells | Proposed | Changes legacy training and inference semantics |
| ARCH-01 | Extract projection and interaction heads into a shared base/composed component | Proposed | Broad refactor; intended to preserve behavior |
| QP-01 | Initialize queries dynamically from high-scoring residues | Proposed | New model only |
| QP-02 | Enhance queries with a coarse per-residue interface estimate | Proposed | New model only |
| QP-03 | Update residue features and queries reciprocally | Proposed | New model only |
| QP-04 | Add query specialization or diversity supervision | Proposed | New model and loss only |
| CONTACT-01 | Predict sparse top-k residue-pair edges before constructing a dense map | Proposed | Alternative contact model |
| CONTACT-02 | Use axial or pair-grid Transformer refinement | Proposed | Alternative contact model |
| HEAD-01 | Predict interaction directly from decoded queries in addition to map pooling | Proposed | Alternative interaction head |

## Deferred Baseline and Pipeline Changes

### BASE-01: Padding-aware D-SCRIPT contact CNN

**Suggestion:** Zero padded residue-pair features before the D-SCRIPT 2D
convolution and after intermediate convolution or normalization blocks.

**Motivation:** Prevent padded values from affecting valid cells near sequence
boundaries through the local convolution.

**Decision required:** Whether this is a baseline bug fix or a new D-SCRIPT
variant. It changes existing predictions and BatchNorm statistics, so it should
not be bundled with an experimental architecture.

**Required validation:** Frozen-input comparison, variable-padding tests, old
checkpoint loading, and a separate baseline training run.

### PIPE-01: Inference-safe interaction pooling mask

**Suggestion:** Make it optional to ignore `interaction_pair_mask` and pool over
all non-padding sequence residue pairs.

**Motivation:** `interaction_pair_mask` is derived from structural annotation
availability, which is unavailable for a new sequence-only inference pair.

**Decision required:** Preserve legacy behavior, replace it, or expose an
explicit CLI option supporting both. Comparisons must use the same policy for
all architectures.

**Required validation:** Legacy-output regression, sequence-only inference, and
an ablation comparing both pooling policies.

### ARCH-01: Shared model components

**Suggestion:** Extract the projection and map-to-interaction heads into reusable
components instead of inheriting the complete D-SCRIPT model.

**Motivation:** Give experimental contact generators an explicit interface while
keeping the baseline class independent.

**Decision required:** Whether the reduction in coupling justifies a broad
refactor. State-dict and Lightning-checkpoint compatibility must be specified
before approval.

## Deferred Query-Patch Variants

### QP-01: Dynamic residue-derived queries

Use a shared proposal scorer on both proteins and initialize queries from
high-scoring residue features rather than only from fixed learned parameters.
This would be closer to QORT's feature-conditioned query proposals. Decisions
are needed for differentiable selection, A/B swap symmetry, and the number of
queries selected from each protein.

### QP-02: Coarse interface-enhanced queries

Predict one coarse interface/contactness value per residue and incorporate that
signal into the query inputs before decoding. This avoids requiring a full
residue-residue contact map before the final contact prediction. Decisions are
needed for supervision, behavior on examples without contact maps, and whether
the coarse head is shared between proteins.

### QP-03: Reciprocal residue-query updates

Add a residue-to-query update followed by a query-to-residue update so both
representations are refined. This is closer to QORT's feature/query
co-optimization but increases memory and compute. The update order, number of
layers, residual connections, and symmetry constraints require approval.

### QP-04: Query specialization losses

Add diversity, coverage, or assignment losses to discourage multiple queries
from selecting the same interface patch. Candidate approaches include mask
overlap penalties and matching supervised contact components to queries. The
loss definition and weight require approval before implementation.

## Deferred Alternative Contact and Interaction Models

### CONTACT-01: Sparse bipartite edge prediction

Score a restricted set of candidate residue pairs and construct the dense map
only for output or supervision. Candidate selection must retain enough positive
contacts and remain trainable without structural input at inference.

### CONTACT-02: Axial or pair-grid Transformer refinement

Create a coarse pair representation and alternate attention along the two
residue axes. This models dependencies directly on the contact grid but can be
substantially more expensive than the low-rank query-patch decoder.

### HEAD-01: Query-level interaction prediction

Pool decoded query states or patch strengths directly into an auxiliary PPI
score. This could complement the established normal/max map heads, but it would
change the intended constraint that interaction evidence is extracted from the
predicted contact map.

## Approval Record

When a proposal is approved or rejected, add a dated entry here with the chosen
scope and rationale. There are currently no approved deferred proposals.
