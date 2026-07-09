## Code style

- Prefer compact, direct code, but do not make code terse at the expense of readability.
- Optimize for human readability: clear names, simple control flow, and code that can be read top-to-bottom.
- Avoid scattering logic across many small helper functions. Keep related logic together when it is easier to understand in one place.
- Extract a function only when it meaningfully improves clarity, removes real duplication, isolates a distinct concept, or makes testing easier.
- Do not add abstractions just to make code look cleaner; favor straightforward, maintainable code.

## Commenting and explainability

- Prefer code that is immediately understandable to a new reader, especially for ML, data loading, loss functions, tensor transformations, and model architecture.
- Add explanatory comments and docstrings for every non-trivial block of code.
- Comments should explain why the code exists, what assumptions it relies on, expected input/output shapes, label encodings, masking/padding behavior, edge cases, and how losses or metrics are computed.
- For ML code, explicitly comment tensor shapes, value ranges, thresholds, unknown/missing-value conventions, and how batches are padded or masked.
- Use highly visible section comments to separate major parts of a file, for example general/shared code, baseline-only code, model-specific code, and training orchestration.
- Avoid comments that only restate obvious syntax. Prefer comments above logical blocks that explain intent and data flow.
- When adding a new model or pipeline component, include a short top-level overview explaining how data moves through it and how it is trained/evaluated.