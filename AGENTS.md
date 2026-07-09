## Code style

- Prefer compact, direct code, but do not make code terse at the expense of readability.
- Optimize for human readability: clear names, simple control flow, and code that can be read top-to-bottom.
- Avoid scattering logic across many small helper functions. Keep related logic together when it is easier to understand in one place.
- Extract a function only when it meaningfully improves clarity, removes real duplication, isolates a distinct concept, or makes testing easier.
- Do not add abstractions just to make code look cleaner; favor straightforward, maintainable code.