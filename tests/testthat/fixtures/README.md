Shared test fixtures live here.

- `data/` is for reusable tabular or serialized datasets used by multiple tests.
- `matrices/` is for reusable sparse or dense matrix fixtures.
- `models/` is for reusable model specs or saved fit objects when inline setup becomes noisy.

Keep fixtures small and deterministic. If a fixture is only used once, prefer building it inline in the test file.
