# Contributing to `mfrmr`

Thanks for helping improve `mfrmr`.

## Before you start

- Search existing issues and pull requests first.
- If behavior changes are large, open an issue before implementation.
- Keep changes focused and reviewable.

## Development setup

```r
# from package root
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::document()
devtools::test()
devtools::check(args = c("--no-manual"), document = FALSE)
```

## Coding guidelines

- Prefer readable, explicit code over compact but opaque code.
- Keep public API names descriptive.
- For user-facing behavior changes, update docs and examples in the same PR.
- Use base R plotting defaults in this package unless there is a strong reason not to.

## Testing expectations

- Add or update tests for every user-visible change.
- Keep tests deterministic (fixed seeds, no network calls).
- Avoid writing plot files during tests (use `draw = FALSE` unless plotting is under test).

## Documentation expectations

- Update roxygen comments when function arguments/returns change.
- Run `devtools::document()` before committing.
- If workflow changes, update README and/or vignette accordingly.

## Pull request checklist

- [ ] Code and docs updated together
- [ ] `devtools::test()` passes
- [ ] `devtools::check(args = c("--no-manual"), document = FALSE)` passes
- [ ] NEWS entry updated when behavior changes materially

## Reporting bugs

Please include:

- minimal reproducible example,
- session info (`sessionInfo()`),
- expected vs observed behavior,
- relevant data schema (without sensitive content).
