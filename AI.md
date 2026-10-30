# AI Guidelines

Applies to any AI assistant working here (Claude, Codex, local models, ...).

## Start of every session

1. Read `docs/overview.md`: goals, status, decisions so far.
2. Read the most recent file in `docs/sessions/`.
3. Ask the user to clarify the session goal if it is not clear.

## During a session

- `docs/overview.md` is the source of truth for project state.
- Architecture changes → `docs/design.md`. New or changed inputs → `docs/data.md`.
- **Never modify anything in `~/projects/Sc-x2`.** It is a read-only input. All
  generated files go in `data/` here (gitignored).
- **Strain/line naming:** in anything a human reads, use AP2 gene names
  (e.g. `AP2XI-2/AP2XII-1`), never shorthand like `2KD`/`Duo`/`TripleKD`.
  Same rules as Sc-x2 (`~/projects/Sc-x2/docs/naming-conventions.md`).
  Pre-existing file/folder names are left as they are.
- Barcodes: the h5ad/CSV barcodes carry a sample suffix (`AAAC...-1-0`); the
  BAM `CB` tag does not (`AAAC...-1`). Always join through an explicit
  dataset + condition → BAM mapping (see `docs/data.md`), never by string match
  on the barcode alone.
- Serving: igv.js needs HTTP Range requests. `python -m http.server` does
  not support them; use the project's server (see `docs/design.md`).

## End of every session

1. Write `docs/sessions/YYYY-MM-DD_<topic>.md` (template below).
2. Update `docs/overview.md`.
3. `git add` + `git commit` tracked changes.
4. Tell the user what was written and what is flagged for next time.

### Session file template

```markdown
# Session — YYYY-MM-DD: <topic>

## Goal
## What was done
## Results / findings
## Decisions made
## Next session
```
