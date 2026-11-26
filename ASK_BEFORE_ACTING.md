# Ask-before-acting rules (for Claude)

These rules apply to every session in this repository. They exist because intermediate
files are huge, cluster time is expensive, and physics results must be checked by Maria
after every step before anything is deleted.

## Always ask first, and wait for an explicit OK
- **Before each workflow step** from `Recipe.md` (McMule integration, evtgen, lhe->root,
  analysis, split, TGEANT, Coral+Phast, merge, ...): say which step, which config/inputs,
  how many events/jobs, expected runtime and disk usage. Then wait.
- **After each step**: show the result (logs, chi^2, event counts, file sizes, errors)
  and ask whether to continue. Never chain steps on my own.
- **Before deleting, moving or overwriting** any data or output file (LHE, ROOT, .mcmule,
  logs) - even if it looks unnecessary. List exactly what would be removed and how big it is.
- **Before submitting or killing cluster jobs** (`qsub`, `qdel`, submit_*.sh scripts).
- **Before changing an input `.conf`** that a running or queued job might read.
- **Before changing physics code** (userfiles, McMule, TGEANT, Phast user code) or
  rebuilding it.

## Never
- Never reuse an event-generation seed (see `inputs/used_seeds_evtgen.txt`).
- Never delete lines from `inputs/used_seeds_evtgen.txt`.
- Never delete large outputs before Maria confirmed the important information is
  stored in the reduced trees.
- Never push to a remote without being asked.

## Fine without asking
- Committing the AMBER_RadCor repo locally at regular points (after a finished step / batch of edits),
  with a clear message. Never commit large data, binaries or secrets. (Maria, 2026-11-26)
- Reading files, logs, `qstat`, `df -h`, `du`, `head`/`tail` of outputs.
- Writing scratch/test files in the session scratchpad.

The technical enforcement lives in `~/.claude/settings.json` (`permissions.ask`):
`rm`, `mv`, `qsub`, `qdel`, `find -delete` and the submit scripts always trigger a prompt.
