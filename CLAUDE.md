# AMBER_RadCor

Radiative corrections for mu-p scattering at AMBER (McMule -> TGEANT -> Coral/Phast).
The step-by-step workflow is in `Recipe.md`; planned improvements in `TODO_future.md`;
current tasks + results in the newest file in `meetings/` (now `meetings/2026-KW48-49_ToDo.md`; keep status,
results and result paths updated there). Findings for the thesis: one file per week in `Erkenntnisse/` (now `Erkenntnisse/2026-KW48_Erkenntnisse.md`);
add every important finding there while working. New week -> new ToDo file in `meetings/` and new Erkenntnisse file.

@ASK_BEFORE_ACTING.md

## Key locations
- Scripts: `commands/` (McMule evtgen: `commands/mcmule_eventgen/`), configs: `inputs/`
- Large data / outputs: `/nfs/momos/user/mleibelt` (check space with `df -h /nfs/momos/user/mleibelt`)
- Cluster: SGE (`qsub`, `qstat`), queues e.g. `long_epyc`
- Used event-generation seeds: `inputs/used_seeds_evtgen.txt`

## Notes
- Claude's shell does not read `~/.bashrc` (it returns early without a tty). Run repo scripts that use
  `module load gcc/...` in a login shell: `bash -lc 'bash <script>'` — `~/.bash_profile` loads spack
  there (since 2026-11-26).
- McMule CU generation only needs the `.mcmule` integration files; the weighted N run is optional.
- For NLO0 and R integration use the `mcmule-event-generator` build, for other pieces `mcmule-devel`.
- Keep `TODO_future.md` up to date and point out when something there becomes relevant.
