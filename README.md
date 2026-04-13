# ISA-NUTS_potentials

## Resource planning from benchmarks

This workflow uses Snakemake benchmark TSV files to estimate per-rule resource
requests.

Each benchmark file stores one summary row for one benchmark path. If the same
job is executed again, that TSV is rewritten. Different wildcard combinations
produce different benchmark files.

For each rule, resource inference follows this order:

1. Use the exact benchmark for that job if it exists.
2. Otherwise use the worst observed benchmark available for that rule.
3. Otherwise use a conservative static default defined in the Snakefile.

The workflow derives:

1. `threads` from `cpu_time / s`, with a small safety margin.
2. `mem_mb` from `max_rss`, with a small safety margin.

## Important distinction: per-job vs global memory

Declaring `resources: mem_mb=...` inside rules does not by itself impose a
global memory cap in local execution. Those values become active for scheduling
only when Snakemake is started with a global resource budget.

Examples:

```bash
pixi run snakemake all --cores 32
```

This limits total CPU usage through `threads`, but does not use `mem_mb` to
limit overall concurrency.

```bash
pixi run snakemake all --cores 32 --resources mem_mb=230000
```

This limits both:

1. Total CPU usage to 32 cores.
2. Total scheduled memory to 230000 MB across concurrent jobs.

For a machine with 252 GB RAM, `mem_mb=230000` is a reasonable starting point
that leaves headroom for the OS and non-workflow processes.

## Note on enforcement

In local execution, Snakemake resources are used for scheduling, not as hard OS
limits on the spawned processes. In other words, `mem_mb` helps avoid launching
too many memory-hungry jobs at once, but it does not sandbox each process.