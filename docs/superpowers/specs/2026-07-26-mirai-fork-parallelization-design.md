# Design: fork-on-Unix / mirai-on-Windows parallelization for `PTE_bootstrap_inference`

## Context

`PTE_bootstrap_inference()` in `PTE/R/bootstrap_inference.R` is the only place in the
package that runs parallel work. It has two `foreach %dopar%` loops (the main
`B`-sample bootstrap loop, and the optional leave-one-out BCA loop), both run on a
single PSOCK cluster created via `makeCluster(num_cores)` / `registerDoParallel()`.

`PTE/CHANGELOG` already lists, under "version 2.0": *"Parallelization improved using
mirai (for windows) and forking (for *nix)"* — this design implements that line.

No other file in the package (`PTE.R` etc.) calls parallel primitives directly; `PTE.R`
only carries the roxygen `@import` tag.

## Goals

- On Unix/Linux/macOS, use a `FORK` cluster so worker processes inherit already-built
  objects (`Xy`, `cutoff_obj`, model-building closures, etc.) via copy-on-write instead
  of paying PSOCK serialization cost.
- On Windows (no fork available), use `mirai::make_cluster()`, which is a drop-in
  cluster object compatible with `doParallel::registerDoParallel()` /
  `parallel::stopCluster()`.
- Keep the existing `%dopar%` loop bodies, error-handling semantics
  (`.errorhandling = "stop"`, the default), and `verbose` dot-printing untouched.
- Fix two incidental bugs uncovered while reading this code:
  1. The current `num_cores` auto-detect reads `Sys.getenv('NUMBER_OF_PROCESSORS')`,
     a Windows-only env var — on Unix this is empty, so the computed default is `NA`.
  2. `stopCluster(cluster)` is only called inside the `if (run_bca_bootstrap)` branch,
     so the cluster (and its worker processes) leaks whenever
     `run_bca_bootstrap = FALSE` (the default) or when an error is thrown mid-run.

## Design

### 1. Core-count default, detected once at package load

`num_cores` stays a normal argument to `PTE_bootstrap_inference()` — callers can still
pass an explicit value to override. When left `NULL`, it now resolves to
`max(cores - 1, 1)` where `cores` is detected once at package load time (not
per-call), via R's standard `mc.cores` option so it composes with any prior user
configuration:

`PTE/R/zzz.R` gains an `.onLoad`:

```r
.onLoad = function(libname, pkgname){
    if (is.null(getOption("mc.cores"))){
        options(mc.cores = max(parallel::detectCores() - 1, 1))
    }
}
```

(guarded so it never clobbers an `mc.cores` option the user set before loading PTE).

In `PTE_bootstrap_inference()`, the existing block:

```r
if (is.null(num_cores)){
	num_cores = as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
	num_cores = max(num_cores - 1, 1)
}
```

becomes:

```r
if (is.null(num_cores)){
	num_cores = getOption("mc.cores")
}
```

The `@param num_cores` roxygen text is updated to describe the new default
(`max(cores - 1, 1)` from cores detected at package load via
`options(mc.cores = ...)`, overridable per-call).

### 2. OS-conditional cluster creation

Replacing the current `cluster = makeCluster(num_cores); registerDoParallel(cluster)`:

```r
cluster = if (.Platform$OS.type == "windows"){
	if (!requireNamespace("mirai", quietly = TRUE)){
		stop("The \"mirai\" package is required to parallelize on Windows. Please install it via install.packages(\"mirai\").")
	}
	mirai::make_cluster(num_cores)
} else {
	makeCluster(num_cores, type = "FORK")
}
registerDoParallel(cluster)
on.exit(stopCluster(cluster), add = TRUE)
```

Placed at the same point in the function (right before the first `%dopar%` loop), so
both the bootstrap loop and (if enabled) the BCA loop share the one cluster, same as
today. The `requireNamespace` check only runs (and only matters) on the Windows
branch — Unix/Linux/macOS callers never touch `mirai` at all.

### 3. Cluster teardown

Delete the existing explicit `stopCluster(cluster)` call inside the
`if (run_bca_bootstrap)` branch. The `on.exit(stopCluster(cluster), add = TRUE)` added
in step 2 now handles cleanup unconditionally — whether or not BCA runs, and even if
an error is thrown partway through — closing the leak.

### 4. Loop bodies

Both existing `%dopar%` blocks are otherwise untouched. They keep running through
`foreach`/`doParallel` exactly as before; only the underlying cluster object changes
per OS.

## Dependencies

- Add `mirai` to `DESCRIPTION` `Suggests:`, not `Imports:` — it's only ever touched on
  the Windows branch, so Unix/Linux/macOS installs of PTE never need it. Guarded at
  runtime with `requireNamespace("mirai", quietly = TRUE)`, erroring with an
  informative message (naming the package and the install command) if a Windows user
  hasn't installed it.
- Call it fully qualified as `mirai::make_cluster()` — no NAMESPACE entry needed for
  it, since it's used for exactly one call (unlike `parallel`/`foreach`/`doParallel`,
  which are blanket-`@import`ed and used unqualified throughout this function).
- `foreach`, `doParallel`, and `parallel` stay exactly as they are today (all remain
  required `Imports:`, used on every platform).

## Error handling

Unchanged: `%dopar%` keeps its default `.errorhandling = "stop"`, so a failure in any
bootstrap iteration aborts the whole call exactly as it does today. The only behavior
change is that cluster cleanup is now guaranteed via `on.exit`.

## Testing

No test suite currently exists in the package. Verification plan:

- On this Linux machine: build/install the package, run `PTE_bootstrap_inference()`
  with a small `B` against the bundled `continuous_example` / `incidence_example` /
  `survival_example` datasets, with `verbose = TRUE`, both with and without
  `run_bca_bootstrap = TRUE`. Confirm results look sane (consistent with pre-change
  runs) and that no leftover worker processes remain after the call returns (check via
  `ps`), including after an artificially forced error mid-run.
- The Windows/`mirai` path cannot be exercised on this machine. This relies on
  `mirai::make_cluster()`'s documented compatibility with
  `registerDoParallel()`/`stopCluster()` as a drop-in cluster type (confirmed via
  mirai's own vignette / R 4.5 official "MIRAI" cluster type documentation). Flagged
  here as an assumption, not verified end-to-end on Windows. On this Linux machine I
  can still verify the `requireNamespace("mirai", ...)` guard itself is unreachable
  dead code on the Unix branch (i.e. Unix runs never require `mirai` to be installed),
  but not the Windows success path or its error message.

## Out of scope

- The `num_cores` parameter's existence/signature otherwise. A stale `CHANGELOG` line
  under "version 1.6" claims it was "removed" in favor of `options(mc.cores = ...)`,
  but it is still present and used in the current code — this design keeps it as an
  override, and now also happens to make the underlying default use the same
  `mc.cores` option, resolving that discrepancy without removing the argument.
- Adding a test suite.
- Any change to `PTE.R` beyond nothing — it doesn't call these functions directly.
