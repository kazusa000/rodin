# L15 graded FEM campaign

The source branch is `codex/p03-fem-domain`, based on the compiler-only Carlos
reference `3ae2c041`. The selected initial sphere has radius one and the full
outer cube has edge 30 (`outer-radius=15`). The geometry and FEM optimization
method are inherited; no BEM, new thickness correction or seed selection is used.

The reference spacing remains exactly `h0=1/6`, matching n13/L2. The same
minimum MMG size `h0/10`, far size `10h0`, Hausdorff tolerance `0.1h0^2`,
regularization `4h0` and step `0.1h0` remain fixed. A new opt-in coordinate
grid is uniform through coordinate radius 2 and then grows by 1.3, capped
at the far size. The identical monotone map on all axes preserves chamber
cut planes. It avoids allocating an n91 uniform cube before adaptation.
This is a seed grading policy, not a guarantee that every adapted edge has
exactly the requested length. Existing distance-weighted MMG adaptation
targets the same near size around the evolving interface and coarse far field.

Use the preset `reproduction/fem-domain-L15.json`. The native diagnostic checks
the original near coordinates, unchanged MMG size targets and monotone finite
graded intervals. The CPU-only runner verifies the unchanged n13/L2 one-update
FEM values against the frozen baseline, and validates L15 geometry before
starting the optimization. All stages execute on 28h in separate output folders.

The actual optimizer requests 21 evaluated states: initial state plus 20 real
level-set/MMG updates. Completed-state CSV, transport increments, reconstruction
timings, commands, source/submodule/binary hashes and actual exit codes are saved.
MMG retries are explicitly disabled to avoid silently halving the mesh scale.
An exception stops the campaign and retains inputs; no tolerance is relaxed.
The existing failed-input exports remain enabled. Full fields stay on 28h.

The launcher `scripts/run-fem-domain.py` requires an exact clean pushed commit,
new build/output directories and the original 28h baseline build and campaign
records. It inherits only the verified dependency paths/configuration, not old
binaries. CPU execution uses four OpenMP threads and one BLAS thread by default.
