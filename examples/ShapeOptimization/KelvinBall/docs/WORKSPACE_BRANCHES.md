# Project 03 baseline and research branches

The workspace organization of 14 October 2026 selects the original Carlos
MMG/SWIFT program as its reproducible comparison baseline. The baseline is
not the coupled method and must not be inferred from the newest branch name.

## Fixed references

- Original upstream source: `79718dc7d018cf41c727dc2df2634628b42612c4`.
- Compiler-only baseline: `3ae2c0414ffebf11f5814a46510ddf4d8375ba0a`, tag
  `p03/carlos-original-20261011`. The entire source difference is one `typename`.
- Integrated source before this organization: `c1fe06404bf693e1d8900b91c50febc0587daba2`,
  tag `p03/bem-integration-20261011`.
- Actually executed native/migration source: `d625db418e4cf7359a044cc1f60913a2e54fdab4`;
  see the preserved migration validation, not a later documentation commit.

## Branch responsibilities

| Branch | Start | Responsibility and status |
| --- | --- | --- |
| `codex/p03-fem-domain` | compiler-only baseline | FEM outer-domain and mesh/cost comparisons; no new campaign executed by this organization |
| `codex/p03-bem-coupled` | integrated source | six BEM wall traces, internal FEM proxy direction, measured BEM acceptance and controlled MMG/SWIFT comparisons |
| `codex/p03-bem-gradient` | integrated source | reserved for exact BEM derivative research; it currently still contains the proxy direction and has no separate checkout |

MMG and SWIFT remain reconstruction choices within the same workflow. The
current native coupled solver is cubic O24; reading old T12 histories is not
an implemented T12 native optimizer. Changing thickness algorithms, migrating
symmetry or adding a long loop requires its own explicit experiment and evidence.

The baseline tag freezes the original n13/L2/thickness-off program. A large
outer domain is a research variation, not a retroactive change to that reference.
Likewise, the FEM proxy direction's known finite-difference discrepancy and
the original-panel open seams remain documented limitations.

## Local paths and provenance

Within the internship workspace, the detached baseline is at
`tmp/project-03/baseline/rodin-carlos-20261011`; active FEM and coupled source
are at `tmp/project-03/branches/{fem-domain,bem-coupled}/rodin`.
Builds, environments, generated fields and results remain outside these sources.
The unified entry remains `scripts/kelvin-bem.py`; `compare` validates a bounded
paired update, not a long optimization.

The workspace index `notes/project-03/workspace.md` defines the current entry
and retained historical result roles. Exact local moves and checksums are in
`metadata/project-03-organization-20261014.json`. All old branch references and
formal result paths are retained. Old absolute-path launch records are evidence;
restore their original paths before replaying them.

No native build, optimization, remote deployment, numerical parameter change,
or source merge is performed by this organization. The shared original Rodin
checkout used by other projects is unchanged.
