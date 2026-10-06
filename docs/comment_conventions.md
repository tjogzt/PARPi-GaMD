# Inline-comment conventions (adopted 2026-10-06)

Adopted from the four-team code-audit review (seat c, Q3); complements the
standardised script headers (Purpose/Created/Depends/Run). Applies to
`code/`, `scripts/` and `common/`.

1. **Units.** Every physical quantity, threshold or coordinate mentioned in a
   comment carries its unit (e.g. `# cutoff = 5.0 Å`, `# temperature = 310 K`).
2. **Function notes.** A custom function or key block gets a one-line comment
   stating inputs, outputs and side effects.
3. **No conversational tone.** No `TODO: fix later` / `I think …`; use
   objective descriptions (`# NOTE: convergence threshold may need
   adjustment for membrane systems.`).
4. **Magic numbers.** Hard-coded numbers must state their source or meaning
   (e.g. `# 0.002 ps timestep (2 fs)`, `# 3-sigma threshold per <ref>`).
5. **Density.** Roughly one comment per 5–10 lines of logic; complex
   algorithms (DCCM matrix maths, GaMD reweighting) are explained block by
   block.

**Spot-check (2026-10-06, mechanical):** zero `TODO`/`FIXME`/conversational
phrases and zero non-ASCII (Chinese) comments across `code/`, `scripts/` and
`common/` — compliant at adoption time.
