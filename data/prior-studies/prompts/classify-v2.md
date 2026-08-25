# Heterogeneous effects code classification — prompt v2

Revision of v1 following hand adjudication of 100 file-chunks by the author.

What the adjudication showed. Precision was 91.5% for gpt-oss-120b and 76.7%
for Llama 3.3 70B, so Llama is dropped. The remaining problem is recall: the
files gpt-oss missed were almost all cases where the heterogeneity is visible
only in the shape of the whole analysis, not in any single line — a loop
fitting one model across subgroups, a `margins, at()` sweep over a moderator,
an interaction variable built in one place and used in another, a `by=` smooth
in a GAM. Meanwhile the over-flagging model fired on data cleaning, plotting,
and file loading, which are local surface patterns.

So v2 reframes the task from line-hunting to reading the analysis as a whole
and asking what it was written to estimate. Evidence is still required, but as
support for a judgment rather than as the basis of one. v1 also contained a
rule treating variable names joining two names with `x` as interactions; that
rule caused `gen linnovation_x_oil3 = .` to be read as an interaction and is
removed.

## SYSTEM

You are reading statistical code from a political science replication archive.

Read the whole file first and work out what the analysis is trying to estimate.
Then answer one question: **does this file estimate heterogeneous effects?**
That is, does it fit or compute at least one thing that lets the effect of some
variable differ across units, subgroups, or contexts — or that compares such an
effect between groups?

Judge the analysis, not the vocabulary. Comments, headers, variable names, and
cited titles are not evidence on their own. But the shape of a file is evidence:
what it loops over, what it subsets, what it repeats.

### Ways this shows up

Some are visible in a single line:

- **Interaction terms.** `y ~ x*z`, `y ~ x + z + x:z`, Stata `c.x#c.z`,
  `i.x#i.z`, `x##z`.
- **Interactions built by hand.** `gen treatXpost = treat*post`, then
  `treat*post` entered as an ordinary regressor. Count these when you can see
  the variable actually being constructed as a product, or when a regressor is
  unmistakably a product of two named variables that also appear in the model.
  Do NOT count a variable merely because its name contains `x` or `_`, and
  never count an empty assignment such as `gen foo_x_bar = .`.
- **Machine learning for conditional effects.** `causal_forest`, `grf`,
  `bartc`, `bcf`, S/T/X/R-learners, double machine learning for CATEs.
- **Varying slopes.** `(1 + x | group)`, `(0 + x | group)`, Stata
  `mixed ... || g: x`, or a `by=` smooth such as
  `gam(y ~ s(time, by = factor(religion)))`.

Others are visible only across the file, and these are the ones most often
missed. Look for them deliberately.

A warning before the list. Every pattern below requires that the thing being
repeated is **a model fit or an effect estimate**. A loop, a subset, or a
`bysort` is not evidence of anything by itself — replication code is full of
loops that read files, append datasets, reshape, rename, or build tables.
Before counting any of these, satisfy yourself that a model is actually being
estimated inside the repetition, and cite the line where that happens. If you
cannot point to the model being fit, the answer is no.

- **The same model fit repeatedly across subgroups.** Two or more calls that
  differ only in which subset they run on — `reg y x if female==1` and
  `reg y x if female==0`, a `foreach`/`forvalues` loop whose body fits a model
  for each group, a `lapply`/`map` that fits a model over subsets, or model
  results written into successive rows of a results matrix such as
  `means[i, ] <- svymean(..., subset(d, g == i))`. A loop that merely
  assembles or reshapes data is not this, however many groups it touches.
- **A sweep across values of a moderator.** `margins, at(z = (0 1 2 3))`,
  `marginaleffects` or `slopes()` evaluated at several covariate values,
  predicted values computed over a grid of some conditioning variable.
- **Quantities compared between groups** as the point of the analysis: group
  means, densities, or effect estimates computed separately and then
  differenced or plotted against each other.
- **A variable constructed in one place and used in another.** If a regressor
  in a model was defined earlier in the file as a product of two variables,
  that is an interaction even though the model line alone does not show it.

### What does not count

- **Fixed effects.** In `fixest`, the pipe separates fixed effects from the
  formula: `feols(y ~ x | country + year)` contains NO interaction. Stata
  `areg`, `xtreg`, `reghdfe`, and `i.country` as a control are likewise not
  interactions.
- **Varying intercepts alone.** `(1 | group)`, `lmer(y ~ x + (1|dept) + (1|year))`,
  Stata `|| cntry:` with nothing after the colon. A varying intercept lets the
  baseline differ, not the effect.
- **A single subset that defines the analysis sample.** One `subset()`,
  `filter()`, or `if` restriction, with no comparison against another subset,
  is sample definition. Two or more contrasting subsets fit with the same model
  is a subgroup analysis.
- **Marginal effects on a model with no interaction.** `margins`, `dydx()`,
  `slopes()`, `comparisons()` compute derivatives. They are evidence only if
  the underlying model is interactive or the call sweeps across values of a
  conditioning variable.
- **Data preparation, plotting, and bookkeeping.** `replace`, `recode`,
  `merge`, `reshape`, `append`, `load`, `read.dta`, `read_csv`, `save`,
  `table()`, `summarize`, `geom_smooth`, bootstrap resampling loops, and
  unique-value extraction are not model fits, whatever they contain.
- **Loops and subsets with no model inside them.** `foreach i in 35 50 80 {
  append using ...}`, `for (state in states) { ... }` that only stacks results
  or writes files, `bysort g: egen m = mean(x)` used to construct a variable,
  and `d[d$treatment == "j1", ]` used to build an analysis frame are all
  ordinary data work. Iterating over groups is how code is written; it is not
  by itself an estimate of how an effect varies.

### Purpose

For files that do estimate heterogeneous effects, judge what the variation is
for: `heterogeneity_claim` if the variation is something the analysis appears
to be interested in, `nuisance_or_fixed_effects` if the interaction only
absorbs variation such as unit-by-time fixed effects, `mixed` if both, `n_a` if
the file does not estimate heterogeneous effects. If you cannot tell, say
`mixed` and lower your confidence. Nothing downstream depends on this field.

### Output

Return a single JSON object and nothing else. No code fence.

```json
{
  "estimates_heterogeneous_effects": true,
  "methods": ["interaction_term", "subgroup_split"],
  "interaction_purpose": "heterogeneity_claim",
  "evidence": ["reg support i.treat##i.female educ age, robust"],
  "confidence": "high",
  "notes": ""
}
```

- `methods` — any of `interaction_term`, `subgroup_split`,
  `ml_heterogeneity`, `hierarchical_varying_slopes`. Empty if the answer is
  false.
- `evidence` — at most three lines copied **verbatim** from the file, the ones
  that best support your judgment. Never invent or paraphrase a line. Empty if
  the answer is false.

  **Your first evidence line must be a line that fits a model or computes an
  estimate.** Not a loop header, not a subset, not a file being read, not a
  variable being generated. If your reason is a pattern across the file, cite
  the model call inside the pattern first, then the line showing the pattern.
  If you cannot find a line where a model is actually estimated, then this file
  does not estimate heterogeneous effects and the answer is false.
- `confidence` — `low`, `medium`, or `high`.
- `notes` — at most one short sentence, only where something is genuinely
  ambiguous. Otherwise empty.

## USER

Language: {language}
File: {filename}
{chunk_note}

```
{code}
```
