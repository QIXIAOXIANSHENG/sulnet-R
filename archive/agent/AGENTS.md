# AGENTS.md — Roxygen documentation for this R package

## Assignment

Add or improve **roxygen2 source comments for the main R functions** in this existing package. This is a documentation-only task. The package and its algorithms are already implemented; do not redesign them. Use the handwritten mathematical specification below to understand the intended meaning of names, but treat the **actual R/Fortran implementation as authoritative** for function signatures, default values, supported inputs, array layouts, edge-case behavior, output objects, and export status.

**Do not run `devtools::document()`, `roxygen2::roxygenise()`, or any equivalent documentation-generation command.** Do not generate or edit `man/*.Rd` or `NAMESPACE`. The user will generate the documentation locally later.

## Mathematical and statistical context

### 1. Original model and coefficient reparameterization

The underlying linear regression model is

\[
y = X\gamma + \varepsilon.
\]

For predictor \(j\), use a supplied or data-derived weight \(\beta_j\) and solve for \(\gamma_j\) through

\[
\gamma_j = \beta_j\theta_j.
\]

Here `beta` is the initial coefficient/weight vector, `theta` is the coefficient vector fitted by the penalized algorithm, and `gamma` is the final coefficient vector on the original predictors. **Do not confuse these three quantities in documentation.**

### 2. Construction of the transformed design matrix `F`

The handwritten procedure is:

1. If no user prior is available for predictor \(j\), obtain \(\beta_j=b_j\), the univariate data-derived coefficient/weight. The note writes \(b_j=\operatorname{cov}(X_j,y)\). Whether the implementation standardizes predictors or computes an ordinary univariate regression slope must be determined from the code: **do not silently replace covariance by covariance/variance or vice versa.**
2. If a user prior \(P_j\) is available, the schematic scaling is
   \[
   \beta_j = \frac{\operatorname{median}_l |r_l|}{\operatorname{median}_l |P_l|}\,P_j,
   \qquad r_l=\operatorname{cor}(X_l,y).
   \]
   Confirm the exact subset of indices, missing-prior behavior, and safeguards against zero or non-finite values in the actual code. A user may supply priors only for some predictors.
3. Compute the per-predictor offset
   \[
   \beta_{0j}=\operatorname{mean}(y-X_j\beta_j),
   \]
   and the transformed predictors
   \[
   F_{ij}=X_{ij}\beta_j+\beta_{0j}.
   \]

The algorithm then fits `theta` using `F`. Explain any centering, scaling, weights, intercept treatment, or sparse-matrix support **only as actually implemented**. `F` is a transformed design matrix, not an unrelated input.

### 3. Penalized objective and parameterizations

The schematic *zeta–tau* objective is

\[
\frac{1}{2n}\lVert y-F\theta\rVert_2^2
+\zeta\sum_{j=1}^p pf1_j|\theta_j|
+\tau\sum_{j=1}^p pf3_j(\theta_j)_-
+\lambda_2\sum_{j=1}^p pf2_j\theta_j^2,
\]

where \((t)_-=\max(-t,0)\) is the negative part. In particular:

- `zeta` multiplies an absolute-value (L1) penalty; it is **not** merely a penalty on positive coefficients.
- `tau` adds an asymmetric penalty for negative coefficient values.
- `lambda2` controls the separately written quadratic penalty, if implemented in the corresponding routine. Do not insert an unobserved factor of `1/2` or conflate `lambda2` with `lambda`.
- `pf1`, `pf3`, and `pf2` are different penalty-factor vectors; keep their indices and meaning distinct. Confirm exact defaults and allowed values from the code.

The *lambda–alpha* version uses

\[
\zeta=\lambda(1-\alpha),\qquad \tau=\lambda\alpha,
\]

so that the L1 and negative-part terms become

\[
\lambda(1-\alpha)\sum_j pf1_j|\theta_j|
+\lambda\alpha\sum_j pf3_j(\theta_j)_-.
\]

The algorithm is intended for `alpha >= 0.5`. The handwritten sum–ratio parameterization is

\[
S=\lambda=\zeta+\tau,\qquad
r=\frac{\alpha}{1-\alpha}=\frac{\tau}{\zeta}\quad (\alpha<1).
\]

Here `r` is a penalty *rate*, distinct from the correlations `r_l` used in prior scaling. `alpha = 1` is a boundary case for the ratio parameterization; document how the package handles it only after inspecting the actual code. Do not assume all three parameterizations are exposed by every function or that they use the same grid traversal.

Do not describe this objective as the ordinary `glmnet` elastic-net objective: the extra negative-part term changes its meaning.

### 4. Returning coefficients to the original scale

The handwritten reconstruction is

\[
\widehat\gamma_j=\widehat\theta_j\beta_j,
\qquad
\widehat\gamma_0=\widehat\theta_0+\sum_{j=1}^p\widehat\theta_j\beta_{0j}.
\]

Use this distinction when documenting fitted values, `coef()`, `predict()`, and returned `gamma`/`gamma0` fields, **if the implementation follows this reconstruction**. Confirm how intercepts are represented and which coefficient scale each output exposes.

## Repository exploration before editing

1. Read `DESCRIPTION`, `NAMESPACE`, and the files in `R/`. Look at relevant `src/` routines and `.Fortran`/`.Call` wrappers only when needed to establish argument meaning, indexing, failure codes, or returned object structure.
2. Identify the actual exported/high-level fitting functions and their associated prediction, coefficient, cross-validation, print, and plot methods **if present**. Also inspect important path/grid and output-conversion helpers when their documentation is useful. Do not assume any named function or API exists.
3. This implementation borrows organizational patterns and variable names from `glmnet` and `gcdnet`. Use that familiarity to interpret code, **not** to copy documentation, invent options, or impose their default settings or objective functions.
4. Preserve existing roxygen comments wherever correct, improving only missing, misleading, redundant, or incomplete parts. Check the existing project's roxygen Markdown convention before using Markdown-specific syntax.

Potential identifiers to search for, **not a mandated function list**, include `zeta`, `tau`, `lambda`, `alpha`, `lambda2`, `sum`, `rate`, `beta`, `beta0`, `theta`, `gamma`, `gamma0`, `f`, `pf1`, `pf2`, `pf3`, `getoutputLA`, `gridIndexLA`, `gridIndexSR`, and `suniZTNET`. Some of these may be Fortran symbols or internal helpers rather than exported R functions; annotate R functions, not Fortran subroutines.

## Roxygen2 authoring rules

- Add comments immediately before each selected R function definition, using `#'` and valid roxygen2 tags.
- For public fitting functions, write a concise, informative title, a one- or two-sentence description, relevant `@details` for the reparameterization/objective, **one accurate `@param` per actual formal argument**, and a concrete `@return` describing the implemented return type, dimensions, and significant named components.
- For applicable methods, explain whether outputs refer to `theta` or the reconstructed `gamma`, the intercept convention, the effect of `type` and `s`/penalty selection, and the shape of returned matrices or arrays. Only document arguments and options that are actually accepted.
- For path/grid functions, document the parameterization (`zeta`–`tau`, `lambda`–`alpha`, or `sum`–`rate`), grid-order semantics, indexing conventions, warm-start meaning, error handling, and output flattening **when verified in source**. For helpers such as `getoutputLA`, inspect how `nlambda`, `nalpha`, `napenalty`, and corresponding paths are handled before documenting them.
- Preserve the existing API, including parameter spellings, defaults, vector recycling/validation, `...`, class names, and return component names. Never substitute more familiar `glmnet` terminology when it disagrees with this package.
- Retain existing `@export` / S3 registration conventions. Add `@export` only if it is already clear that the function or method is intended to be exported; never change export policy as part of documentation. Use `@noRd` for significant internal helpers only when documenting them is useful; do not automatically document every helper.
- Use `@inheritParams` only if an appropriate, existing roxygen source block makes the inherited parameter descriptions correct. Use `@seealso` only for verified functions. Add `@importFrom` only for genuinely used imports and only when needed by established package conventions.
- Include `@examples` only when a small, valid, self-contained example is readily supported by the actual API. Avoid costly examples, unsupported inputs, or dependency on external files/network services. An accurate skeleton without examples is better than a fabricated example.
- Prefer compact, specific language over filler. Explain what a parameter **does** and its interpretation, not merely its data type. Use placeholder descriptions such as `TODO`, `parameter x`, or `a numeric value` when theoretical details are missing or uncertain.
- Use valid roxygen math markup (`\eqn{}` or `\deqn{}` inside the R source where appropriate). Write the mathematics accurately; do not paste raw Markdown LaTeX delimiters into roxygen tags.
- If an implementation detail is ambiguous or inconsistent with this handwritten design, **do not guess or silently change code**. Document what is verified and flag unresolved questions in the final summary.

## Editing constraints

- Only edit `.R` source files in `R/` to add or revise roxygen comments, except when a different location is absolutely necessary and justified. Do not change executable R code, C/Fortran code, tests, installed data, dependencies, or package configuration.
- Do not rename functions or variables, reorder arguments, add checks, fix algorithms, or change numerical behavior—even if you notice a possible bug. Report suspected issues separately.
- Do not run `devtools::document()`, `roxygen2::roxygenise()`, `R CMD roxygen`, or any other step that generates Rd files or rewrites `NAMESPACE`.
- You may use read-only repository searches and `git diff --check` to review the documentation edits. Do not run package checks, tests, installations, or numerical benchmarks unless explicitly requested.
- Respect pre-existing uncommitted changes; do not overwrite unrelated edits.

## Expected result and report

Make the roxygen edits directly in the existing R files. At the end, briefly report:

1. Which functions/files received documentation.
2. Any important mathematical or API ambiguities discovered from the source (particularly covariance versus univariate slope, optional priors, coefficient scales, `alpha = 1`, intercept handling, or penalty factors).
3. Confirmation that **no code behavior was modified and `document()` was not run**.

Do not generate the `.Rd` files; leave that final step to the user.
