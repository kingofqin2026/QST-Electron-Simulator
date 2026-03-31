# QST v9 Appendix-Style Audit Statement

## Appendix D.X — Audit-Certified Operating Region for the Electron Reduced-Residue Chain

### D.X.1 Scope and status

This appendix records the current numerical-audit status of the electron reduced-form chain constructed across the V0–V8 simulator sequence.

The chain audited here is

\[
\text{selected branch}
\;\Rightarrow\;
\text{dominant shell}
\;\Rightarrow\;
J^\mu_{\rm eff,e}
\;\Rightarrow\;
R^{(2,e)}_{LR}(s_*)_{\rm reduced}
\;\Rightarrow\;
\mathcal{N}_e^{\rm toy,de}(s_*)
\]

and its audit-compressed output is an **operating envelope**, not a full one-loop theorem.

### D.X.2 Main numerical verdict

The V8 audit extracts a certified operating region with the following summary:

- total sampled points: **12**
- robust-band passes: **8**
- certified points: **8**
- excluded points: **4**

Across the certified subset:

- selected branch remains **filled_core**
- dominant shell remains **n_* = 4**

Hence the currently supported audit statement is:

\[
\boxed{
\text{Within the certified operating region, the reduced electron chain selects}
\; b_* = \text{filled_core}
\;\text{and}
\; n_* = 4.
}
\]

### D.X.3 Certified operating region

The V8 extractor yields the following operating envelope:

\[
\boxed{
\sigma_{\rm init} \in \{0.42,\;0.44,\;0.46\},
\qquad
\text{overlap\_width\_factor} \in \{1.8,\;2.1\},
\qquad
\text{strict\_channel\_floor} \in \{0.05,\;0.07\}
}
\]

with the branch/dominant-shell lock

\[
\boxed{
b_* = \text{filled_core},
\qquad
n_* = 4.
}
\]

### D.X.4 Recommended audit point

The median-anchored certified recommendation from V8 is:

- \(\sigma_{\rm init} = 0.42\)
- overlap width factor \(= 1.8\)
- strict channel floor \(= 0.05\)
- selected branch \(=\) **filled_core**
- dominant shell \(n_* = 4\)
- emitted coherence \(\sigma_{\rm out} \approx 0.742621\)
- raw reduced residue \(\approx -2.210659e-05\)
- clipped reduced residue \(\approx -2.210659e-05\)

This point is not claimed as unique or theorem-level optimal. It is the **recommended certified operating point** under the present audit rule.

### D.X.5 Formal closure classification

The present closure classification is:

- macroscopic foundation: **closed**
- branch-to-field lowest-order working closure: **pass**
- dominant-shell payload emission: **pass**
- reduced residue injection layer: **pass**
- median-anchored robust operating region: **pass**
- full \(R^{(2,e)}_{LR}(s_*)\) one-loop normalized theorem: **open**
- full \(N(s_*)\) derivation: **open**
- oscillatory amplitude \(A\): **open**

Equivalently,

\[
\boxed{
\text{QST v9 currently supports an audit-certified operating region for the reduced electron chain,}
\\
\text{but does not yet support a full first-principles one-loop normalized closure theorem.}
}
\]

### D.X.6 Machine-readable certified subset

The certified subset contains the following parameter points:
- (0.42, 1.8, 0.05) -> clipped residue -2.210659e-05
- (0.42, 1.8, 0.07) -> clipped residue -2.210659e-05
- (0.44, 1.8, 0.05) -> clipped residue -2.348685e-05
- (0.44, 1.8, 0.07) -> clipped residue -2.348685e-05
- (0.44, 2.1, 0.05) -> clipped residue -1.126400e-05
- (0.44, 2.1, 0.07) -> clipped residue -1.126400e-05
- (0.46, 2.1, 0.05) -> clipped residue -1.222687e-05
- (0.46, 2.1, 0.07) -> clipped residue -1.222687e-05

### D.X.7 Conservative interpretation

The strongest defensible claim at this stage is not “full closure,” but rather:

\[
\boxed{
\text{There exists a numerically stable, audit-certified operating region in which the reduced electron chain}
\\
\text{reproducibly selects the filled-core branch and the same dominant shell } n_* = 4.
}
\]

This statement is strong enough to support appendix-level reporting, reduced-form numerical citations, and further side-constraint development, but not yet strong enough to collapse the remaining microscopic bridge into theorem status.