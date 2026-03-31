# QST Electron Integrated Audit Simulator v0

This is a **working zero-calibration skeleton** for the QST v7.5 electron branch program.

It implements the program chain:

continuous vortex / radial BVP
-> DSI shell projection
-> E8-channel graph classification
-> QST-Matrix sigma self-consistency loop
-> dominant-shell audit payload extraction
-> effective current emission

## What it does

- solves a **minimal radial stationary envelope** problem for an electron filled-core branch
- emits `D(r)` from `u'(r)/u(r)` and `J_SC(r)` from `kappa * C_core * u(r)^2`
- projects the solution onto DSI shells
- classifies channels from a sparse E8-like selection graph into
  `core-like`, `edge-like`, and `flow-arm/chiral`
- builds a block-tridiagonal QST-Matrix
- emits `sigma` from cross-shell locking:
  `sigma_new = tanh(phi * C_lock)`
- audits candidate branches and extracts a dominant-shell payload
- constructs the minimal electron effective current
  `J_eff,e = (1 + sigma^2 D_e) * \bar Psi_e gamma^mu Psi_e`
  in reduced scalar form

## What it does NOT honestly close

This simulator does **not** claim final first-principles closure of the following open items:

- half-spin admissibility -> branch inevitability
- unique practical E8 basis `chi_a`
- full selected-branch PDE/BVP closure beyond lowest-order working closure
- exact one-loop normalization `N(s_*)`
- exact oscillatory residue amplitude `A`

Those remain explicit placeholders / audit flags in the output.

## Files

- `qst_electron_simulator.py` main simulator
- `run_demo.py` minimal demo run
- `README.md` this file

## Run

```bash
python run_demo.py
```

The demo writes a JSON report in the same folder.
