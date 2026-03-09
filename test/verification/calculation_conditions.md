# CCEX Verification Tests — Calculation Conditions

---

## Common Conditions (All Tests)

| Parameter | Value | Description |
|-----------|-------|-------------|
| **Central spin** | S = 1/2 | Electron spin qubit |
| **α state** | ms = +1/2 | Upper qubit level |
| **β state** | ms = −1/2 | Lower qubit level |
| **Bath spin** | S = 1/2 (e spin) | Electron spin bath |
| **Bath density** | 1 × 10¹⁶ cm⁻³ | Random bath, seed = 2 |
| **Gyromagnetic ratio** | −17608.597 rad/(ms·G) | Electron spin |
| **Magnetic field** | B = 300 G (z-axis) | |
| **r_bath** | 500 Å | Bath inclusion radius |
| **r_dipole** | 300 Å | Cluster pair radius |
| **CCE order** | 2 | Pair cluster |
| **Pulse sequence** | Hahn echo (npulse = 1) | π pulse at t/2 |
| **nstate** | 0 | Ensemble average |
| **addsubclus** | True | Include all sub-clusters |
| **Quantity** | coherence | Re[L(t)] |

---

## Test 1 — ME-CCE Verification

**Goal:** Verify ME-CCE(rate=0) = CCE, and ME-CCE(rate>0) shows additional decoherence.

**Directory:** `1_verify_mecce/`

| Parameter | Value |
|-----------|-------|
| **Time step** (Δt) | 0.04 ms |
| **# of steps** (nstep) | 50 |
| **Time range** | 0 ~ 2.0 ms |

### Runs

| File | Method | Jump operator rate |
|------|--------|--------------------|
| `ccex_cce.json` | CCE | — |
| `ccex_mecce_r0.json` | ME-CCE | 0.0 MHz (rate = 0) |
| `ccex_mecce.json` | ME-CCE | 1.592 × 10⁻⁴ MHz = **1 rad/ms** |

### Jump Operators (ME-CCE with rate)

| Bath | Operator | Rate (MHz) | Rate (rad/ms) |
|------|----------|------------|----------------|
| e (electron) | S⁺ (op: `"-+"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |
| e (electron) | S⁻ (op: `"+-"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |

### Verification Criterion
- ME-CCE(rate=0) vs CCE: **max|diff| < 1×10⁻⁶** → PASS

---

## Test 2 — ME-gCCE Verification

**Goal:** Verify ME-gCCE(rate=0) = gCCE, and ME-gCCE(rate>0) shows additional decoherence.

**Directory:** `2_verify_megcce/`

| Parameter | Value |
|-----------|-------|
| **Time step** (Δt) | 0.1 ms |
| **# of steps** (nstep) | 20 |
| **Time range** | 0 ~ 2.0 ms |

### Runs

| File | Method | Jump operator rate |
|------|--------|--------------------|
| `ccex_gcce.json` | gCCE | — |
| `ccex_megcce_r0.json` | ME-gCCE | 0.0 MHz (rate = 0) |
| `ccex_megcce.json` | ME-gCCE | 1.592 × 10⁻⁴ MHz = **1 rad/ms** |

### Jump Operators (ME-gCCE with rate)

| Bath | Operator | Rate (MHz) | Rate (rad/ms) |
|------|----------|------------|----------------|
| e (electron) | S⁺ (op: `"-+"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |
| e (electron) | S⁻ (op: `"+-"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |

### Verification Criterion
- ME-gCCE(rate=0) vs gCCE: **max|diff| < 1×10⁻⁶** → PASS

---

## Test 3 — pmeCCE Verification

**Goal:** Verify pmeCCE(rate=0) = pCCE, and pmeCCE(rate>0) shows additional decoherence.

**Directory:** `3_verify_pmecce/`

| Parameter | Value |
|-----------|-------|
| **Time step** (Δt) | 0.04 ms |
| **# of steps** (nstep) | 50 |
| **Time range** | 0 ~ 2.0 ms |
| **sK** | 2 | Spins per cluster |
| **max_trial** | 100 | k-means trials |
| **max_iter** | 100 | k-means iterations |

### Runs

| File | Method | Jump operator rate |
|------|--------|--------------------|
| `ccex_cce.json` | CCE | — (reference) |
| `ccex_mecce.json` | ME-CCE | 1.592 × 10⁻⁴ MHz (reference) |
| `ccex_pcce.json` | pCCE | — |
| `ccex_pmecce_r0.json` | pmeCCE | 0.0 MHz (rate = 0) |
| `ccex_pmecce.json` | pmeCCE | 1.592 × 10⁻⁴ MHz = **1 rad/ms** |

### pCCE Clustering Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| **sK** | 2 | Number of bath spins per cluster group |
| **max_trial** | 5000 | Number of k-means restarts |
| **max_iter** | 5000 | Max iterations per k-means run |
| **kmeans_pp** | True (default) | k-means++ initialization |

### Jump Operators (pmeCCE with rate)

| Bath | Operator | Rate (MHz) | Rate (rad/ms) |
|------|----------|------------|----------------|
| e (electron) | S⁺ (op: `"-+"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |
| e (electron) | S⁻ (op: `"+-"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |

### Verification Criterion
- pmeCCE(rate=0) vs pCCE: **max|diff| < 1×10⁻⁶** → PASS

---

## Test 4 — pgCCE Verification

**Goal:** Verify pgCCE runs correctly (t=0 → 1.0), and compare gCCE vs pgCCE behavior.

**Directory:** `4_verify_pgcce/`

| Parameter | Value |
|-----------|-------|
| **Time step** (Δt) | 0.04 ms |
| **# of steps** (nstep) | 50 |
| **Time range** | 0 ~ 2.0 ms |
| **sK** | 2 | Spins per cluster |
| **max_trial** | 100 | k-means trials |
| **max_iter** | 100 | k-means iterations |

### Runs

| File | Method | Description |
|------|--------|-------------|
| `ccex_cce.json` | CCE | Standard CCE (reference) |
| `ccex_gcce.json` | gCCE | Full Hamiltonian CCE (reference) |
| `ccex_pcce.json` | pCCE | pCCE clustering + CCE kernel |
| `ccex_pgcce.json` | pgCCE | pCCE clustering + gCCE kernel |

### pCCE Clustering Parameters (for pCCE and pgCCE)

| Parameter | Value | Description |
|-----------|-------|-------------|
| **sK** | 2 | Number of bath spins per cluster group |
| **max_trial** | 5000 | Number of k-means restarts |
| **max_iter** | 5000 | Max iterations per k-means run |
| **kmeans_pp** | True (default) | k-means++ initialization |

### Verification Criterion
- pgCCE at t = 0: **|L(0) − 1.0| < 1×10⁻⁶** → PASS

---

## Test 5 — pmegCCE Verification

**Goal:** Verify pmegCCE(rate=0) = pgCCE, and pmegCCE(rate>0) shows additional decoherence.

**Directory:** `5_verify_pmegcce/`

| Parameter | Value |
|-----------|-------|
| **Time step** (Δt) | 0.04 ms |
| **# of steps** (nstep) | 50 |
| **Time range** | 0 ~ 2.0 ms |

### Runs

| File | Method | Jump operator rate |
|------|--------|--------------------|
| `ccex_pgcce.json` | pgCCE | — (reference) |
| `ccex_pmegcce_r0.json` | pmegCCE | 0.0 MHz (rate = 0) |
| `ccex_pmegcce.json` | pmegCCE | 1.592 × 10⁻⁴ MHz = **1 rad/ms** |

### pCCE Clustering Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| **sK** | 2 | Number of bath spins per cluster group |
| **max_trial** | 5000 | Number of k-means restarts |
| **max_iter** | 5000 | Max iterations per k-means run |
| **kmeans_pp** | True (default) | k-means++ initialization |

### Jump Operators (pmegCCE with rate)

| Bath | Operator | Rate (MHz) | Rate (rad/ms) |
|------|----------|------------|----------------|
| e (electron) | S⁺ (op: `"-+"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |
| e (electron) | S⁻ (op: `"+-"`) | 1/(2π × 1000) ≈ 1.592 × 10⁻⁴ | 1.0 |

### Verification Criterion
- pmegCCE(rate=0) vs pgCCE: **max|diff| < 1×10⁻⁶** → PASS

---

## Rate Conversion Reference

CCEX uses **MHz** as the unit for jump operator rates.
PyCCE uses **rad/ms**.

$$\text{rate}_\text{CCEX} [\text{MHz}] = \frac{\text{rate}_\text{PyCCE} [\text{rad/ms}]}{2\pi \times 1000}$$

| PyCCE rate (rad/ms) | CCEX rate (MHz) |
|---------------------|-----------------|
| 0.0 | 0.0 |
| 1.0 | 1.592 × 10⁻⁴ |

---

## Method Summary Table

| Method | Clustering | Computation kernel | Dissipation |
|--------|------------|--------------------|-------------|
| CCE | Hash (dipole cutoff) | Product formula, unitary | ✗ |
| gCCE | Hash (dipole cutoff) | Full Hamiltonian, unitary | ✗ |
| pCCE | k-means (sK spins/group) | Product formula, unitary | ✗ |
| pgCCE | k-means (sK spins/group) | Full Hamiltonian, unitary | ✗ |
| ME-CCE | Hash (dipole cutoff) | Product formula, Lindblad | ✓ |
| ME-gCCE | Hash (dipole cutoff) | Full Hamiltonian, Lindblad | ✓ |
| pmeCCE | k-means (sK spins/group) | Product formula, Lindblad | ✓ |
| pmegCCE | k-means (sK spins/group) | Full Hamiltonian, Lindblad | ✓ |

---

## Bath File Info

Bath files are located at:
```
test/mecce_test/comparison_megcce_S05/
├── bath_espin     # xyz positions + spin name ("e") for all bath spins
├── qubit_espin    # qubit position: (0.0, 0.0, 0.0)
└── gyro_espin     # gyromagnetic ratio: e  0.5  -17608.597050
```

Bath generated by PyCCE:
```python
pycce.random_bath('e', [5e3, 5e3, 5e3], density=1e16, density_units='cm-3', seed=2)
# → ~1250 electron spins in a 50×50×50 Å³ box (r_bath=500 Å cuts to ~476 spins)
```
