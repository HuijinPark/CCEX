# DEER bath-pulse regression test

Six short CCEX runs that check the bath-pulse path against answers that are known rather than
against a previous run. Twelve spins, cluster order 1, non-sampled — a few seconds in total.

```sh
sh run_test.sh                     # uses the binary at the repo root
sh run_test.sh /path/to/main.out
NP=4 sh run_test.sh                # more ranks
```

Exit status is 0 only if every check passes, so it can be wired into CI.

## Why order 1

At cluster order 1 the bath spins do not interact, so the echo is a plain product over them and
can be written down:

```
V(t) = prod_j cos( a_j · t · min(f, 1-f) ),      a_j = D (1 - 3cos²θ_j) / r_j³
```

with `D/2π = 52.04 MHz·nm³` and `f` the fraction of the total time at which the bath π lands.
At `f = 0.5` the pulse sits on the echo centre and the argument becomes the familiar `a_j t/2`.

So the test does not need a stored reference curve — it compares against the analytic answer.

## What is checked

| check | why it exists |
|---|---|
| no pump → `V ≡ 1` exactly | order 1 refocuses every static coupling; any drift here means the qubit ZFS or the sequence is wrong |
| off-resonant pump ≡ no pump, bit for bit | the frequency selection rule must not touch spins outside the matching window |
| multi-tone `[resonant, far]` ≡ single resonant | array form of `bpulse_energy_shift` must behave as "any tone in range drives the spin" |
| multi-tone `[far, far]` ≡ no pump | and must not drive anything when no tone matches |
| pump at the echo centre ≡ closed form | the actual physics, with no free parameter |
| **pump next to the qubit π applied once, not twice** | regression for the timing-tolerance bug, see below |
| `wiDiv ≡ noDiv` | at order 1 every inclusion–exclusion coefficient is +1, so no division exists; if these ever differ the cluster bookkeeping has changed |

## The timing-tolerance regression

`build_tot_sequence()` splits the sequence at pulse times and then asks, at each segment edge,
whether a pulse belongs there. The qubit branch matched to `1e-12`, but the bath branch used
`0.01` — a hundredth of the whole sequence, 75 ns when 2τ = 7.5 µs.

A bath π landing that close to the qubit π therefore matched at **both** edges and was applied
twice. Two π pulses cancel, so the pump silently did nothing and the result collapsed onto a
plain Hahn echo — no error, no warning, just a wrong number.

Measured on a real bath before the fix:

```
NV π offset      0 ns     15 ns     38 ns     50 ns     75 ns
before         0.4931    0.9999    0.9999    0.9999    0.5081
after          0.4931    0.4961    0.5006    0.5031    0.5081
```

The `nearpulse` case puts the bath π at `f = 0.505` — inside the old blind spot — and requires
it to match the closed form and to differ from the no-pump curve.

## Files

```
bath_DEERtest_12spin.txt   12 P1 spins on diamond lattice sites, 65-120 A, fixed
avaax_DEERtest             axis assignment; 5 spins on the pumped line
{case}/ccein.js            one input per case (GYROFILE is substituted at run time)
run_test.sh                runs the six cases and calls the checker
check_deer.py              the checks; exits non-zero on failure
```

Couplings of the driven spins are 0.01–0.30 MHz, so the echo decays over a few microseconds and
the 0–10 µs window in the inputs covers it.
