# Pulse-related tags
- [npulse](#tag_npulse)       
- [pulsename](#tag_pulsename)    
- [sequence](#tag_sequence)
- [pulseiter](#tag_pulseiter)
- [pulse_time](#tag_pulse_time)
- [detuning_factor](#tag_detuning_factor)

<a id="tag_npulse"></a>
## npulse
<b>Description:</b><br>
The number of pulses. Without [sequence](#tag_sequence) the run builds the default sequence
from it: `0` is a Ramsey experiment (free induction), `1` a Hahn echo, `n` > 1 CPMG with n
equally spaced π pulses. With [pulsename](#tag_pulsename) = `Manual` it is the number of
entries of `sequence` minus one. `method` = `CCE` applies every pulse as an ideal π flip about
X or Y; other angles need gCCE. Required. <br>

<b>Data Type:</b> Integer

<b>Default:</b> - (required)

<b>Example:</b> ``` "npulse" : 1 ```

<b>Related variable:</b> [Pulse](#Pulse)::npulse

<a id="tag_pulsename"></a>
## pulsename
<b>Description:</b><br>
The name of the pulse sequence. Without it the run derives the sequence from
[npulse](#tag_npulse) and reports the name it chose (`Ramsey`, `HahnEcho` or `CPMG`); the
run log prints the sequence as a table of intervals, axes and angles. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `Ramsey`, `HahnEcho`, `CPMG` | the default sequences; which one runs follows `npulse` |
| `WAHUHA` | the four-pulse WAHUHA sequence (90° pulses; gCCE) |
| `Manual` | the sequence written out in [sequence](#tag_sequence) |

<b>Default:</b> - (derived from `npulse`)

<b>Example:</b> ``` "pulsename" : "Manual" ```

<b>Related variable:</b> [Pulse](#Pulse)::pulsename

<a id="tag_sequence"></a>
## sequence
<b>Description:</b><br>
The pulse sequence written out, one entry per pulse: `[fraction, axis, angle]`, where
`fraction` is the free evolution before the pulse as a fraction of the total time (a string
such as `"1/6"` or a number; the fractions must sum to 1), `axis` is `X`, `Y`, `Z` or `I` (no
pulse) and `angle` is the rotation in degrees. Required, with `npulse` + 1 entries, when
[pulsename](#tag_pulsename) = `Manual`; the last entry is the free evolution after the last pulse.
Without it the run builds the default sequence of `npulse` equally spaced π pulses (Hahn echo,
CPMG). `example/tutorial/06_Pulse_Manual/` has Hahn echo, WAHUHA and CPMG written this way. <br>

<b>Data Type:</b> List of [String or Double, String, Double]

<b>Default:</b> - (the default sequence)

<b>Example:</b> ``` "sequence" : [["1/6", "X", 90.0], ["1/6", "Y", 270.0], ["2/6", "Y", 90.0], ["1/6", "X", 270.0], ["1/6", "I", 0.0]] ```

<b>Related variable:</b> [Pulse](#Pulse)::sequence, [Pulse](#Pulse)::pulse_axes, [Pulse](#Pulse)::pulse_angles

<a id="tag_pulseiter"></a>
## pulseiter
<b>Description:</b><br>
`false`: a pulse rotates the qubit array as a whole. `true`: the pulse is applied to each qubit
separately and the rotations are combined. With a finite pulse ([detuning_factor](#tag_detuning_factor)
other than 1) `true` is not implemented. <br>
**Known defect (D47 in (internal notes), not fixed):** with gCCE, `true` also stops with an internal
matrix error for ideal Hahn-echo and CPMG pulses (the per-qubit X operator is left empty) and for a `Manual` sequence (a
matrix of the wrong size); Ramsey (no pulse) is not affected. Use the default `false` until it is fixed. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "pulseiter" : true ```

<b>Related variable:</b> [Pulse](#Pulse)::pulseiter

<a id="tag_pulse_time"></a>
## pulse_time
<b>Description:</b><br>
The duration of a pulse. Used only by the finite-pulse propagation, i.e. with
[detuning_factor](#tag_detuning_factor) other than 1. <br>

<b>Data Type:</b> Double

<b>Unit:</b> ns

<b>Default:</b> `0.0`

<b>Example:</b> ``` "pulse_time" : 10.0 ```

<b>Related variable:</b> [Pulse](#Pulse)::pulse_time

<a id="tag_detuning_factor"></a>
## detuning_factor
<b>Description:</b><br>
gCCE. `1.0` applies every pulse as an ideal instantaneous rotation. Any other value propagates
each pulse as a finite drive of duration [pulse_time](#tag_pulse_time): the drive frequency is
`detuning_factor` times the qubit transition frequency, the detuning is their difference, and the
Rabi frequency follows from the angle and the duration. <br>

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "detuning_factor" : 0.98 ```

<b>Related variable:</b> [Pulse](#Pulse)::detuning_factor
