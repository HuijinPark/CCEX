# Pulse-related tags
The `pulse` block (new on 2026-10-13: the keys `npulse`, `pulsename`, `sequence`, `pulseiter`, `pulse_time` and
`detuning_factor` were top-level until then, see the table of [InputFile.md](InputFile.md#OldNames)). The block and its
`count` are required.

``` "pulse" : {"count": 1} ```

- [pulse](#tag_pulse)
    - [count](#tag_pulse_count)
    - [name](#tag_pulse_name)
    - [sequence](#tag_sequence)
    - [per_qubit](#tag_pulse_per_qubit)
    - [duration](#tag_pulse_duration)
    - [detuning_factor](#tag_detuning_factor)

<a id="tag_pulse"></a>
## pulse
<b>Description:</b><br>
The pulse sequence of the run. Without it the run stops (`pulse is not found in the input file`). <br>
<b>Data Type:</b> Object

<b>Example:</b> ``` "pulse" : {"count": 2, "name": "Manual", "sequence": [["1/4", "I", 0], ["1/2", "X", 180], ["1/4", "X", 180]]} ```

<a id="tag_pulse_count"></a>
## pulse.count
<b>Description:</b><br>
The number of pulses (`npulse` until 2026-10-13). Without [sequence](#tag_sequence) the run builds the default sequence
from it: `0` is a Ramsey experiment (free induction), `1` a Hahn echo, `n` > 1 CPMG with n
equally spaced π pulses. With [name](#tag_pulse_name) = `Manual` it is the number of
entries of `sequence` minus one. `method` = `cCCE` applies every pulse as an ideal π flip about
X or Y; other angles need gCCE. Required. <br>

<b>Data Type:</b> Integer

<b>Default:</b> - (required)

<b>Example:</b> ``` "pulse" : {"count": 1} ```

<b>Related variable:</b> [Pulse](#Pulse)::npulse

<a id="tag_pulse_name"></a>
## pulse.name
<b>Description:</b><br>
The name of the pulse sequence (`pulsename` until 2026-10-13). Without it the run derives the sequence from
[count](#tag_pulse_count) and reports the name it chose (`Ramsey`, `HahnEcho` or `CPMG`); the
run log prints the sequence as a table of intervals, axes and angles. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `Ramsey`, `HahnEcho`, `CPMG` | the default sequences; which one runs follows `count` |
| `WAHUHA` | the four-pulse WAHUHA sequence (90° pulses; gCCE) |
| `Manual` | the sequence written out in [sequence](#tag_sequence) |

<b>Default:</b> - (derived from `count`)

<b>Example:</b> ``` "pulse" : {"count": 2, "name": "Manual", "sequence": [...]} ```

<b>Related variable:</b> [Pulse](#Pulse)::pulsename

<a id="tag_sequence"></a>
## pulse.sequence
<b>Description:</b><br>
The pulse sequence written out, one entry per pulse: `[fraction, axis, angle]`, where
`fraction` is the free evolution before the pulse as a fraction of the total time (a string
such as `"1/6"` or a number; the fractions must sum to 1), `axis` is `X`, `Y`, `Z` or `I` (no
pulse) and `angle` is the rotation in degrees. Required, with `count` + 1 entries, when
[name](#tag_pulse_name) = `Manual`; the last entry is the free evolution after the last pulse.
Without it the run builds the default sequence of `count` equally spaced π pulses (Hahn echo,
CPMG). `example/tutorial/06_Pulse_Manual/` has Hahn echo, WAHUHA and CPMG written this way. <br>

<b>Data Type:</b> List of [String or Double, String, Double]

<b>Default:</b> - (the default sequence)

<b>Example:</b> ``` "pulse" : {"count": 4, "name": "Manual", "sequence" : [["1/6", "X", 90.0], ["1/6", "Y", 270.0], ["2/6", "Y", 90.0], ["1/6", "X", 270.0], ["1/6", "I", 0.0]]} ```

<b>Related variable:</b> [Pulse](#Pulse)::sequence, [Pulse](#Pulse)::pulse_axes, [Pulse](#Pulse)::pulse_angles

<a id="tag_pulse_per_qubit"></a>
## pulse.per_qubit
<b>Description:</b><br>
(`pulseiter` until 2026-10-13.) `false`: a pulse rotates the qubit array as a whole. `true`: the pulse is applied to each qubit
separately and the rotations are combined. With a finite pulse ([detuning_factor](#tag_detuning_factor)
other than 1) `true` is not implemented. <br>
**Known defect (D47 in (internal notes), not fixed):** with gCCE, `true` also stops with an internal
matrix error for ideal Hahn-echo and CPMG pulses (the per-qubit X operator is left empty) and for a `Manual` sequence (a
matrix of the wrong size); Ramsey (no pulse) is not affected. Use the default `false` until it is fixed. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "pulse" : {"count": 1, "per_qubit": true} ```

<b>Related variable:</b> [Pulse](#Pulse)::pulseiter

<a id="tag_pulse_duration"></a>
## pulse.duration
<b>Description:</b><br>
The duration of a pulse (`pulse_time` until 2026-10-13). Used only by the finite-pulse propagation, i.e. with
[detuning_factor](#tag_detuning_factor) other than 1. <br>

<b>Data Type:</b> Double

<b>Unit:</b> ns

<b>Default:</b> `0.0`

<b>Example:</b> ``` "pulse" : {"count": 1, "duration": 10.0, "detuning_factor": 0.98} ```

<b>Related variable:</b> [Pulse](#Pulse)::pulse_time

<a id="tag_detuning_factor"></a>
## pulse.detuning_factor
<b>Description:</b><br>
gCCE. `1.0` applies every pulse as an ideal instantaneous rotation. Any other value propagates
each pulse as a finite drive of [duration](#tag_pulse_duration): the drive frequency is
`detuning_factor` times the qubit transition frequency, the detuning is their difference, and the
Rabi frequency follows from the angle and the duration. <br>

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "pulse" : {"count": 1, "duration": 10.0, "detuning_factor": 0.98} ```

<b>Related variable:</b> [Pulse](#Pulse)::detuning_factor
