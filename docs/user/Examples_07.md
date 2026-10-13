# 07 h-BN V_B⁻

The boron vacancy in hexagonal BN with natural-abundance ¹⁰B, ¹¹B and ¹⁴N, and
the DFT hyperfine and quadrupole tensors in `bath/flat/` (made with
`example/tools/SpinTensorReader`). `ensemble_CCE/` reads the tensors
(`hamiltonian.common.hyperfine`, `hamiltonian.common.quadrupole`) in three combinations; `hyperfine_mediated_term/`
turns `hamiltonian.cCCE.hfmedi` on and off; `single-sample_gCCE_DPS/` is gCCE2 with and without
the `OH` term. 0.4 to 5 s each (4 ranks, one configuration).

The hyperfine file is `hBN_VB_flat_ATensor_HSE_CellPBE_v3`, in tensor format v3 (a `VASP` first line and the
provenance line `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = VASP OUTCAR, S = 1`), so
CCEX reads it as written; the nearest ¹⁴N has A_zz = 48.6 MHz (measured 47 MHz). Until 2026-10-13 the inputs read
`hBN_VB_flat_ATensor_HSE_CellPBE_tag`, which held twice these numbers and was read 2x; the figures below were made
again with the v3 file. The file's S is checked against the qubit spin, 1
([hyperfine.file](InputFile_General.md#tag_hyperfine_file)).

![07: CCE2 with the DFT hyperfine and quadrupole tensors](../images/tutorials/07_hBN_tensors.png)

![07: the hyperfine-mediated term](../images/tutorials/07_hfmedi.png)

![07: gCCE2 with and without the OH term](../images/tutorials/07_DPS_OH.png)

## The input

`example/tutorial/07_hBN_VB_natab/ensemble_CCE/ccein_hf3_qd2.json`; the other two leave out the hyperfine block (`ccein_hf0_qd2.json`, point dipole) or the quadrupole block (`ccein_hf3_qd0.json`). `run.sh` in the same folder runs it as
`mpirun -n $NP bin/main.out -f <this file> -I <bath> -o results/<name>` and then the analyzer.

```
{
    # Config (General)
    "method"        : "cCCE", 
    "quantity"      : "coherence",
    "Qubit"         : {
        "file"          : "../bath/flat/bath_hBN_VB_flat_Defect",
        # Qubit
        "alphams"      : 1,
        "betams"       : 0
    },
    "gyrofile"      : "../bath/hBN_natab_gyro",
    "bathfile"      : ["../bath/flat/bath_hBN_VB_flat_1"],
    "order"         : 2,
    "bfield"        : 30000,
    "rbath"         : 20,
    "rdip"          : 10,
    "deltat"        : 0.003,
    "nstep"         : 50,
    "nstate"        : 0,

    # Pulse
    "pulse"         : {"count": 1},
    "hamiltonian"   : {
        "common"        : {
            # Bath HF
            "hyperfine"     : {
                "dipole"        : "dft",
                "contact"       : true,
                "file"          : "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3",
                "missing_row"   : "stop"
            },
            # Bath QD
            "quadrupole"  : {
                "structure"     : "plain",
                "file"          : "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"
            } #0 (X), exp, dft(2 widefect) dft(wi+wo defect), dft(wi+wo defect)
        }
    },

    # Output
    "outfile"     : "./results/gCCE2_hBN_VB_natab_1"

    # Bath
    #"Bath" : {
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] },
    #    [ 'bath1', {'xyz': [], 'gyro': '' , 'A': [tensor], 'Q': [tensor] }
    #}

}
```
