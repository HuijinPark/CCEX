/**
 * @file example_input.cpp
 * @brief The example input (ccex_example.jsonc) that main.out -e writes.
 */
// The input template `main.out -e` writes, the one place it is kept. The active part is the calculation of
// example/tutorial/00_quickstart (the NV-13C cCCE2 Hahn echo, ensemble, rbath 50 A, 2 ms) without its files: gyrofile is
// commented (the built-in isotope table) and bathfile is a name the user replaces, so the template does not run as it
// is. Every other input key the readers accept (src/input/option.cpp, the adaptive options table, the automatic-
// convergence runner) is a commented line with its allowed values and default. Comments are // (vim and VS Code show
// a .jsonc file's // as comments), every trailing comment starts at column 41, and the ':' of one object's keys share a
// column ("key" : value). test/code_verification checks it: with a real bath it runs, no accepted key is missing from
// it, the comment column and the ':' column.
#include "../../include/input/example_input.h"

const char* const kExampleInputName = "ccex_example.jsonc";

const char* const kExampleInput = R"CCEX(// ccex_example.jsonc - an input template written by `main.out -e` (CCEX)
//
// This is a template: it does not run until "bathfile" below names a real bath file. A ready bath for the qubit
// written here is example/tutorial/00_quickstart/bath/bath_Diamond_NV_natab_1. With a bath file in place, run
//
//     mpirun -n 4 main.out -f ccex_example.jsonc -o results
//
// The uncommented part is the calculation of example/tutorial/00_quickstart: the Hahn echo of an NV centre (spin 1)
// in diamond with 13C, conventional CCE (cCCE) to order 2, the bath in its unpolarised (ensemble) state, 13C within
// 50 Angstrom, 2 ms. It writes ccex_example_noDiv and ccex_example_wiDiv (output.file defaults to the input file name
// without .jsonc or .json, with -o DIR or output.path they go to that folder). Every other key the input accepts is a
// commented line with its allowed values and its default: remove the leading // to use one. A //, # or ! outside a
// string starts a comment to the end of the line. Keys are read without regard to case. A key no part of the run
// reads is reported, and a likely typo of a key stops the run. The keys are described in docs/user/InputFile.md and
// the pages it links.
//
// Units: a number is in the unit its key names below. The keys marked (unit) also take {"value": ..., "unit": "..."}
// (the value a number or a list, as the plain form): time s, ms, us, ns. Field G, mT, T. Length A (Angstrom), nm.
// Frequency Hz, kHz, MHz, GHz. "deltat": {"value": 5, "unit": "us"} is "deltat": 0.005.
{
    // ---- method and quantity ----------------------------------------------------------------------------------------
    "method"           : "cCCE",        // cCCE (default), gCCE, meCCE, megCCE (noise is a reserved name, refused)
    "quantity"         : "coherence",   // coherence (default) or dm (the qubit density matrix, gCCE)

    // ---- the qubit --------------------------------------------------------------------------------------------------
    // Form 2 (used here): the qubits written out. Form 1 (at the end of the block, in place of "qubit") reads one
    // qubit's position from a file.
    "Qubit"            : {
        "qubit"       : [
            {
                "name"     : "NV",      // default "q0", "q1", ... by position in the list
                "spin"     : 1.0,       // default 1.0
//              "gyro"     : -17608.597050,
                                        // rad kHz / G, default -17608.597050 (the electron)
                "xyz"      : [73.73895, 75.72808, 73.97033],
                                        // position, Angstrom (unit), in the frame of the bath file (required)
//              "detuning" : 0.0,       // kHz (unit), default 0.0
                "alphams"  : 1,         // ms of the alpha state, default 1.0
                "betams"   : 0          // ms of the beta state, default 0.0
            }
        ]
//      "intmap"      : [               // qubit-qubit tensors, default zero self terms, point dipole between qubits
//          {
//              "between"      : ["NV", "NV"],
                                        // two qubit names, the same name twice is a self term (zfs)
//              "tensor"       : [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]],
                                        // 3x3 (or 9 numbers)
//              "unit"         : "kHz", // Hz, kHz (default), MHz, GHz
//              "tensor_frame" : "bath" // bath or qubit: the basis of the components (with the rotation on)
//          }
//      ],
//      "alpha_state" : {"index": 0},   // {"index": n} (n-th eigenstate by descending energy) or
//      "beta_state"  : {"vector": [0, 1, 0]},
                                        // {"vector": [...]}, default: the product of the qubits' alphams/betams
//      "psi0"        : [0.7071, 0.7071, 0.0],
                                        // initial qubit state, default (|alpha> + |beta>) normalized
//      "file"        : "defect_position",
                                        // form 1: a file whose first line is the qubit position x y z (Angstrom)
//      "spin"        : 1.0,            // form 1: qubit spin, default 1.0
//      "zfs"         : [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                                        // form 1: 3x3 zero-field splitting in kHz, default zero tensor
//      "alphams"     : 1,              // form 1: ms of the alpha state, default 1.0
//      "betams"      : 0               // form 1: ms of the beta state, default 0.0
    },

    // ---- the bath configuration -------------------------------------------------------------------------------------
    "configuration"    : {
//      "path"              : "data",   // folder of the file names of this block (relative to the working directory, or
                                        // absolute, used as it is), default "."
//      "k_range"           : [1, 10],  // [first, last]: every k is a configuration of its own, outputs
                                        // <output.file>_k<k>_..., {k} in bathfile, bath_state_file, defect_axis_file,
                                        // defect_state_file and pcce.cache.file is k ({{k}} a literal {k}). Without
                                        // k_range {k} is literal
//      "gyrofile"          : "gyro.txt",
                                        // isotope, spin, gyromagnetic ratio per line, default: the built-in isotope
                                        // table
        "bathfile"          : ["my_bath.txt"]
                                        // FILL IN: one or more bath files (merged), a count line, then x y z isotope
                                        // per line, in the frame of the qubit xyz above. A ready one for this NV:
                                        // example/tutorial/00_quickstart/bath/bath_Diamond_NV_natab_1
//      "bath_shift"        : [[0.0, 0.0, 0.0]],
                                        // Angstrom (unit), one [x, y, z] per bath file, default [0, 0, 0] each
//      "bath_state_file"   : ["state_"],
                                        // one prefix per bath file of the state files <prefix>1, <prefix>2, ...
                                        // (nstate > 0), default: absent or a missing file draws them from seed
//      "defect_axis_file"  : ["axis"], // one file per bath file: the axis index of each Defect spin, default: drawn
                                        // from seed
//      "defect_state_file" : ["stateEx_"],
                                        // one prefix per bath file of the Defect onsite-spin state files, default:
                                        // drawn from seed
//      "nstate"            : 0,        // 0 (default): ensemble, > 0: that many sampled bath states
//      "seed"              : 1         // random seed shared by every rank, default: generated on rank 0 and reported
    },
    "order"            : 2,             // largest cluster size (required)
    "bfield"           : 500,           // Gauss (unit): Bz, or the list Bx, By, Bz (required, 0 for no field)
    "rbath"            : 50,            // Angstrom (unit): the bath radius, or the list inner, outer (inner default
                                        // 0) (required)
    "rdip"             : 10,            // Angstrom (unit): the pair distance, or the list inner, outer (inner default
                                        // 0) (required, except with cluster.mode bundle)
    "deltat"           : 0.001,         // ms (unit): the time step (required)
    "nstep"            : 2000,          // number of time steps (required)

    // ---- the pulse sequence -----------------------------------------------------------------------------------------
    "pulse"            : {
        "count"           : 1           // number of pi pulses: 0 Ramsey, 1 HahnEcho, > 1 CPMG (required)
//      "name"            : "Manual",   // Ramsey, HahnEcho, CPMG, WAHUHA, Manual, default None (derived from count)
//      "sequence"        : [["1/2", "X", 180.0], ["1/2", "I", 0.0]],
                                        // Manual: count + 1 rows [fraction, axis X|Y|Z|I, angle in degrees], the
                                        // fractions sum to 1, default: the default sequence
//      "per_qubit"       : false,      // pulses applied qubit by qubit, default false
//      "duration"        : 0.0,        // pulse duration, ns (unit), default 0.0 (ideal)
//      "detuning_factor" : 1.0,        // default 1.0
//      "bath"            : {           // DEER bath pulses (gCCE, megCCE, docs/methods/deer.md)
//          "count"        : 1,         // default 0
//          "sequence"     : [["1/2", "X", 180.0], ["1/2", "I", 0.0]],
                                        // count + 1 rows, as pulse.sequence
//          "defect"       : "P1",      // the driven Defect species name (required when count > 0)
//          "energy_shift" : [0.0],     // MHz (unit), a number or a list of tones (required when count > 0)
//          "bandwidth"    : 0.0015915494309189536,
                                        // MHz (unit) > 0: the half width of the frequency window of the bath pulse (a
                                        // spin is driven when its detuning is this close to a tone), default 10/(2 pi
                                        // 1000) MHz = 10 rad kHz, the fixed value before 2026-10-16. The physical
                                        // value is to be decided with the DEER author
//          "spin"         : 0.5,       // default 0.5
//          "alphams"      : 0.5,       // default 0.5
//          "betams"       : -0.5       // default -0.5
//      }
    }                                   // the last active key: a key uncommented below needs a ',' after this '}'

    // ---- the Hamiltonian --------------------------------------------------------------------------------------------
//  "hamiltonian"      : {
//      "common" : {
//          "qubit_bath_mean_field" : false,
                                        // the Overhauser field of the bath on the qubits, default false
//          "bath_bath_mean_field"  : "z",
                                        // off, z (default) or xyz: static field of the other bath spins (nstate > 0)
//          "hyperfine"             : { // default {"dipole": "point", "contact": false}
//              "dipole"      : "point",
                                        // point (default) or dft (from file)
//              "contact"     : false,  // Fermi contact term from file, default false
//              "file"        : "hyperfine_file",
                                        // required with dipole dft or contact true (CCEX-tensor format)
//              "missing_row" : "stop"  // stop (default) or point_dipole for a spin the file does not hold
//          },
//          "quadrupole"            : { // "off" (default) or this object
//              "structure"           : "plain",
                                        // plain, strained or bubble
//              "file"                : "efg_file",
                                        // required with every structure
//              "file_without_defect" : "efg_file_pristine",
                                        // required with strained and bubble
//              "cell"                : [10.0, 10.0, 10.0]
                                        // required with strained and bubble
//          }
//      },
//      "cCCE"   : {                    // cCCE terms (meCCE reads them too)
//          "knight"  : false,          // default false
//          "hfmedi"  : false,          // hyperfine-mediated term, default false
//          "delta_A" : true            // default true, false replaces a cluster's hyperfine tensors by their mean
//      },
//      "gCCE"   : {},                  // no settings yet (also "noise": {}, a reserved method)
//      "meCCE"  : {                    // or "megCCE": master-equation dissipators (docs/methods/master-equation.md)
//          "jump_operators" : [
//              {"bath_name": "13C", "operator": "z", "rate": 0.001}
                                        // operator - + z x y, rate MHz (unit) >= 0
//          ]
//      }
//  },

    // ---- the coordinate frame ---------------------------------------------------------------------------------------
//  "frame"            : {
//      "rotation" : {                  // the qubit-aligned frame: on when the block is given, absent = no rotation
                                        // (the new z axis is qubit_axis, bath_axis only fixes the azimuth of the new x
                                        // axis)
//          "enabled"              : true,
                                        // default true, false turns a written block off
//          "bath_axis"            : [0, 0, 1],
                                        // in the bath file's frame, default [0, 0, 1]
//          "qubit_axis"           : [1, 1, 1],
                                        // the reference qubit's symmetry axis in the bath file's frame, default
                                        // [0, 0, 1]
//          "hf_tensor_frame"      : "bath",
                                        // bath or qubit: basis of the hyperfine file (needed when it is read)
//          "qd_tensor_frame"      : "bath",
                                        // bath or qubit: basis of the quadrupole file (needed when it is read)
//          "qubit_position_frame" : "bath",
                                        // bath (the only value), frame of every Qubit xyz
//          "reference_qubit"      : "NV"
                                        // default: the first qubit (required with more than one qubit)
//      }
//  },

    // ---- paramagnetic defects in the bath ---------------------------------------------------------------------------
//  "Defect"           : {
//      "species" : [
//          {
//              "name"               : "P1",
                                        // the species name in the bath file
//              "n_onsite_spins"     : 1,
                                        // spins on the defect site (nuclear spins), default 0
//              "n_axis_states"      : 4,
                                        // axis states of the defect (its possible orientations), default 0
//              "coordinate_frame"   : "bath",
                                        // bath or qubit, default qubit for the legacy forms
//              "axis"               : [1.0, 1.0, 1.0],
                                        // symmetry axis (needs coordinate_frame), default absent
//              "types"              : ["14N"],
                                        // isotope of each nuclear spin
//              "spins"              : [1.0],
//              "gyros"              : {"values": [3.0777], "unit": "MHz/T"},
                                        // a list (radkHz/G) or {values, unit}: radkHz/G, rad/ms/G, kHz/G, MHz/G, MHz/T,
                                        // Hz/T, radMHz/T, rad/us/T, rad/s/T
//              "quadrupole_moments" : {"values": [20.44], "unit": "mbarn"},
                                        // a list (fm2) or {values, unit}: fm^2 fm2 m^2 m2 cm^2 cm2 barn b mbarn
                                        // millibarn mb, default 0
//              "nuclear_positions"  : [[1, "14N", [0.0, 0.0, 1.0]]],
                                        // [configuration, type, xyz], default 0
//              "hyperfine"          : {"values": [[1, "14N", [81.3, 0, 0, 0, 81.3, 0, 0, 0, 114.0]]], "unit": "MHz"},
                                        // a list (MHz) or {values, unit}: Hz kHz MHz GHz
//              "efg"                : {"values": [[1, "14N", [0, 0, 0, 0, 0, 0, 0, 0, 0]]], "unit": "Hartree/Bohr^2"},
                                        // a list (au) or {values, unit}: Hartree/Bohr^2 Ha/Bohr^2 Ha/a0^2 a.u. au
                                        // atomic_unit V/m^2 V/A^2 ...
//              "zfs"                : {"D": 2870.0, "E": 0.0, "unit": "MHz"},
                                        // {D, E, unit}, {tensor, unit} or the legacy list [[configuration, "e", 9
                                        // numbers]], default zero tensor
//              "detuning"           : [[1, "e", 84.9]]
                                        // MHz per configuration, or {"values": [...], "unit": "MHz"}, default 0
//          }
//      ]
//  },

    // ---- the clusters -----------------------------------------------------------------------------------------------
//  "cluster"          : {
//      "mode"      : "standard",       // standard (default), pcce (cCCE, gCCE, meCCE, megCCE), bundle (cCCE, gCCE),
                                        // dsj, itb, dsjitb (cCCE), the pcce and bundle blocks hold their settings
//      "selection" : {
//          "add_subclusters" : true,   // true (default) or false
//          "max_per_order"   : [0, 0, 0]
                                        // order + 1 entries, the most clusters kept per order, default 0 = all
//      },
//      "pcce"      : {
//          "partition_size" : 2,       // spins per partition (required with mode pcce)
//          "trials"         : 30000,   // default 30000
//          "max_iterations" : 30000,   // default 30000
//          "kmeans_pp"      : true,    // default true
//          "print_trials"   : false,   // default false
//          "print_limit"    : 32,      // rank-zero trial details, 0 disables, default 32
//          "partition"      : "deterministic",
                                        // deterministic (default) or legacy
//          "connections"    : "centers",
                                        // centers (default) or spins
//          "cache"          : {"mode": "off", "file": "partition_cache.json"}
                                        // mode off (default), read or readwrite, with configuration.k_range the file
                                        // needs {k} (one cache per configuration)
//      },
//      "bundle"    : {
//          "max_dim"           : 64,   // largest bundle dimension (required)
//          "merge_cutoff"      : 0,    // default 0
//          "connections"       : "distance",
                                        // distance (default) or auto, without rdip and this key: checked by the run
//          "connection_budget" : 1000, // default: none
//          "max_clusters"      : 200000,
                                        // default 200000
//          "max_work"          : 20000000
                                        // default 20000000
//      },
//      "adaptive"  : {                 // docs/user/InputFile_Cluster.md. input, binary, cwd, mpiexec and ranks
                                        // are decided by the run, work_file is written by the driver into a step input
//          "output"                 : "./adaptive",
                                        // folder of the steps (required)
//          "probe_memory_gb"        : 0,
                                        // node memory of the probes running at once, default 0 (automatic)
//          "atol"                   : 1e-3,
                                        // absolute tolerance (required)
//          "rtol"                   : 0,
                                        // relative tolerance (required)
//          "policy"                 : "merge",
                                        // merge (the only policy, default)
//          "max_rounds"             : 8,
                                        // default 8
//          "max_patches"            : 64,
                                        // default 64
//          "max_probes"             : 4096,
                                        // default 4096
//          "max_candidates"         : 100,
                                        // default 16
//          "neighbors"              : 8,
                                        // default 8
//          "cache_mb"               : 256,
                                        // default 256
//          "timeout"                : 600,
                                        // seconds, default 600
//          "raw_cache"              : "raw_cache",
                                        // default none
//          "probe_batch_clusters"   : 0,
                                        // default 0
//          "merge_priority"         : "gain-per-cost",
                                        // gain-per-cost (default) or worst-error
//          "max_total_work"         : 8e12,
                                        // default 8000000000000.0
//          "prefer_physical"        : false,
                                        // switch, default false
//          "hub_candidates"         : false,
                                        // switch, default false
//          "regroup"                : false,
                                        // switch, default false
//          "regroup_max_regions"    : 64,
                                        // default 64
//          "regroup_max_candidates" : 256,
                                        // default 256
//          "regroup_factors"        : "probes",
                                        // probes (default) or factors
//          "regroup_whole"          : "probes"
                                        // probes (default) or except
//      },
//      "adaptive"  : {"work_file": "./work.json"},
                                        // (written by the driver into a step input, beside probe_memory_gb only)
//  },

    // ---- where and how it is computed -------------------------------------------------------------------------------
//  "execution"        : {
//      "backend"          : "cpu",     // cpu (default), cuda, auto
//      "gpu_options"      : {
//          "min_dim"               : 768,
                                        // positive, auto: smallest cluster dimension sent to the GPU, default 768
//          "packed_max_dim"        : 4,
                                        // 0 (off), 2 or 4, default 4 in main_gpu.out with cuda/auto and cCCE, else 0
//          "warp_max_dim"          : 27,
                                        // 0 (off) to 32, default 27 in main_gpu.out with cuda/auto and
                                        // cCCE/gCCE, else 0
//          "evolution_batch_size"  : 8,
                                        // default 8
//          "evolution_batch_steps" : 64,
                                        // default 64
//          "cluster_batch_mb"      : 64,
                                        // MiB per rank, default 64
//          "cache_mb"              : 64
                                        // nonnegative, MiB per rank, 0 disables, default 64
//      },
//      "propagator"       : "eigen",   // eigen (default), sparse (one state, cCCE, gCCE vector), expm (gCCE), or
                                        // {"method": "sparse", "min_dim": 4096, "tol": 1e-12} (min_dim, tol: sparse
                                        // only)
//      "evolution"        : "vector",  // gCCE: vector (default where it can run) or matrix
//      "evolution_update" : "recursive",
                                        // exact or recursive, default recursive where it exists, exact elsewhere
//      "gcce_refinement"  : {"max_sweeps": 32},
                                        // gCCE, off by default, writing the block turns it on ({"enabled": false}
                                        // keeps a written block off), max_sweeps default 32
//      "column_evolution" : {"min_dim": 256, "columns": 64},
                                        // one-state cCCE, on by default with these values, {"enabled": false} turns it
                                        // off
//      "profile"          : false      // per-dimension counts and timing, default false
//  },

    // ---- output -----------------------------------------------------------------------------------------------------
//  "output"           : {
//      "path"                  : "results",
                                        // the folder of the outputs, made when missing, -o DIR takes its place,
                                        // default: the working directory
//      "file"                  : "ccex_example",
                                        // output name, default: the input file name without .jsonc or .json (an
                                        // absolute name is refused with -o DIR or path)
//      "savemode"              : {"mode": "coherence", "nodiv": true},
                                        // or a string: coherence (default), coherence_tagged, cluster_coherence,
                                        // cluster_density_matrix, nodiv default true (writes _noDiv)
//      "cluster_contributions" : {"enabled": true, "max_mb": 256},
                                        // default off, max_mb default 256
//      "cluster_only"          : false // true: write <file>.clusters.json and stop before the dynamics (cCCE, gCCE,
                                        // and cluster.mode pcce plans with cCCE only. Not with cluster.mode bundle
                                        // without rdip (checked connections) or the adaptive probe, factors and except
                                        // work), default false
//  },

    // ---- automatic convergence --------------------------------------------------------------------------------------
//  "auto_convergence" : {              // main.out -f runs the automatic convergence
                                        // (docs/user/InputFile_AutoConvergence.md). Its runs use the mpirun -n N
                                        // main.out was started with. The bath configurations are the configuration
                                        // block's (with k_range those k, without it the k whose files exist)
//      "enabled"     : true,           // default true, false runs the ordinary calculation
//      "resume"      : "continue",     // continue (default) or fresh
//      "workdir"     : "auto",         // required
//      "time_window" : {"converge": true, "value": 1.0, "unit": "ms", "threshold": 0.01,
//                       "revival_check_factor": 2.5, "production_factor": 1.5},
                                        // unit s, ms (default), us, ns. Or "value": {"value": 1.0, "unit": "ms"} (unit)
                                        // without "unit" beside it
//      "rdip"        : {"converge": true, "value": 8, "candidates": [4, 6, 8], "num_candidates": 7, "threshold": 0.02},
//      "rbath"       : {"converge": true, "value": 15, "candidates": [10, 15], "num_candidates": 4, "threshold": 0.02},
//      "t2_fit"      : {"points": "peaks", "floor": 0.05, "n_range": [0.3, 4.0]},
                                        // points peaks (default) or all
//      "sampling"    : {
//          "test"        : [1, 2], "production": [3, 4],
                                        // default: found
//          "converge"    : true, "number": 10, "states": 5, "threshold": 0.05,
//          "test_number" : 10, "test_states": 5, "states_threshold": 0.05,
//          "max_states"  : 30, "min_configurations": 10, "max_configurations": 50
//      }
//  }
}
)CCEX";
