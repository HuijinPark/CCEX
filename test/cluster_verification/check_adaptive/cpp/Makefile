# Tests for the refinement port (plan 0B). Test only: no MPI, no MKL.
# Tests run with THIS directory as the working directory: fixtures are ../reference/...
CXX ?= g++
BIN ?= bin
CXXFLAGS ?= -I../../../../external/sqlite -std=c++17 -O2 -g -Wall -Wextra -MMD -MP -I../../../../include/clustering/algo_adaptive
# Floating point stays exact. Two separate things are needed:
#   -ffp-contract=off   no fused multiply-add, on every compiler (design 3.3).
#                       -fp-model=precise alone does NOT disable FMA: Intel documents
#                       precise as fp-contract=on.
#   -fp-model=precise   Intel's default is fast, where std::isinf() is false and
#                       inf == inf can be false. On idm that stopped the non-finite
#                       strings from being written at all (measured 2026-09-16).
# The compiler is asked, not guessed from its name: CCEX builds through the mpiicpc
# wrapper with I_MPI_CXX=icpx, where the name says nothing (10th review).
# Intel Classic (icpc, __INTEL_COMPILER) as well as icpx: do_compile.sh asks the compiler whether
# isinf/isnan survive and adds -fp-model=precise for either (D15, 2026-10-02).
INTEL := $(shell $(CXX) -dM -E -x c++ /dev/null 2>/dev/null | grep -cE '__INTEL_LLVM_COMPILER|define __INTEL_COMPILER ')
#   -fmath-errno        CPython reads errno after a libm call and turns ERANGE into
#                       OverflowError - cmath.phase(1e300+1e-300j) raises rather
#                       than returning a tiny number. icpx defaults to not setting
#                       errno at all, which lost 112 of those cases on idm
#                       (measured 2026-09-17); gcc sets it by default, and asking
#                       for it explicitly says that the port depends on it.
FPFLAGS := -ffp-contract=off -fmath-errno $(if $(filter-out 0,$(INTEL)),-fp-model=precise,)
TESTS := test_counts test_binary_hash test_harness test_numeric test_scaled test_expansion test_merge test_manifest test_probes test_native test_errors test_int_division test_number_text test_complex test_fsum test_float_text test_probe_cache test_options test_json_write test_json_read test_json_config test_regroup
# Sources every test links against: everything in src/clustering/algo_adaptive but its entry point, and the expansion
# the two halves share (clustering/bundle/bundle_powers.cpp, integer code).
SHARED := ../../../../src/clustering/algo_adaptive/base_numeric.cpp ../../../../src/clustering/algo_adaptive/base_json.cpp ../../../../src/clustering/algo_adaptive/base_run.cpp ../../../../src/clustering/algo_adaptive/base_bundles.cpp ../../../../src/clustering/algo_adaptive/policy_merge.cpp ../../../../src/clustering/algo_adaptive/policy_regroup.cpp ../../../../src/clustering/algo_adaptive/options.cpp \
          ../../../../src/clustering/bundle/bundle_powers.cpp
SHARED_OBJ := $(addprefix $(BIN)/,$(notdir $(SHARED:.cpp=.o)))
PROGRAMS := $(addprefix $(BIN)/,$(TESTS))

all: $(PROGRAMS) $(BIN)/ccex-adaptive $(BIN)/regroup_replay
	@status=0; for t in $(PROGRAMS); do ./$$t || status=1; done; exit $$status

# The driver itself, built here so the same objects and flags are used. It links
# no MPI and no MKL: the binary it starts does that (design, structure C).
$(BIN)/ccex-adaptive: $(BIN)/adaptive_main.o $(SHARED_OBJ) $(BIN)/sqlite3.o
	$(CXX) $(CXXFLAGS) $(FPFLAGS) $^ -o $@ -lpthread -ldl

$(BIN)/adaptive_main.o: ../../../../src/clustering/algo_adaptive/adaptive_main.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) $(FPFLAGS) -c $< -o $@

# One object per source. Compiling several sources in one command writes ONE .d
# holding only the last source's headers, so a header change stops triggering a
# rebuild and a stale binary passes (measured 2026-09-16, 8th review).
# sqlite's amalgamation is C, and it is compiled as C: -x c keeps it on the path
# its authors support, and -w keeps its warnings out of ours. SQLITE_SANITIZE carries
# verify.sh sanitize's instrumentation here, since CXXFLAGS holds C++ flags (0ef89c R3).
SQLITE_SANITIZE ?=
$(BIN)/sqlite3.o: ../../../../external/sqlite/sqlite3.c | $(BIN)
	$(CXX) -x c -std=c11 -O2 -w $(SQLITE_SANITIZE) -c $< -o $@

$(BIN)/%.o: %.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) $(FPFLAGS) -c $< -o $@

$(BIN)/%.o: ../../../../src/clustering/algo_adaptive/%.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) $(FPFLAGS) -c $< -o $@

$(BIN)/bundle_powers.o: ../../../../src/clustering/bundle/bundle_powers.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) $(FPFLAGS) -c $< -o $@

$(BIN)/%: $(BIN)/%.o $(SHARED_OBJ) $(BIN)/sqlite3.o
	$(CXX) $(CXXFLAGS) $(FPFLAGS) $^ -o $@ -lpthread -ldl


$(BIN):
	mkdir -p $(BIN) out

.SECONDARY:            # keep the objects, or every build recompiles everything

clean:
	rm -rf $(BIN) out

-include $(wildcard $(BIN)/*.d)
