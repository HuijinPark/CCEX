"""auto_convergence setup driver for CCEX (plan: (internal notes)). This first version analyses
existing CCEX runs: the time window (time_window step) and the numbers of configurations and states (sampling step).
Submitting runs to a server is the next step."""
import sys as _sys

if _sys.version_info < (3, 10):     # the modules use X | None annotations (review 01ecda #4)
    raise SystemExit(f"autoconv needs Python 3.10 or newer with numpy, scipy and matplotlib (this is {_sys.version.split()[0]}); "
                     "set CCEX_PYTHON to such an interpreter")
from .config import ConfigError, Settings, load_config, parse_auto_convergence  # noqa: F401
from .curves import read_run, read_runs  # noqa: F401
from .sampling import allocate, bootstrap_norm_error, sampling_rms_proxy, variance_components  # noqa: F401
from .window import choose_window, pooled_mean  # noqa: F401
