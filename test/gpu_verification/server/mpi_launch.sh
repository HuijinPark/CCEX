#!/usr/bin/env bash
set -euo pipefail
exec mpirun --bind-to core --map-by core "$@"

