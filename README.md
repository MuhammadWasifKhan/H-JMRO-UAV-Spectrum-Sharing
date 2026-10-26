# H-JMRO: heterogeneous UAV spectrum sharing

MATLAB source code for the proposed H-JMRO scheme. The original two implementations have been separated into function files without changing their algorithm bodies.

## Requirements
- MATLAB (with `table`, `writetable` and Optimization Toolbox functionality as used by the original code)
- CVX installed and initialized with `cvx_setup`

## Running
Open MATLAB in the repository root and run either:

```matlab
ICC_Proposed_3UAV_Mode_Selection
ICC_HJMRO_Proposed_Only_40Iter_Tol1e3
```

Each entry point adds `shared/` and its experiment-specific folder to the MATLAB path. The two experiments should be run independently.

- `shared/`: functions identical in both supplied implementations.
- `flight_time/`: flight-time/mode-selection-specific functions.
- `convergence/`: convergence-specific functions.

The original scripts use 30 spatial drops and 300 Monte Carlo realizations and may take substantial time. They write results into the current MATLAB working directory.

**Validation:** Both MATLAB implementations were successfully executed using MATLAB and CVX. The flight-time, mode-selection, and convergence results were confirmed to reproduce the original implementations under the specified simulation settings. Numerical results were verified by the researcher.
