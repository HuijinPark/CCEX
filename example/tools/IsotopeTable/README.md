# IsotopeTable

`make_isotope_table.py` writes `src/input/isotope_table.cpp`, the nuclear isotopes CCEX uses for a bath species the
gyrofile does not name (or for every species when the input has no gyrofile; `docs/user/InputFile_General.md#gyrofile`).
The source is the EasySpin isotope database already in the repository
(`example/tools/SpinTensorReader/spHamilReader/spin_database/isotopedata.txt`; magnetic moments from N. Stone,
IAEA INDC(NDS)-0658, 2014); γ = g_n μ_N / ħ in rad/ms/G (CODATA 2018), 154 isotopes (mass number > 0, positive spin; the database's element rows with mass number 0 and spin -1 are left out).

```sh
python3 example/tools/IsotopeTable/make_isotope_table.py      # from the repository root, after the database changes
```

Against the gyro files of the group: 13C 6.72828544 (file 6.728284, 2e-7 relative), 11B 8.58470641 (8.585), 10B 2.87467967
(2.875). An NV-13C cCCE2 run without a gyrofile differs from the one with `13C 0.5 6.728284` by 1.0e-5 in |L| (the gyro
difference only: the table's value written in a gyrofile gives 0).
