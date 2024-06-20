# ωFQ
The frequency-dependent fluctuating charge model (ωFQ) for calculating plasmonic response of metal nanoparticles. Include both frequency-domain and time-domain implementations.

# Build
```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/wFQ Na515rod.xyz
```

# References
* 10.1039/C8NR09134J The original paper of ωFQ
* 10.3389/fchem.2020.00340 Comparison between ωFQ and continuum model
* 10.1021/acs.jpcc.1c04716 ωFQ applied to graphene
* 10.1021/acsphotonics.2c00761 The original paper of ωFQFμ
* 10.1021/acs.jctc.3c00177 Using QM/ωFQFμ to calculate Raman spectra
* 10.1039/D4NA00080C Using QM/ωFQFμ to calculate fluorescence spectra
* 10.1063/5.0205845 RT-TDDFT/TD-ωFQ method (my work)
* 10.1021/acs.jpclett.4c01337 Using RT-TDDFT/TD-ωFQ to study strong coupling (my work)
* 10.48550/arXiv.2406.10926 Real-time ωFQ and ωFQFμ
