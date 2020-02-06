# Constant-Q-PSTD

These codes are used to simulate viscoacoustic wave propagation where the frequency-dependent attenuation and dispersion obey the frequency-independent Q model. The wave equation to be solved is a fractional Laplacian viscoacoustic wave equation. Two pseudo-spectral time-domain (PSTD) numerical schemes are involved in this package. The first one is used to solve the forward modeling problem and viscoacoustic common-shot gathers and wavefield snapshots can be generated. The second one is used to simulate wave propagation in inverse-Q media where the wavefield amplitudes increase exponentially. To avoid numerical instability caused by the boosting of high-frequency noise, a robust time-variant filter is incorporated into the PSTD scheme. The amplitude-boosted PSTD scheme can be used as wave propagation engines in attenuation-compensated reverse-time migration (Q-RTM) and viscoacoustic full waveform inversion (FWI).

All the codes are written with the standard C language under a linux system. The MPI/GPU parallel programming techniques are used. Therefore, one should first install MPICH and CUDA C compilers into your laptop or multi-node cluster. Of course, before that, one should download a compatible NVIDIA driver to support your NVIDIA graphics cards.

To compile the codes, one should 
