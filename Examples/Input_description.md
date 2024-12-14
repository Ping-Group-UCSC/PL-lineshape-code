# Input description

input can be comment out.

### Basic input:

- path_gs:
    
    path of ground state calculation. 
    
    for qe it should be the header of scf or relax file. Not “.in” or “.out” included
    
    e.g “./relax-gs/scf”, “./relax-gs/relax”
    
    for VASP it should be the directory path that include the POSCAR or CONTCAR
    
- path_ex: path if excited state
    
    same as path_gs
    
- phonon_interface:
    
    the default is qe . if phonopy then read phonon from the phonopy input.
    
    Notice we don’t need quote “” for this input
    
    The phnonpy file support .yaml or .hdf5 format.
    
- file_phonon:
    
    file path of the phonon file. 
    
    for qe it should be : e.g “./relax-gs/PH/dynmat.mold”
    
    for phonopy it should be : e.g “./phonopy_calc/band.yaml” or “./band.hdf5”
    
- zpl:
    
    Zero phonon line in eV. If it is not provided then it will be calculated from path_gs and path_ex.
    

### Integral control

- skfile:
    
    file to read previously calculated wk, Sk. No need to use it. I don’t see a interface for this function. 
    
- smear:
    
    smearing for S(hw) for the gaussian function. default is 0.006 eV.
    
    - smear_end
        
        FFT method implemented the frequency dependent smearing. As it was discussed (Jin, Yu, Marco Govoni, et al, Giulia Galli.  *Physical Review Materials* 5, no. 8 (August 24, 2021): 084603. [https://doi.org/10.1103/PhysRevMaterials.5.084603](https://doi.org/10.1103/PhysRevMaterials.5.084603).) 
        
        the frequency dependent broadening of gaussian function is account for the continuum of vibrational modes.
        
- gamma:
    
    The gamma factor of PL with unit of s^-1(Hz). gamma represent the broadening PL, It is a free parameter in this model. In real situation this broadening have two contributions: 1. the homogeneous broadening due to anharmonic phonon interactions. 2.the inhomogeneous broadening due to ensemble averaging.
    
    In NV- case we choose 2.9e12=0.01eV ; a value around this order should be tested for a reasonable broadening. default=0.28571428571e14
    
    - alternatively, use `gamma_ev` to assign gamma in eV unit. this will directly relate to the peak broadening.
- hw_min and hw_max, hw_steps:
    
    hw_min and hw_max is the minimum and maximum energy of the PL spectrum; 
    
    hw_steps is the sampling steps. The unit is eV
    
    default = 1eV, 2.1eV, 600 steps.
    
- resolution:
    
    resolution = number of steps in 1eV; 
    
    The hw_steps is more prioritized, so in order to set the sampling with resolution, do not spedify hw_steps in the input.
    
- cal_method:
    
    default = ‘integral’; it can be either ‘integral’ or ‘FFT’. 
    
    This is the method used to compute the integral from S(hw) to S(t), and from G(t) to A(hw)
    
- integral_method = "quad_vec” or 'romberg’
    
    Parameter used when cal_method=”integral” .  It select the integral algorithm . default=’quad_vec’
    
- limit:
    
    Parameter used when cal_method=”integral” . determine the integral upper limit of integral of G(t) . unit of s. the should be large enough to approximate integral to infinite (larger than the region where G(t) >0) .  default=3.5e-13
    
- tolerance
    
    Parameter used when cal_method=”integral” . The error tolerance of integral. 1e-18 is default value, and 1e-23 should be well enough.
    
- Plot option:
    
    PL_hw_xmin, PL_hw_xmax : the plot range of PL, in eV
    
    S_hw_xmin, S_hw_xmax : the plot range of S(hw) in eV
    

### Output file:

- ‘'sk.dat’’ :
    
    1st column: $\hbar\omega_k$ , the energy of phonon mode k (meV)
    
    2nd column: The Sk value (partial HR’s factor contributed from phonon mode k)