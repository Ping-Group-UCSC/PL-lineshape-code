NV-Diamond Example
===================================

What is computed here
-----------------------------------
In this example we calculate the photoluminescence of the negatively charge NV center in Diamond

Basic setup:
-----------------------------------
In this directory are two folders `relax-gs` and `relax-cdftup1` containing Quantum ESPRESSO input and output for running relaxation of NV- center diamond in its ground state and its first excited state, respectively. Within the `relax-gs` directory is another subdirectory `PH`, which contains the input and output of `ph.x` and `matdyn.x`. However, the PL-code only needs to read one file `*.mold` (here it is called `dynmat.mold`) in order to read the phonon modes and frequencies. To set up this calculation from scratch one would:

1. Run `pw.x` for relaxing ground state `relax-gs/relax.in`
2. Run `pw.x` for relaxing excited state `relax-cdftup1/relax.in`
3. Run `ph.x` for phonon calculation of ground state `relax-gs/PH/ph.in`
4. Run `dynmat.x` for applying acoustic sum rule `relax-gs/PH/dynmat.in`

Instructions:
-----------------------------------
0. In this folder all steps above have already been executed. Please return to the root of the example `Examples/NV-Diamond`.

1. Run the pl calculation. Note running the pl code without a specified input file results in several default settings. (Note that minimally one can simply execute `../../code/main.py`, the command below is good for running the calculation and wishing to save the output to file while also viewing the progress in the console.)

```bash
python3 -u ../../code/main.py | tee pl.out
```

2. The first thing the code will output is the parameters to be used in the calculation. Each of these parametrs can be changed in input file named `pl.in` (as done in step 10).

```
     No 'pl.in': Default parameters are assumed
         path_to_qe   = .
         zpl (eV)     = from input
         skfile       = None
         smear (eV)   = 0.006
         limit (s)    = 3.5e-13
         gamma (s^-1) = 2857142857100.0
         tolerance    = 1e-18
         hw_min (eV)  = 0.0
         hw_min (eV)  = 2.0
         hw_steps     = 300
```

3. Next (if not specified by the input) the code calculates the zero-phonon line (ZPL) reading the total energies from `relax-gs/relax.out` and `relax-cdfupt1/relax.out`.

```
Zero-Phonon Line calculated from QE output
     ZPL =   1.569642
```

4. Next the code calculates and plots S(hw), S(t), G(t). Refer to the cited work in the root repository `README.md` for explanation of each of these functions. The plots for each of these functions are automatically saved to png files `S_hw.png`, `S_t.png`, `G_t.png`, respectively. Also during this section the code prints some other useful information including the total Huang-Rhys factor (HR), a check that S(t->0) = HR, and G(t->inf) as shown below.

```
     Huang-Rhys factor =   2.199241             # total HR factor
     S(t->0) == HR?   True                      # this should always be the case
     G(t->inf) =   0.110887 +  -0.000000 i      # the limit of G(t) as t-> infinity
```

5. Finally the code begins the laborious task of integrating the generating function G(t) to calculate the final photoluminescence. This code uses an integration from the scipy package which will prints messages such as the one below until it reaches the convergence threshold set by `tolerance`.

```
Calculating A(ZPL - E)
     time now:  162.7176229953766  sec
/usr/local/lib/python3.7/site-packages/scipy/integrate/quadrature.py:802: AccuracyWarning: divmax (10) exceeded. Latest difference = 1.443711e-14
  AccuracyWarning)
/usr/local/lib/python3.7/site-packages/scipy/integrate/quadrature.py:802: AccuracyWarning: divmax (10) exceeded. Latest difference = 7.614841e-15
  AccuracyWarning)
/usr/local/lib/python3.7/site-packages/scipy/integrate/quadrature.py:802: AccuracyWarning: divmax (10) exceeded. Latest difference = 1.269757e-15
...
```

6. This integration will take approximately 45 minutes, and afterwards a final PL spectra will be produced with the data saved in `pl.dat` and plot saved to `pl.png`


7. Enter `Better` and run another calculation but this time specifying an input file `pl.in`. The code will automatically the code will read from its contents.

```bash
cd Better
python3 -u ../../../code/main.py | tee pl.out
```

-----------------------------------
Note from newer version:
The path_to_qe has been replaced by : "path_gs" and "path_ex"
- path_gs:
    
    The path of ground state calculation. 
    With qe interface it should be the header of scf or relax file. No “.in” or “.out” included
    
    e.g “./relax-gs/scf”, “./relax-gs/relax”
    
    With VASP interface it should be the directory path that include the POSCAR or CONTCAR
    
- path_ex: path if excited state
    
    same as path_gs
    
- phonon_interface:
    
    the default is qe . if phonopy then use phonopy for the phonon;
    (Notice we don’t need quote “” for this input)
    
- file_phonon:
    
    file path of the phonon file. 
    for qe it should be the paht of dynmat.mold file. e.g “./relax-gs/PH/dynmat.mold”
    for phonopy it should be the path of band.yaml file. e.g “./phonopy_calc/band.yaml”