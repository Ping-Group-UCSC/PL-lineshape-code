NV-Diamond Example
===================================

What is computed here
-----------------------------------
In this example we calculate the photoluminescence of the negatively charge NV center in Diamond; We use Quantum ESPRESSO to compute the ground state and excited state electronic structure. We use vasp as the force calculator of phonoy and compute phonon with phonopy.  (phonopy referece:https://phonopy.github.io/phonopy/index.html)

Basic setup:
-----------------------------------
In this directory are two folders `relax-gs` and `cdftup` containing Quantum ESPRESSO input and output for running relaxation of NV- center diamond in its ground state and its first excited state, respectively.The phonon is calculated in `phonon`, PL code only need the `band.yaml` in order to read the phonon modes and frequencies. To set up this calculation from scratch one would:

1. Run vasp for relaxing ground state 
2. Run vasp for relaxing excited state 
3. Go to the three directory in `phonon` to compute the phonon calculation:
    (a). Run `1_Gen_displacement.sh` in `0_Generate_inputs` to generate displaced geometries 
    (b). Run `0_Generate_vasp_inp.sh` in `1_Run_forces` to generate the scf inputs. Run all the inputs.
    (c). Run `0_Generate_Forceset.sh` and `1_Run_pp.sh` in `2_Calc_phonon` phonon
4. Run `local_run_pl.sh` in `Run_pl` to get PL result.

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
