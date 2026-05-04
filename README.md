# 1D-SIM — 1D McGill Sea Ice Model
 
1D-SIM is a simple one-dimensional Viscous-Plastic (VP) sea ice model designed to study sea ice dynamics in idealised cases.
 
## Background
 
The model has been used in several published studies, including:
 
- Auclair, J. P., Lemieux, J. F., Tremblay, L. B. & Ritchie, H. (2017). Implementation of Newton's method with an analytical Jacobian to solve the 1D sea ice momentum equation. *Journal of Computational Physics*, 340, 69–84.
- Lemieux, J.-F., Knoll, D. A., Losch, M. & Girard, C. (2014). A second-order accurate in time IMplicit–EXplicit (IMEX) integration scheme for sea ice dynamics. *Journal of Computational Physics*, 263, 375–392.

## About This Fork
 
This fork extends the original model with two main objectives:
 
1. **Well-posedness of the VP Rheology**: Investigate the well-posedness of the Viscous-Plastic Rheology and propose a new parametrization to address it.
2. **Granular-Compressible Rheology**: Develop a new Granular-Compressible Rheology and study its well-posedness.


# Running the Model
 
## Architecture
 
All model source files are located in the `src/` directory.
 
## Initialisation
 
Experiments and model parameters are configured in `ice.f90`. Output frequency is controlled via `out_step`: each entry specifies a time level at which the fields are saved to the output directory.
 
For example, with `Deltat=3600`, setting `out_step(1)=24` saves the fields after 1 day. Additional output times can be defined with `out_step(2)`, `out_step(3)`, and so on.
 
## Compilation
 
To compile the model, run:
 
```bash
make
```
 
This will create a `build/` directory containing all compiled `*.mod` files.
 
To clean the build and executable, run:
```bash
make clean
```

## Setup
 
Edit `run_1dsim.sh` to set the experiment number and target directory.

Change also the experiment number in `input_1dsim`. You can also choose to read the namelist `namelistSIM` or go with the default parameters. The option for restart can also be set to `.true.` in this file. 

The script will automatically create the experiment folder if it does not already exist, and copy the namelist into it.

 
## Running
 
To launch the simulation, you have two options:

 ### Option 1
```bash
./run_1dsim.sh
```
Which will create a directory for you in a specified path to store the experiment namelist and other files related to this experiment. 

### Option 2
```bash
./zoupa < input_1dsim > 'output_file' 
```
Which will run the model in the main directory and place the log file in it. 



