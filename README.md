# SCAcode
A steady-state, 1-D, nuclear reactor thermal hydraulics single channel analysis code 
  
 ```
##############################################################################
####       #######      ####       ###########################################
##    ###   ####    #######   ###   #############################  ###########
#    ##########   ########   #####   ############################  ###########
#    #########   #########  #######  ############################  ###########
##    ########   #########  #######  ############################  ###########
####     #####   #########           ###       ###      ####       ###      ##
######    ####   #########    ###    ##   #######   ##   ##   ###  ##   ##   #
######    #####   ########   #####   ##  ########  ####  ##  ####  ##   ##   #
#   #    #######    ######   #####   ##   #######   ##   ##   ###  ###     ###
##     ###########      ##   #####   ###       ###      ####       ####      #
##############################################################################
 ```

## Introduction
SCAcode is intended to determine temperature, both radially throughout a fuel pin and of the moderator (water), and pressure loss within an ideal single channel within a nuclear reactor fuel assembly unit cell. It assumes steady state conditions and a nominal (cosine) flux shape throughout the core.  
  
The fundamental principles of this work can be found in: Nuclear Systems I: Thermal Hydraulic Fundamentals, N.E. Todreas and M.S. Kazimi, 2011 (2nd edition), (ISBN: 9781439808870).
  

## Getting Started
### Requirements
See also requirements.txt  
- numpy
- matplotlib
- pandas
- scipy

### Inputs
Input files must be python files. All variables within the provided input files (./inputs/inp_Benchmark_1_PWR_SCB_True.py and ./inputs/inp_Benchmark_2_BWR_SCB_True.py) must be defined. The values included are reasonable for each design. g and SB_constant should not be modified unless you are modeling a reactor on a different planet or in a different universe.  

### Outputs
Results will be exported to ./outputs/ into a new folder for each run. Data csv file and a single png plot are created every run. 

### Useage
#### If within current directory
```
python SCAcode_Setup.py -i <./local/path/to/input/file>.py  
```
#### If not
```
python SCAcode_Setup.py -i </entire/path/to/input/file>.py  
```

## Known Issues and Limitations
The following limitations apply based on physical limitation of the systems the code is intended to model (e.g., the fuel can not have melted for this model to be valid and water must boil in a boiling water reactor) and the limitations of mathematical equations/models used:  

| Type     | Limit     |
| :---     | :---      |
| Physical | Fuel melting |
| Physical | Cladding phase transition |
| Physical | Surpassing critical heat flux (CHF) |
| Physical | BWR-specific, boiling must occur ($0 < x_e < 1$) |
| Physical | PWR-specific, amount of permitted boiling ($x_e < 0$) |
| Physical | PWR-specific, power minimum ($q\prime \geq 0$) |
| Model    | Decoupled momentum and energy equations ($\Delta P_{tot}/P_{abs} \leq 10 \%$) |
| Model    | Cheng-Todreas, turbulent flow (Re $> 10^4$) required |
| Model    | Weisman, turbulent flow ($Re > 10^4$) required |
| Model    | HEM, mass flux upper limit |
| Model    | Schrock and Grossman, turbulent flow ($Re > 10^4$) required |
| Model    | Schrock and Grossman, $q\prime < CHF$ |
| Model    | Bowring, turbulent flow ($Re > 10^4$) required |
| Model    | Bowring, mass flux limit |


## Version History
- v1.0.0: initial release, 230512 BAT

## License
GNU GPL v3.0

## Authors
Brice Turner  
  
## Acknowledgments  
This work would not have been possible without the teachings of DuWayne Scubring, PhD, of the University of Florida.

