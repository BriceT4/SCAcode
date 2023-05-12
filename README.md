# SCAcode
A steady-state, 1-D, nuclear reactor thermal hydraulics single channel analysis code 


## Introduction
- Purpose of the project
SCAcode is intended to determine temperature, both radially throughout a fuel pin and of the moderator (water), and pressure loss within an ideal single channel within a nuclear reactor fuel assembly unit cell. It assumes steady state conditions and a nominal (cosine) flux shape throughout the core.  
  
The fundamental principles of this work can largely be found in: Nuclear Systems I: Thermal Hydraulic Fundamentals, N.E. Todreas and M.S. Kazimi, 2011 (2nd edition). (ISBN: 9781439808870)  
  

## Getting Started
### Requirements
See also requirements.txt  
- numpy
- matplotlib
- pandas
- scipy

### Inputs
Input files must be python files. The below variables must all be defined. Reasonable values for both PWR and BWR designs are provided in /inputs/inp_Benchmark_1_PWR_SCB_True.py and /inputs/inp_Benchmark_2_BWR_SCB_True.py, respectively. g and SB_constant should not be modified unless you are modeling a reactor on a different planet or in a different universe.  
- L: channel length (m)
- L_e: effective length (m)
- num_CV: number of control volumes along length.  

- D: (m), Diameter, converted from (mm)
- D_ci: (m), Diameter of cladding, inner, converted from (mm) 
- D_fo: (m), Diameter of fuel, outer, converted from (mm) 
- Pitch: (m), converted from (mm)  

- T_m_in: Moderator temperature at bottom of channel, (Celcius)
- P_nom: nominal pressure, (Pa)
- qp_max: Maximum q', (W/m) converted from (kW/m) 
- m_dot: moderator mass flow rate, (kg/s)
- k_cl: cladding conductivity (W/m-K)  

- g = 9.81, gravitational acceleration on Earth Z(m/s^2)
- SB_constant = 5.67E-8, Stefan-Boltzmann constant (W m^-2 K^-4)  

- Reactor_Type: must be 'PWR' or 'BWR'
- SCB_flag: must be True or False
- CHFR_crit_limi: must be > 1.0
- convergence_crit_gap: should be <= 0.001, but not much need to go smaller
- convergence_crit_fuel: should be <= 0.0001, but not much need to go smaller

### Outputs
Results will be exported to SCAcode/outputs/ into a new folder for each run. Data csv file and a single png plot are created every run. 

### Useage
#### If within current directory
```
python SCAcode_Setup.py -i &lt;.\path\to\input\file&gt;.py  
```
#### If not
```
python SCAcode_Setup.py -i &lt;\entire\path\to\input\file&gt;.py  
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
| Model    | Weisman, turbulent flow (Re $> 10^4$) required |
| Model    | HEM, mass flux upper limit |
| Model    | Schrock and Grossman, turbulent flow (Re $> 10^4$) required |
| Model    | Schrock and Grossman, $q\prime <$ CHF |
| Model    | Bowring, turbulent flow (Re $> 10^4$) required |
| Model    | Bowring, mass flux limit |


## Version History
N/A

## License
GNU GPL v.3.0

## Authors
Brice Turner  

