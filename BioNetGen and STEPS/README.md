#  BioNetGen translation of SBtab [Viswan_2018](https://github.com/jpgsantos/Model_Viswan_2018) model 
##### The model was translated from [SBtab](https://github.com/tlubitz/SBtab) model format to rule-based [BioNetGen language](http://bionetgen.org/) for the simulation with [STEPS](http://steps.sourceforge.net/STEPS/default.php) and [NFsim](http://michaelsneddon.net/nfsim/) solvers embedded in the [subcellular web app](https://subcellular.humanbrainproject.eu/ ) and with the [RuleBender](https://github.com/RuleWorld/rulebender)
#
### Conversion steps 

- Run ***convert_Viswan_2018_for_STEPS_optimised_from_SBTAB_to_SBML.R*** in [RStudio](https://www.rstudio.com/products/rstudio/download/) to translate SBtab model to SBML. This step depends on [SBtab to SBML converter](https://github.com/a-kramer/SBtabVFGEN)
- Run ***convert_Viswan_2018_for_STEPS_optimised_from_SBML_to_BNGL.ipynb*** [jupyter notebook](https://jupyter.org/) to translate SBML to BioNetGen language
- Import the BioNetGen model (***SBTAB_Viswan_2018_alternative.bngl***) to the **[subcellular web app](https://subcellular.humanbrainproject.eu/ )**. Add spine geometry ( ***.json***, ***.node***, ***.ele***, ***.face*** files) and stimulation pattern (***stim_E0.tsv***). See the **[subcellular web app help](https://humanbrainproject.github.io/hbp-sp6-guidebook/online_usecases/subcellular_level/subcellular_app/subcellular_app.html)** for details
- Simulate final model (***SBTAB_Viswan_2018_alternative.ebngl***) in the [subcellular web app](https://subcellular.humanbrainproject.eu/ ) using **STEPS** or **NFsim** solvers
- Simulate the BioNetGen model with the [RuleBender](https://github.com/RuleWorld/rulebender)

### Files and folders

***SBTAB_Viswan_2018_alternative.ebngl*** - extended BNGL model corresponding to the ***Viswan_2018_for_STEPS_optimised.xlsx*** SBTAB model with added geometry and stimulation patterns. Can be imported and simulated in the **[subcellular web app](https://subcellular.humanbrainproject.eu/ )**


***Viswan_2018_for_STEPS.xlsx*** - SBtab model equivalent to the original [Viswan_2018](https://github.com/jpgsantos/Model_Viswan_2018) model. It was obtained from ***[Viswan_2018.xlsx](https://github.com/jpgsantos/Model_Viswan_2018/blob/master/Viswan_2018.xlsx)*** by the modification of the model features incompatible with BNGL.


***Viswan_2018_for_STEPS_optimised.xlsx*** - SBtab model equivalent to the optimized [Viswan_2018](https://github.com/jpgsantos/Model_Viswan_2018) model. It was obtained from ***[Viswan_2018_optimized.xlsx](https://github.com/jpgsantos/Model_Viswan_2018/blob/master/Viswan_2018_optimized.xlsx)*** by the modification of the model features incompatible with BNGL.


***SBTAB_Viswan_2018.xml*** - SBML model translated from ***Viswan_2018_for_STEPS_optimised.xlsx*** model by ***convert_Viswan_2018_for_STEPS_optimised_from_SBTAB_to_SBML.R*** script based on [SBtab to SBML converter](https://github.com/a-kramer/SBtabVFGEN)


***SBTAB_Viswan_2018_alternative.bngl*** - BioNetGen model obtained from ***SBTAB_Viswan_2018.xml*** by ***convert_Viswan_2018_for_STEPS_optimised_from_SBML_to_BNGL.ipynb*** [jupyter notebook](https://jupyter.org/) which is based on ***sbml_to_bngl.py*** conversion tool 

***cell.ele***, ***cell.face***, ***cell.node***, ***cell.json*** - these files specify [TetGen](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) meshes and model geometry needed for the [subcellular web app](https://subcellular-bsp-epfl.apps.hbp.eu/model/meta/) *[Geometry](https://subcellular-bsp-epfl.apps.hbp.eu/model/geometry)* section and STEPS solver (see the **[subcellular web app help](https://humanbrainproject.github.io/hbp-sp6-guidebook/online_usecases/subcellular_level/subcellular_app/subcellular_app.html)** for details)


***stim_E0.tsv***, ***stim_E1.tsv*** - these files specify the stimulation pattern in *[Simulations](https://subcellular-bsp-epfl.apps.hbp.eu/model/simulations)* section of the [subcellular web app](https://subcellular-bsp-epfl.apps.hbp.eu/model/meta/) (corresponds to the experiments [E0](https://github.com/jpgsantos/Model_Viswan_2018/blob/master/tsv/Viswan_2018_optimized/E0.tsv) and [E1](https://github.com/jpgsantos/Model_Viswan_2018/blob/master/tsv/Viswan_2018_optimized/E1.tsv) of the SBtab model).


***Viswan_2018_alternative_RuleBender.bngl*** - BNGL model corresponding to the ***Viswan_2018_for_STEPS_optimised.xlsx*** SBtab model with additional section specifying stimulation and BioNetGen solver. Can be imported and simulated in the [RuleBender](https://github.com/RuleWorld/rulebender/releases/tag/RuleBender-2.3.2)

***SBtabVFGEN-master*** - the folder containing copy of [SBtab to SBML converter](https://github.com/a-kramer/SBtabVFGEN)

***SBTAB_Viswan_2018_for_STEPS_optimised*** - the folder containing **tsv** tables of ***Viswan_2018_for_STEPS_optimised.xlsx***

***sbml_to_bngl.py*** - the python tool for conversion of SBML models to BioNetGen language.

### Conversion of SBML to BioNetGen language

The conversion is implemented in ***sbml_to_bngl.py*** python module.
Two approaches are supported by *sbml_to_bngl.transform()* function:
- if *converter='pysb'* - the convertor based on the [Atomizer](https://ruleworld.github.io/atomizer/blog/basic/bng.html) implemented by pysb.importers.sbml.sbml_translator() function will be used. The [Atomizer](https://ruleworld.github.io/atomizer/blog/basic/bng.html) will try to modify the set of model molecules and reactions to convert them from reaction network to rule-based BioNetGen format. 
- *if converter='plain'* - a libsbml based convertor for sbml level 2, version 4 will be used. This converter produces a bngl approximation to reaction network format of a model. It is assumed that sbml models were obtained by importing a MATLAB simbiology model to sbml, or from SBTAB model by [SBtab to SBML converter](https://github.com/a-kramer/SBtabVFGEN).

Models expressed by SBML and SBTAB often are not fully competible with BNGL.
Additional model adaptation steps are required in this case to obtain a working BNGL model. 
These steps will be partially automatized by *sbml_to_bngl.transform()* function
if *adapt_steps* argument dictionary *'list_of_steps'* is nonempty.

The adaptation steps include:

1) The STEPS and NFsim solvers require different units for species quantities an kinetic rates. An adapted BNGL model provides modifed expressions for all species concentrations and kinetic rates and provides an easy way for units changing by specification of auxilary bngl model parameters: *Na* and *V_comparment_name*. These parameters should be selected to: *Na=6.022e23* and *V_comparment_name* = volume of corresponding compartment in liters for NFsim and to: *Na=1* and *V_comparment_name=1* for STEPS. 

2) Species with fixed concentrations are not supported by NFsim solver. The BNGL model adaptation will modify model reactions such that a fixed species concentration became a model parameter. This parameter can be used for the clamping of species concentration or for the stimulation pattern application

3) If SBML to BNGL convertor implemented in Atomizer is selected then additional transformation steps include renaming of duplicated molecule sites and reparing incorrect molecule names and kinetic rate transformations

4) In case when MATLAB simbiology is used for SBML model creation, adapt_steps will repare incorrect molecule and parameter names

5) Compartmental model of BNGL assumes tree structure of the set of model comartments. This assumption is often incompatible with mesh geometries supported by STEPS. The model adaption steps can produce bngl models with flexible compartmental structure

6) There is a number of incompatibilities between SBTAB and BNGL which still require manual correction. These include concentrations fixed to an expression, functional expressions for reaction rates in case of STEPS solver, some nonstandard types of reaction kinetic functions etc. The adapt_steps detects the cases of known incompatibilities and produces corresponding warning messages 
