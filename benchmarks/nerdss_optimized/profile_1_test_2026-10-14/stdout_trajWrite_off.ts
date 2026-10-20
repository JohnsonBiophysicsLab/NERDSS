# t0_monotonic 841194.275889
    0.007268 Command: /private/tmp/claude-501/-Users-yue-Workspace-NERDSS/795ae25f-3942-4661-85be-f52029fde836/scratchpad/nerdss_opt_g -f parms_0.inp
    0.007625 Start date: 2026-10-14 13:45:00
    0.007639 RNG Seed: 1932554700
    0.007653 
    0.007662 Parsing Input:
    0.007668 This is a new simulation with input file: parms_0.inp
    0.007882 Parsing simulation parameters...
    0.007895 Read in nItr: 100000 timeSteps
    0.007900 Read in timeStep: 1 us
    0.007909 Read in timeWrite: 1000 timeSteps
    0.007915 Read in trajWrite: 1000000 timeSteps
    0.007920 Read in restartWrite: 2000000 timeSteps
    0.007928 Read in checkPoint: -1 timeSteps
    0.007933 Read in pdbWrite: -1 timeSteps
    0.007937 Parsing simulation boundary conditions...
    0.007945 Read in waterBox: [125.387 nm, 125.387 nm, 5470.63 nm]
    0.007950 Parsing mol file information...
    0.007957 lipid.mol
    0.008128 Read in isLipid: true
    0.008135 Read in checkOverlap: false
    0.008143 Read in D: [1um^2s^-1, 1um^2s^-1, 0um^2s^-1]
    0.008152 Read in Dr: [0rad^2s^-1, 0rad^2s^-1, 0rad^2s^-1]
    0.008171 Coordinates:
    0.008177 COM: [0nm, 0nm, 0nm]
    0.008185 head: [0nm, 0nm, 1nm]
    0.008190 Bonds:
    0.008198 com-head,
    0.008205 States:
    0.008210 head has 2 state(s): U	P
    0.008218 radius calculated from the ifacesCoord: 1 nm
    0.008223 mass auto calculated from the radius: 1
    0.008227 
    0.008234 adaptor.mol
    0.008349 Read in checkOverlap: false
    0.008359 Read in D: [25um^2s^-1, 25um^2s^-1, 25um^2s^-1]
    0.008366 Read in Dr: [0.5rad^2s^-1, 0.5rad^2s^-1, 0.5rad^2s^-1]
    0.008381 Coordinates:
    0.008386 COM: [0nm, 0nm, 0nm]
    0.008393 lbs: [0nm, 0nm, -2nm]
    0.008398 ebs: [0nm, 1nm, -1nm]
    0.008405 Bonds:
    0.008410 com-lbs, com-ebs, Invalid atom in adaptor mol file. Ignoring bonds...
    0.008421 radius calculated from the ifacesCoord: 2 nm
    0.008426 mass auto calculated from the radius: 2
    0.008430 
    0.008435 kinase.mol
    0.008541 Read in checkOverlap: false
    0.008550 Read in D: [25um^2s^-1, 25um^2s^-1, 25um^2s^-1]
    0.008558 Read in Dr: [0.5rad^2s^-1, 0.5rad^2s^-1, 0.5rad^2s^-1]
    0.008571 Coordinates:
    0.008576 COM: [0nm, 0nm, 0nm]
    0.008584 abs: [0nm, 0nm, 0.5nm]
    0.008588 lbs: [0nm, 0nm, -0.5nm]
    0.008596 Bonds:
    0.008600 com-abs, com-lbs,
    0.008608 radius calculated from the ifacesCoord: 0.5 nm
    0.008613 mass auto calculated from the radius: 0.5
    0.008618 
    0.008622 phosphatase.mol
    0.008786 Read in checkOverlap: false
    0.008796 Read in D: [25um^2s^-1, 25um^2s^-1, 25um^2s^-1]
    0.008803 Read in Dr: [0.5rad^2s^-1, 0.5rad^2s^-1, 0.5rad^2s^-1]
    0.008811 Coordinates:
    0.008819 COM: [0nm, 0nm, 0nm]
    0.008823 abs: [0nm, 0nm, 0.5nm]
    0.008830 lbs: [0nm, 0nm, -0.5nm]
    0.008838 Bonds:
    0.008843 com-abs, com-lbs,
    0.008850 radius calculated from the ifacesCoord: 0.5 nm
    0.008855 mass auto calculated from the radius: 0.5
    0.008859 
    0.008866 Parsing reactions...
    0.008877 Parsing new reaction with form:
    0.008882 	kinase(lbs) + lipid(head~U) <-> kinase(lbs!1).lipid(head~U!1)
    0.008887 Reversible reaction.
    0.008895 Bimolecular reaction detected.
    0.008900 Parsing reactants...
    0.008904 Parsing products...
    0.008909 Read in value of onRate3DMacro: 0.00137276uM^-1s^-1
    0.008916 Read in value of offRateMacro: 3.31104s^-1
    0.008921 Read in value of kcat: 0.0894988
    0.008925 Read in value of norm1: [0i + 0j + 1k]
    0.008932 Read in value of norm2: [0i + 0j + 1k]
    0.008998 Read in value of assocAngles: [theta1: 1.5708, theta2: 1.5708, phi1: nan, phi2: nan, omega: 3.14159]
    0.009006 Read in value of sigma: 1nm
    0.009013 Read in value of length3Dto2D: 1nm
    0.009018 Read in value of coupledRxn: label, PhosphorylationByKinase
    0.009025 Finished reading reaction parameters, constructing reaction.
    0.009030 Checking for valid states for reactant kinase
    0.009034 Checking for valid states for reactant lipid
    0.009039 Found state, looking for state change...
    0.009061 Found corresponding interfaces: React: kinase(lbs) -- Prod: kinase(lbs)
    0.009069 Found corresponding interfaces: React: lipid(head) -- Prod: lipid(head)
    0.009073 Allow Michaelis-Menten for this reaction. Probability of catalysis: kcat/(kcat+offRatekb): 0.0263186 total rate: 3.40059
    0.009078 Parameter set complete for reaction [kinase(lbs)+lipid(head~U)<->kinase(lbs!1).lipid(head~U!1)].
    0.009081 
    0.009089 Parsing new reaction with form:
    0.009093 	lipid(head~U) -> lipid(head~P)
    0.009097 Irreversible reaction.
    0.009101 UniMolStateChange reaction detected
    0.009105 Parsing reactants...
    0.009108 Parsing products...
    0.009112 Read in value of rate: 0
    0.009115 Read in value of rxnLabel: PhosphorylationByKinase
    0.009126 Finished reading reaction parameters, constructing reaction.
    0.009132 Checking for valid states for reactant lipid
    0.009138 Found state, looking for state change...
    0.009143 State change found on species 0. Interface head changes state from U to P.
    0.009147 Parameter set complete for reaction [lipid(head~U)->lipid(head~P)].
    0.009151 
    0.009154 Parsing new reaction with form:
    0.009157 	phosphatase(lbs) + lipid(head~P) <-> phosphatase(lbs!1).lipid(head~P!1)
    0.009161 Reversible reaction.
    0.009165 Bimolecular reaction detected.
    0.009169 Parsing reactants...
    0.009172 Parsing products...
    0.009176 Read in value of onRate3DMacro: 0.0691851uM^-1s^-1
    0.009185 Read in value of offRateMacro: 0.00443464s^-1
    0.009191 Read in value of kcat: 0.0276077
    0.009199 Read in value of norm1: [0i + 0j + 1k]
    0.009205 Read in value of norm2: [0i + 0j + 1k]
    0.009209 Read in value of assocAngles: [theta1: 1.5708, theta2: 1.5708, phi1: nan, phi2: nan, omega: 3.14159]
    0.009213 Read in value of sigma: 1nm
    0.009217 Read in value of length3Dto2D: 1nm
    0.009221 Read in value of coupledRxn: label, DephosphorylationByPhosphatase
    0.009232 Finished reading reaction parameters, constructing reaction.
    0.009236 Checking for valid states for reactant phosphatase
    0.009240 Checking for valid states for reactant lipid
    0.009244 Found state, looking for state change...
    0.009248 Found corresponding interfaces: React: phosphatase(lbs) -- Prod: phosphatase(lbs)
    0.009252 Found corresponding interfaces: React: lipid(head) -- Prod: lipid(head)
    0.009259 Allow Michaelis-Menten for this reaction. Probability of catalysis: kcat/(kcat+offRatekb): 0.861516 total rate: 0.0320455
    0.009263 Parameter set complete for reaction [phosphatase(lbs)+lipid(head~P)<->phosphatase(lbs!1).lipid(head~P!1)].
    0.009268 
    0.009272 Parsing new reaction with form:
    0.009276 	lipid(head~P) -> lipid(head~U)
    0.009284 Irreversible reaction.
    0.009288 UniMolStateChange reaction detected
    0.009292 Parsing reactants...
    0.009295 Parsing products...
    0.009298 Read in value of rate: 0
    0.009302 Read in value of rxnLabel: DephosphorylationByPhosphatase
    0.009309 Finished reading reaction parameters, constructing reaction.
    0.009313 Checking for valid states for reactant lipid
    0.009317 Found state, looking for state change...
    0.009320 State change found on species 0. Interface head changes state from P to U.
    0.009324 Parameter set complete for reaction [lipid(head~P)->lipid(head~U)].
    0.009328 
    0.009332 Parsing new reaction with form:
    0.009335 	lipid(head~P) + adaptor(lbs) <-> lipid(head~P!1).adaptor(lbs!1)
    0.009340 Reversible reaction.
    0.009345 Bimolecular reaction detected.
    0.009350 Parsing reactants...
    0.009354 Parsing products...
    0.009357 Read in value of onRate3DMacro: 1.16403uM^-1s^-1
    0.009361 Read in value of offRateMacro: 0.0865906s^-1
    0.009368 Read in value of norm1: [0i + 0j + 1k]
    0.009372 Read in value of norm2: [0i + 0j + 1k]
    0.009375 Read in value of assocAngles: [theta1: 1.5708, theta2: 1.5708, phi1: nan, phi2: nan, omega: 3.14159]
    0.009379 Read in value of length3Dto2D: 1nm
    0.009382 Read in value of sigma: 1nm
    0.009399 Finished reading reaction parameters, constructing reaction.
    0.009403 Checking for valid states for reactant lipid
    0.009407 Found state, looking for state change...
    0.009411 Checking for valid states for reactant adaptor
    0.009415 Found corresponding interfaces: React: lipid(head) -- Prod: lipid(head)
    0.009419 Found corresponding interfaces: React: adaptor(lbs) -- Prod: adaptor(lbs)
    0.009423 Parameter set complete for reaction [lipid(head~P)+adaptor(lbs)<->lipid(head~P!1).adaptor(lbs!1)].
    0.009426 
    0.009431 Parsing new reaction with form:
    0.009437 	kinase(abs) + adaptor(ebs,lbs!*) <-> kinase(abs!1).adaptor(ebs!1,lbs!*)
    0.009440 Reversible reaction.
    0.009444 Bimolecular reaction detected.
    0.009448 Parsing reactants...
    0.009451 Parsing products...
    0.009455 Read in value of onRate3DMacro: 5.49406uM^-1s^-1
    0.009459 Read in value of offRateMacro: 0.00612925s^-1
    0.009462 Read in value of norm1: [0i + 0j + 1k]
    0.009466 Read in value of norm2: [0i + 0j + 1k]
    0.009470 Read in value of assocAngles: [theta1: 1.5708, theta2: 2.35619, phi1: nan, phi2: 3.14159, omega: 3.14159]
    0.009476 Read in value of length3Dto2D: 1nm
    0.009480 Read in value of sigma: 1nm
    0.009483 Finished reading reaction parameters, constructing reaction.
    0.009487 Checking for valid states for reactant kinase
    0.009491 Checking for valid states for reactant adaptor
    0.009495 Found corresponding interfaces: React: kinase(abs) -- Prod: kinase(abs)
    0.009498 Found corresponding interfaces: React: adaptor(ebs) -- Prod: adaptor(ebs)
    0.009502 Parameter set complete for reaction [kinase(abs)+adaptor(ebs,lbs!*)<->kinase(abs!1).adaptor(ebs!1,lbs!*)].
    0.009506 
    0.009509 Parsing new reaction with form:
    0.009513 	phosphatase(abs) + adaptor(ebs,lbs!*) <-> phosphatase(abs!1).adaptor(ebs!1,lbs!*)
    0.009516 Reversible reaction.
    0.009520 Bimolecular reaction detected.
    0.009524 Parsing reactants...
    0.009527 Parsing products...
    0.009531 Read in value of onRate3DMacro: 0.153275uM^-1s^-1
    0.009534 Read in value of offRateMacro: 0.00108378s^-1
    0.009538 Read in value of norm1: [0i + 0j + 1k]
    0.009542 Read in value of norm2: [0i + 0j + 1k]
    0.009550 Read in value of assocAngles: [theta1: 1.5708, theta2: 2.35619, phi1: nan, phi2: 3.14159, omega: 3.14159]
    0.009553 Read in value of length3Dto2D: 1nm
    0.009557 Read in value of sigma: 1nm
    0.009560 Checking for valid states for reactant phosphatase
    0.009564 Checking for valid states for reactant adaptor
    0.009567 Found corresponding interfaces: React: phosphatase(abs) -- Prod: phosphatase(abs)
    0.009571 Found corresponding interfaces: React: adaptor(ebs) -- Prod: adaptor(ebs)
    0.009574 Parameter set complete for reaction [phosphatase(abs)+adaptor(ebs,lbs!*)<->phosphatase(abs!1).adaptor(ebs!1,lbs!*)].
    0.009579 
    0.009584 
    0.009609 Input file parsing complete
    0.009620 NUMBER OF MOLECULE TYPES: 4NUMBER OF INTERFACES PLUS STATES, including PRODUCTS: 13
    0.009962 NUMBER OF MOLECULES IN GEN COORDS: 895
    0.009970 Finding and fixing overlapping proteins.
    0.011632 Iteration: 1
    0.011639 Overlapping Proteins (false): 0
    0.011644 No overlapping proteins found.
    0.020570 
    0.020588 Partitioning simulation box into sub-boxes...
    0.020595 Rmaxlimit: 55.2036
    0.020668 Simulation volume parameters:
    0.020677 Total sub-volumes: 396
    0.020682 	Dimensions: [2, 2, 99]
    0.020686 	Maximum sub-volume neighbors: 13
    0.020690 	Sub-volume size: [62.6935, 62.6935, 55.2589]
    0.024790  Total species calculated from init_speciesFile.cpp: 13
    0.024808 Initialized copyNumSpecies. Printing all values. total species: 13
    0.024815  absolute index: 0 Copy nums: 278
    0.024822  absolute index: 1 Copy nums: 0
    0.024827  absolute index: 2 Copy nums: 100
    0.024839  absolute index: 3 Copy nums: 100
    0.024845  absolute index: 4 Copy nums: 296
    0.024849  absolute index: 5 Copy nums: 296
    0.024855  absolute index: 6 Copy nums: 221
    0.024865  absolute index: 7 Copy nums: 221
    0.024869  absolute index: 8 Copy nums: 0
    0.024876  absolute index: 9 Copy nums: 0
    0.024880  absolute index: 10 Copy nums: 0
    0.024885  absolute index: 11 Copy nums: 0
    0.024890  absolute index: 12 Copy nums: 0
    0.024901 Pro pair:  p1: 0 p2: 1 index: 1 lipid,adaptor'	Pro pair:  p1: 0 p2: 2 index: 2 lipid,kinase'	Pro pair:  p1: 0 p2: 3 index: 3 lipid,phosphatase'	Pro pair:  p1: 1 p2: 2 index: 6 adaptor,kinase'	Pro pair:  p1: 1 p2: 3 index: 7 adaptor,phosphatase'	size of proPairlist: 5 Number of pairs: 16 numMolTypes: 4
    0.025007 Spherical shell index: inactive.
    0.025018 
    0.025023 Simulation Parameters
    0.025028 Number of iterations: 100000 timesteps
    0.025036 Timestep: 1 us
    0.025041 Timestep log interval (timeWrite): 1000 timesteps
    0.025045 Restart file write interval: 2000000 timesteps
    0.025050 Coordinate write interval (trajWrite): 1000000 timesteps
    0.025054 PDB Coordinate write interval: -1 timesteps
    0.025058 Checkpoint write interval: 10000 timesteps
    0.025062 overlapSepLimit: 0.1 nm
    0.025067 Transition matrix write interval: 10000 timesteps
    0.025071 ClusterOverlapCheck: false
    0.025075 RNGwrite: false
    0.025078 bondedComplexWrite: -1
    0.025082 Restart file format: JSON
    0.025086 Numerical settings:
    0.025090   2D table integration abs/rel error: 1e-07 / 1e-07
    0.025103   2D table integration fallback/tail cutoff: 1e-06 / 1e-10
    0.025108   normalization integration abs/rel error: 1e-06 / 1e-06
    0.025112   table rate abs/rel tolerance: 1e-08 / 0
    0.025117   table diffusion abs/rel tolerance: 0.0001 / 0
    0.025121   explicit/implicit lipid flat-diffusion cutoff: 1e-10 / 1e-15
    0.025125   association same-angle abs/rel tolerance: 0.0001 / 0
    0.025129   association rotation/endpoint-sign tolerance: 1e-08 / 1e-11
    0.025133   Vec3D coordinate precision: 10000
    0.025138 Molecule specific parameters:
    0.025142 Number of unique molecule types: 4
    0.025146 Total number of unique interfaces and states, including product states: 13
    0.025150 Total number of complexes in system at start: 895
    0.025154 Total number of units (molecules + interfaces) in system at start: 2407
    0.025158 Maximum allowed number of unique 2D reactions: 10000
    0.025162  isSphere? false
    0.025167  sphere Radius 0
    0.025171  hasCompartment? false
    0.025175  compartment Radius 0
    0.025178  BOX geometry, dimensions:
    0.025182 125.387 125.387 5470.63
    0.025185  hasImplicitLipid? false
    0.025189 
    0.025192 Molecule Information
    0.025196 Molecule template 0
    0.025199 Name: lipid
    0.025203 Copy number:278
    0.025206  Diffusion trans: 1 1 0
    0.025220  Diffusion Rot: 0 0 0
    0.025225 Is a lipid: true
    0.025229 Is a promoter: false
    0.025233 Is a rod: true
    0.025237 Is a point: false
    0.025240 Radius: 1
    0.025244 
    0.025248 Interfaces:
    0.025252 Name: head
    0.025256 Relative index: 0
    0.025259 States:
    0.025263 (identity: U, absolute index: 0)
    0.025267 
    0.025271 (identity: P, absolute index: 1)
    0.025274 
    0.025278 Molecule template 1
    0.025282 Name: adaptor
    0.025286 Copy number:100
    0.025289  Diffusion trans: 25 25 25
    0.025293  Diffusion Rot: 0.5 0.5 0.5
    0.025296 Is a promoter: false
    0.025300 Is a rod: false
    0.025303 Is a point: false
    0.025307 Radius: 2
    0.025310 
    0.025313 Interfaces:
    0.025317 Name: lbs
    0.025320 Relative index: 0
    0.025324 Absolute index: 2
    0.025327 Forward Reactions:  (4)
    0.025331 
    0.025334 Name: ebs
    0.025337 Relative index: 1
    0.025346 Absolute index: 3
    0.025350 Forward Reactions:  (5) (6)
    0.025354 
    0.025357 Molecule template 2
    0.025361 Name: kinase
    0.025364 Copy number:296
    0.025368  Diffusion trans: 25 25 25
    0.025371  Diffusion Rot: 0.5 0.5 0.5
    0.025374 Is a promoter: false
    0.025378 Is a rod: true
    0.025381 Is a point: false
    0.025385 Radius: 0.5
    0.025388 
    0.025391 Interfaces:
    0.025401 Name: abs
    0.025404 Relative index: 0
    0.025407 Absolute index: 4
    0.025411 Forward Reactions:  (5)
    0.025416 
    0.025420 Name: lbs
    0.025424 Relative index: 1
    0.025427 Absolute index: 5
    0.025431 Forward Reactions:  (0)
    0.025434 
    0.025437 Molecule template 3
    0.025441 Name: phosphatase
    0.025444 Copy number:221
    0.025447  Diffusion trans: 25 25 25
    0.025455  Diffusion Rot: 0.5 0.5 0.5
    0.025459 Is a promoter: false
    0.025463 Is a rod: true
    0.025467 Is a point: false
    0.025471 Radius: 0.5
    0.025474 
    0.025478 Interfaces:
    0.025482 Name: abs
    0.025485 Relative index: 0
    0.025489 Absolute index: 6
    0.025493 Forward Reactions:  (6)
    0.025497 
    0.025500 Name: lbs
    0.025504 Relative index: 1
    0.025508 Absolute index: 7
    0.025511 Forward Reactions:  (2)
    0.025514 
    0.025518 
    0.025521 Reactions
    0.025525 Forward Rxn: 0
    0.025528 Absolute index: 0
    0.025531 Type: Bimolecular association
    0.025535 Reactants:
    0.025538  [iface name: lbs, iface index: 5, required Free]
    0.025542  [iface name: head, iface index: 0, requires state U, required Free]
    0.025546 
    0.025549 Products:
    0.025552  [iface name: lbs, iface index: 8, requires interaction]
    0.025556  [iface name: head, iface index: 8, requires state U, requires interaction]
    0.025559 
    0.025563 Rate(s):
    0.025566 Rate 0: 0.00227955
    0.025570 Reactant 1 requires interfaces:
    0.025573 Reactant 2 requires interfaces:
    0.025576 Sigma: 1
    0.025580 Association angles:
    0.025583 Theta 1    1.5708
    0.025586 Theta 2    1.5708
    0.025590 Phi 1         nan
    0.025593 Phi 2         nan
    0.025596 Omega     3.14159
    0.025600 
    0.025603 label: nonspecified
    0.025606 Is On Membrane: false
    0.025610 bindRadSameCom 1.1
    0.025613 loopCoopFactor 1
    0.025617 length3Dto2D 1
    0.025620 area3Dto1D 12.5664
    0.025623 isCoupled? true
    0.025631  coupledRxn Number: 2 type: Unimolecular state change prob to perform coupled: 0.0263186
    0.025635 
    0.025638 Back Rxn: 0
    0.025642 Absolute index: 1
    0.025646 Type: Bimolecular association
    0.025650 Reactants:
    0.025654  [iface name: lbs, iface index: 8, requires interaction]
    0.025657  [iface name: head, iface index: 8, requires state U, requires interaction]
    0.025661 
    0.025665 Products:
    0.025669  [iface name: lbs, iface index: 5, required Free]
    0.025673  [iface name: head, iface index: 0, requires state U, required Free]
    0.025677 
    0.025680 
    0.025684 Rate(s):
    0.025687 Rate 0: 3.40059
    0.025691 Reactant 1 requires interfaces:
    0.025694 Reactant 2 requires interfaces:
    0.025697 On Membrane? false
    0.025701 Forward Rxn: 1
    0.025704 Absolute index: 2
    0.025707 Type: Unimolecular state change
    0.025711 State Change Reactant: [iface name: head, iface index: 0, requires state U, required Free]
    0.025714 State Change Product: [iface name: head, iface index: 1, requires state P, required Free]
    0.025718 Rate(s):
    0.025721 Rate 0: 0
    0.025725 label: PhosphorylationByKinase
    0.025728 Forward Rxn: 2
    0.025731 Absolute index: 3
    0.025734 Type: Bimolecular association
    0.025738 Reactants:
    0.025741  [iface name: lbs, iface index: 7, required Free]
    0.025744  [iface name: head, iface index: 1, requires state P, required Free]
    0.025748 
    0.025751 Products:
    0.025758  [iface name: lbs, iface index: 9, requires interaction]
    0.025762  [iface name: head, iface index: 9, requires state P, requires interaction]
    0.025766 
    0.025770 Rate(s):
    0.025774 Rate 0: 0.114966
    0.025778 Reactant 1 requires interfaces:
    0.025781 Reactant 2 requires interfaces:
    0.025785 Sigma: 1
    0.025789 Association angles:
    0.025792 Theta 1    1.5708
    0.025796 Theta 2    1.5708
    0.025800 Phi 1         nan
    0.025803 Phi 2         nan
    0.025807 Omega     3.14159
    0.025811 
    0.025814 label: nonspecified
    0.025817 Is On Membrane: false
    0.025821 bindRadSameCom 1.1
    0.025828 loopCoopFactor 1
    0.025831 length3Dto2D 1
    0.025834 area3Dto1D 12.5664
    0.025838 isCoupled? true
    0.025841  coupledRxn Number: 5 type: Unimolecular state change prob to perform coupled: 0.861516
    0.025845 
    0.025848 Back Rxn: 2
    0.025851 Absolute index: 4
    0.025855 Type: Bimolecular association
    0.025858 Reactants:
    0.025861  [iface name: lbs, iface index: 9, requires interaction]
    0.025864  [iface name: head, iface index: 9, requires state P, requires interaction]
    0.025868 
    0.025871 Products:
    0.025874  [iface name: lbs, iface index: 7, required Free]
    0.025878  [iface name: head, iface index: 1, requires state P, required Free]
    0.025881 
    0.025884 
    0.025887 Rate(s):
    0.025891 Rate 0: 0.0320455
    0.025894 Reactant 1 requires interfaces:
    0.025897 Reactant 2 requires interfaces:
    0.025905 On Membrane? false
    0.025909 Forward Rxn: 3
    0.025912 Absolute index: 5
    0.025916 Type: Unimolecular state change
    0.025919 State Change Reactant: [iface name: head, iface index: 1, requires state P, required Free]
    0.025923 State Change Product: [iface name: head, iface index: 0, requires state U, required Free]
    0.025927 Rate(s):
    0.025931 Rate 0: 0
    0.025935 label: DephosphorylationByPhosphatase
    0.025939 Forward Rxn: 4
    0.025942 Absolute index: 6
    0.025946 Type: Bimolecular association
    0.025950 Reactants:
    0.025954  [iface name: head, iface index: 1, requires state P, required Free]
    0.025958  [iface name: lbs, iface index: 2, required Free]
    0.025961 
    0.025965 Products:
    0.025968  [iface name: head, iface index: 10, requires state P, requires interaction]
    0.025971  [iface name: lbs, iface index: 10, requires interaction]
    0.025975 
    0.025978 Rate(s):
    0.025981 Rate 0: 1.95636
    0.025985 Reactant 1 requires interfaces:
    0.025988 Reactant 2 requires interfaces:
    0.025991 Sigma: 1
    0.025995 Association angles:
    0.025998 Theta 1    1.5708
    0.026001 Theta 2    1.5708
    0.026004 Phi 1         nan
    0.026008 Phi 2         nan
    0.026011 Omega     3.14159
    0.026014 
    0.026018 label: nonspecified
    0.026021 Is On Membrane: false
    0.026024 bindRadSameCom 1.1
    0.026028 loopCoopFactor 1
    0.026031 length3Dto2D 1
    0.026034 area3Dto1D 12.5664
    0.026045 isCoupled? false
    0.026049 
    0.026052 Back Rxn: 4
    0.026056 Absolute index: 7
    0.026059 Type: Bimolecular association
    0.026063 Reactants:
    0.026067  [iface name: head, iface index: 10, requires state P, requires interaction]
    0.026071  [iface name: lbs, iface index: 10, requires interaction]
    0.026074 
    0.026078 Products:
    0.026082  [iface name: head, iface index: 1, requires state P, required Free]
    0.026086  [iface name: lbs, iface index: 2, required Free]
    0.026090 
    0.026093 
    0.026097 Rate(s):
    0.026100 Rate 0: 0.0876411
    0.026104 Reactant 1 requires interfaces:
    0.026107 Reactant 2 requires interfaces:
    0.026110 On Membrane? false
    0.026114 Forward Rxn: 5
    0.026117 Absolute index: 8
    0.026120 Type: Bimolecular association
    0.026124 Reactants:
    0.026127  [iface name: abs, iface index: 4, required Free]
    0.026130  [iface name: ebs, iface index: 3, required Free]
    0.026134 
    0.026137 Products:
    0.026141  [iface name: abs, iface index: 11, requires interaction]
    0.026144  [iface name: ebs, iface index: 11, requires interaction]
    0.026147 
    0.026150 Rate(s):
    0.026154 Rate 0: 9.25751
    0.026157 Reactant 1 requires interfaces:
    0.026160 Reactant 2 requires interfaces:
    0.026163  [iface name: lbs, iface index: 2, requires interaction]
    0.026167 Sigma: 1
    0.026170 Association angles:
    0.026173 Theta 1    1.5708
    0.026177 Theta 2   2.35619
    0.026184 Phi 1         nan
    0.026188 Phi 2     3.14159
    0.026191 Omega     3.14159
    0.026195 
    0.026199 label: nonspecified
    0.026202 Is On Membrane: false
    0.026206 bindRadSameCom 1.1
    0.026210 loopCoopFactor 1
    0.026214 length3Dto2D 1
    0.026221 area3Dto1D 12.5664
    0.026225 isCoupled? false
    0.026228 
    0.026232 Back Rxn: 5
    0.026236 Absolute index: 9
    0.026239 Type: Bimolecular association
    0.026244 Reactants:
    0.026248  [iface name: abs, iface index: 11, requires interaction]
    0.026274  [iface name: ebs, iface index: 11, requires interaction]
    0.026282 
    0.026287 Products:
    0.026294  [iface name: abs, iface index: 4, required Free]
    0.026298  [iface name: ebs, iface index: 3, required Free]
    0.026303 
    0.026307 
    0.026311 Rate(s):
    0.026315 Rate 0: 0.00621957
    0.026319 Reactant 1 requires interfaces:
    0.026324 Reactant 2 requires interfaces:
    0.026330  [iface name: lbs, iface index: 2, requires interaction]
    0.026335 On Membrane? false
    0.026341 Forward Rxn: 6
    0.026347 Absolute index: 10
    0.026353 Type: Bimolecular association
    0.026357 Reactants:
    0.026361  [iface name: abs, iface index: 6, required Free]
    0.026364  [iface name: ebs, iface index: 3, required Free]
    0.026368 
    0.026371 Products:
    0.026374  [iface name: abs, iface index: 12, requires interaction]
    0.026378  [iface name: ebs, iface index: 12, requires interaction]
    0.026381 
    0.026385 Rate(s):
    0.026388 Rate 0: 0.254622
    0.026399 Reactant 1 requires interfaces:
    0.026404 Reactant 2 requires interfaces:
    0.026438  [iface name: lbs, iface index: 2, requires interaction]
    0.026492 Sigma: 1
    0.026510 Association angles:
    0.026517 Theta 1    1.5708
    0.026523 Theta 2   2.35619
    0.026528 Phi 1         nan
    0.026534 Phi 2     3.14159
    0.026539 Omega     3.14159
    0.026553 
    0.026560 label: nonspecified
    0.026566 Is On Membrane: false
    0.026571 bindRadSameCom 1.1
    0.026577 loopCoopFactor 1
    0.026582 length3Dto2D 1
    0.026587 area3Dto1D 12.5664
    0.026593 isCoupled? false
    0.026600 
    0.026607 Back Rxn: 6
    0.026613 Absolute index: 11
    0.026618 Type: Bimolecular association
    0.026624 Reactants:
    0.026628  [iface name: abs, iface index: 12, requires interaction]
    0.026634  [iface name: ebs, iface index: 12, requires interaction]
    0.026640 
    0.026645 Products:
    0.026651  [iface name: abs, iface index: 6, required Free]
    0.026656  [iface name: ebs, iface index: 3, required Free]
    0.026662 
    0.026669 
    0.026676 Rate(s):
    0.026682 Rate 0: 0.00108422
    0.026688 Reactant 1 requires interfaces:
    0.026693 Reactant 2 requires interfaces:
    0.026699  [iface name: lbs, iface index: 2, requires interaction]
    0.026704 On Membrane? false
    0.033831 *************** BEGIN SIMULATION ****************
    0.353506 -------------------
    0.353535 End iteration: 1000, simulation time: 1.000000e-03 seconds.
    0.353631 System time: 2026-10-14 13:45:00
    0.353639 Elapsed time: 0 minutes
    0.353645 Number of molecules: 895
    0.353649 Number of complexes: 895
    0.353656 Total reaction matches: 142490
    0.353664 Avg timestep duration: 3.195960e-04, iterations remaining: 99000, Time left: 0 minutes
    0.353669 Estimated end time: 2026-10-14 13:45:32
    0.353674 -------------------------------------------------
    0.695900 -------------------
    0.695929 End iteration: 2000, simulation time: 2.000000e-03 seconds.
    0.696019 System time: 2026-10-14 13:45:00
    0.696028 Elapsed time: 0 minutes
    0.696034 Number of molecules: 895
    0.696039 Number of complexes: 895
    0.696044 Total reaction matches: 398710
    0.696053 Avg timestep duration: 3.421750e-04, iterations remaining: 98000, Time left: 0 minutes
    0.696058 Estimated end time: 2026-10-14 13:45:34
    0.696062 -------------------------------------------------
    1.019787 -------------------
    1.019826 End iteration: 3000, simulation time: 3.000000e-03 seconds.
    1.019914 System time: 2026-10-14 13:45:01
    1.019925 Elapsed time: 0 minutes
    1.019930 Number of molecules: 895
    1.019934 Number of complexes: 895
    1.019939 Total reaction matches: 615488
    1.019948 Avg timestep duration: 3.236690e-04, iterations remaining: 97000, Time left: 0 minutes
    1.019962 Estimated end time: 2026-10-14 13:45:32
    1.019967 -------------------------------------------------
    1.308232 -------------------
    1.308261 End iteration: 4000, simulation time: 4.000000e-03 seconds.
    1.308366 System time: 2026-10-14 13:45:01
    1.308375 Elapsed time: 0 minutes
    1.308380 Number of molecules: 895
    1.308385 Number of complexes: 895
    1.308390 Total reaction matches: 733956
    1.308398 Avg timestep duration: 2.882240e-04, iterations remaining: 96000, Time left: 0 minutes
    1.308403 Estimated end time: 2026-10-14 13:45:29
    1.308407 -------------------------------------------------
    1.669888 -------------------
    1.669915 End iteration: 5000, simulation time: 5.000000e-03 seconds.
    1.670026 System time: 2026-10-14 13:45:01
    1.670036 Elapsed time: 0 minutes
    1.670041 Number of molecules: 895
    1.670046 Number of complexes: 895
    1.670050 Total reaction matches: 1057573
    1.670059 Avg timestep duration: 3.614250e-04, iterations remaining: 95000, Time left: 0 minutes
    1.670064 Estimated end time: 2026-10-14 13:45:36
    1.670069 -------------------------------------------------
    1.978666 -------------------
    1.978691 End iteration: 6000, simulation time: 6.000000e-03 seconds.
    1.978797 System time: 2026-10-14 13:45:02
    1.978806 Elapsed time: 0 minutes
    1.978811 Number of molecules: 895
    1.978815 Number of complexes: 895
    1.978820 Total reaction matches: 1242434
    1.978828 Avg timestep duration: 3.085540e-04, iterations remaining: 94000, Time left: 0 minutes
    1.978833 Estimated end time: 2026-10-14 13:45:31
    1.978838 -------------------------------------------------
    2.299043 -------------------
    2.299081 End iteration: 7000, simulation time: 7.000000e-03 seconds.
    2.299159 System time: 2026-10-14 13:45:02
    2.299169 Elapsed time: 0 minutes
    2.299174 Number of molecules: 895
    2.299178 Number of complexes: 895
    2.299183 Total reaction matches: 1498500
    2.299192 Avg timestep duration: 3.201350e-04, iterations remaining: 93000, Time left: 0 minutes
    2.299197 Estimated end time: 2026-10-14 13:45:32
    2.299202 -------------------------------------------------
    2.630564 -------------------
    2.630593 End iteration: 8000, simulation time: 8.000000e-03 seconds.
    2.630698 System time: 2026-10-14 13:45:02
    2.630708 Elapsed time: 0 minutes
    2.630713 Number of molecules: 895
    2.630718 Number of complexes: 895
    2.630722 Total reaction matches: 1785830
    2.630730 Avg timestep duration: 3.312990e-04, iterations remaining: 92000, Time left: 0 minutes
    2.630735 Estimated end time: 2026-10-14 13:45:33
    2.630742 -------------------------------------------------
    3.020544 -------------------
    3.020571 End iteration: 9000, simulation time: 9.000000e-03 seconds.
    3.020659 System time: 2026-10-14 13:45:03
    3.020669 Elapsed time: 0 minutes
    3.020674 Number of molecules: 895
    3.020679 Number of complexes: 895
    3.020683 Total reaction matches: 2256722
    3.020695 Avg timestep duration: 3.897420e-04, iterations remaining: 91000, Time left: 0 minutes
    3.020701 Estimated end time: 2026-10-14 13:45:38
    3.020705 -------------------------------------------------
    3.373782 -------------------
    3.373819 End iteration: 10000, simulation time: 1.000000e-02 seconds.
    3.373916 System time: 2026-10-14 13:45:03
    3.373930 Elapsed time: 0 minutes
    3.373937 Number of molecules: 895
    3.373941 Number of complexes: 895
    3.373946 Total reaction matches: 2561147
    3.373954 Avg timestep duration: 3.529560e-04, iterations remaining: 90000, Time left: 0 minutes
    3.373959 Estimated end time: 2026-10-14 13:45:35
    3.373963 -------------------------------------------------
    3.717251 -------------------
    3.717305 End iteration: 11000, simulation time: 1.100000e-02 seconds.
    3.717363 System time: 2026-10-14 13:45:03
    3.717373 Elapsed time: 0 minutes
    3.717378 Number of molecules: 895
    3.717383 Number of complexes: 895
    3.717398 Total reaction matches: 2862086
    3.717407 Avg timestep duration: 3.432320e-04, iterations remaining: 89000, Time left: 0 minutes
    3.717415 Estimated end time: 2026-10-14 13:45:34
    3.717419 -------------------------------------------------
    4.077049 -------------------
    4.077076 End iteration: 12000, simulation time: 1.200000e-02 seconds.
    4.077144 System time: 2026-10-14 13:45:04
    4.077154 Elapsed time: 0 minutes
    4.077159 Number of molecules: 895
    4.077164 Number of complexes: 895
    4.077168 Total reaction matches: 3215341
    4.077177 Avg timestep duration: 3.595770e-04, iterations remaining: 88000, Time left: 0 minutes
    4.077182 Estimated end time: 2026-10-14 13:45:35
    4.077187 -------------------------------------------------
    4.465473 -------------------
    4.465504 End iteration: 13000, simulation time: 1.300000e-02 seconds.
    4.465575 System time: 2026-10-14 13:45:04
    4.465585 Elapsed time: 0 minutes
    4.465592 Number of molecules: 895
    4.465598 Number of complexes: 895
    4.465604 Total reaction matches: 3596015
    4.465613 Avg timestep duration: 3.882400e-04, iterations remaining: 87000, Time left: 0 minutes
    4.465620 Estimated end time: 2026-10-14 13:45:38
    4.465626 -------------------------------------------------
    4.852458 -------------------
    4.852496 End iteration: 14000, simulation time: 1.400000e-02 seconds.
    4.852581 System time: 2026-10-14 13:45:04
    4.852595 Elapsed time: 0 minutes
    4.852601 Number of molecules: 895
    4.852607 Number of complexes: 895
    4.852612 Total reaction matches: 3968292
    4.852622 Avg timestep duration: 3.867730e-04, iterations remaining: 86000, Time left: 0 minutes
    4.852627 Estimated end time: 2026-10-14 13:45:38
    4.852632 -------------------------------------------------
    5.164562 -------------------
    5.164590 End iteration: 15000, simulation time: 1.500000e-02 seconds.
    5.164691 System time: 2026-10-14 13:45:05
    5.164700 Elapsed time: 0 minutes
    5.164706 Number of molecules: 895
    5.164711 Number of complexes: 895
    5.164716 Total reaction matches: 4121603
    5.164724 Avg timestep duration: 3.118870e-04, iterations remaining: 85000, Time left: 0 minutes
    5.164729 Estimated end time: 2026-10-14 13:45:31
    5.164734 -------------------------------------------------
    5.455981 -------------------
    5.456013 End iteration: 16000, simulation time: 1.600000e-02 seconds.
    5.456126 System time: 2026-10-14 13:45:05
    5.456136 Elapsed time: 0 minutes
    5.456142 Number of molecules: 895
    5.456148 Number of complexes: 895
    5.456152 Total reaction matches: 4218734
    5.456161 Avg timestep duration: 2.911840e-04, iterations remaining: 84000, Time left: 0 minutes
    5.456166 Estimated end time: 2026-10-14 13:45:30
    5.456171 -------------------------------------------------
    5.789361 -------------------
    5.789401 End iteration: 17000, simulation time: 1.700000e-02 seconds.
    5.789485 System time: 2026-10-14 13:45:05
    5.789496 Elapsed time: 0 minutes
    5.789502 Number of molecules: 895
    5.789507 Number of complexes: 895
    5.789511 Total reaction matches: 4444432
    5.789520 Avg timestep duration: 3.331260e-04, iterations remaining: 83000, Time left: 0 minutes
    5.789525 Estimated end time: 2026-10-14 13:45:33
    5.789530 -------------------------------------------------
    6.149150 -------------------
    6.149180 End iteration: 18000, simulation time: 1.800000e-02 seconds.
    6.149285 System time: 2026-10-14 13:45:06
    6.149295 Elapsed time: 0 minutes
    6.149300 Number of molecules: 895
    6.149305 Number of complexes: 895
    6.149309 Total reaction matches: 4766692
    6.149318 Avg timestep duration: 3.595660e-04, iterations remaining: 82000, Time left: 0 minutes
    6.149323 Estimated end time: 2026-10-14 13:45:35
    6.149327 -------------------------------------------------
    6.489015 -------------------
    6.489047 End iteration: 19000, simulation time: 1.900000e-02 seconds.
    6.489129 System time: 2026-10-14 13:45:06
    6.489152 Elapsed time: 0 minutes
    6.489158 Number of molecules: 895
    6.489162 Number of complexes: 895
    6.489167 Total reaction matches: 5010556
    6.489176 Avg timestep duration: 3.396010e-04, iterations remaining: 81000, Time left: 0 minutes
    6.489181 Estimated end time: 2026-10-14 13:45:34
    6.489185 -------------------------------------------------
    6.841228 -------------------
    6.841265 End iteration: 20000, simulation time: 2.000000e-02 seconds.
    6.841410 System time: 2026-10-14 13:45:06
    6.841419 Elapsed time: 0 minutes
    6.841424 Number of molecules: 895
    6.841429 Number of complexes: 895
    6.841433 Total reaction matches: 5230967
    6.841442 Avg timestep duration: 3.519430e-04, iterations remaining: 80000, Time left: 0 minutes
    6.841448 Estimated end time: 2026-10-14 13:45:35
    6.841452 -------------------------------------------------
    7.214701 -------------------
    7.214737 End iteration: 21000, simulation time: 2.100000e-02 seconds.
    7.214842 System time: 2026-10-14 13:45:07
    7.214852 Elapsed time: 0 minutes
    7.214858 Number of molecules: 895
    7.214862 Number of complexes: 895
    7.214866 Total reaction matches: 5576952
    7.214874 Avg timestep duration: 3.731900e-04, iterations remaining: 79000, Time left: 0 minutes
    7.214879 Estimated end time: 2026-10-14 13:45:36
    7.214886 -------------------------------------------------
    7.607956 -------------------
    7.607983 End iteration: 22000, simulation time: 2.200000e-02 seconds.
    7.608063 System time: 2026-10-14 13:45:07
    7.608072 Elapsed time: 0 minutes
    7.608077 Number of molecules: 895
    7.608081 Number of complexes: 895
    7.608085 Total reaction matches: 6004454
    7.608093 Avg timestep duration: 3.930000e-04, iterations remaining: 78000, Time left: 0 minutes
    7.608098 Estimated end time: 2026-10-14 13:45:38
    7.608102 -------------------------------------------------
    7.967787 -------------------
    7.967814 End iteration: 23000, simulation time: 2.300000e-02 seconds.
    7.968017 System time: 2026-10-14 13:45:08
    7.968058 Elapsed time: 0 minutes
    7.968100 Number of molecules: 895
    7.968111 Number of complexes: 895
    7.968117 Total reaction matches: 6328055
    7.968123 Avg timestep duration: 3.596360e-04, iterations remaining: 77000, Time left: 0 minutes
    7.968141 Estimated end time: 2026-10-14 13:45:35
    7.968147 -------------------------------------------------
    8.317529 -------------------
    8.317563 End iteration: 24000, simulation time: 2.400000e-02 seconds.
    8.317659 System time: 2026-10-14 13:45:08
    8.317669 Elapsed time: 0 minutes
    8.317675 Number of molecules: 895
    8.317680 Number of complexes: 895
    8.317685 Total reaction matches: 6620848
    8.317694 Avg timestep duration: 3.493820e-04, iterations remaining: 76000, Time left: 0 minutes
    8.317700 Estimated end time: 2026-10-14 13:45:34
    8.317705 -------------------------------------------------
    8.644517 -------------------
    8.644562 End iteration: 25000, simulation time: 2.500000e-02 seconds.
    8.644658 System time: 2026-10-14 13:45:08
    8.644669 Elapsed time: 0 minutes
    8.644674 Number of molecules: 895
    8.644679 Number of complexes: 895
    8.644683 Total reaction matches: 6846817
    8.644692 Avg timestep duration: 3.267220e-04, iterations remaining: 75000, Time left: 0 minutes
    8.644697 Estimated end time: 2026-10-14 13:45:33
    8.644701 -------------------------------------------------
    8.986926 -------------------
    8.986959 End iteration: 26000, simulation time: 2.600000e-02 seconds.
    8.987046 System time: 2026-10-14 13:45:09
    8.987056 Elapsed time: 0 minutes
    8.987062 Number of molecules: 895
    8.987066 Number of complexes: 895
    8.987071 Total reaction matches: 7080169
    8.987080 Avg timestep duration: 3.421740e-04, iterations remaining: 74000, Time left: 0 minutes
    8.987085 Estimated end time: 2026-10-14 13:45:34
    8.987090 -------------------------------------------------
    9.306608 -------------------
    9.306656 End iteration: 27000, simulation time: 2.700000e-02 seconds.
    9.306742 System time: 2026-10-14 13:45:09
    9.306752 Elapsed time: 0 minutes
    9.306757 Number of molecules: 895
    9.306761 Number of complexes: 895
    9.306766 Total reaction matches: 7269723
    9.306774 Avg timestep duration: 3.194560e-04, iterations remaining: 73000, Time left: 0 minutes
    9.306780 Estimated end time: 2026-10-14 13:45:32
    9.306784 -------------------------------------------------
    9.643858 -------------------
    9.643892 End iteration: 28000, simulation time: 2.800000e-02 seconds.
    9.643984 System time: 2026-10-14 13:45:09
    9.643995 Elapsed time: 0 minutes
    9.644000 Number of molecules: 895
    9.644005 Number of complexes: 895
    9.644010 Total reaction matches: 7579100
    9.644018 Avg timestep duration: 3.370080e-04, iterations remaining: 72000, Time left: 0 minutes
    9.644024 Estimated end time: 2026-10-14 13:45:34
    9.644029 -------------------------------------------------
    9.980906 -------------------
    9.980931 End iteration: 29000, simulation time: 2.900000e-02 seconds.
    9.981023 System time: 2026-10-14 13:45:10
    9.981032 Elapsed time: 0 minutes
    9.981037 Number of molecules: 895
    9.981042 Number of complexes: 895
    9.981047 Total reaction matches: 7857377
    9.981055 Avg timestep duration: 3.368280e-04, iterations remaining: 71000, Time left: 0 minutes
    9.981060 Estimated end time: 2026-10-14 13:45:34
    9.981064 -------------------------------------------------
   10.289439 -------------------
   10.289473 End iteration: 30000, simulation time: 3.000000e-02 seconds.
   10.289572 System time: 2026-10-14 13:45:10
   10.289581 Elapsed time: 0 minutes
   10.289587 Number of molecules: 895
   10.289591 Number of complexes: 895
   10.289595 Total reaction matches: 7998578
   10.289604 Avg timestep duration: 3.083170e-04, iterations remaining: 70000, Time left: 0 minutes
   10.289609 Estimated end time: 2026-10-14 13:45:31
   10.289613 -------------------------------------------------
   10.584916 -------------------
   10.584961 End iteration: 31000, simulation time: 3.100000e-02 seconds.
   10.585076 System time: 2026-10-14 13:45:10
   10.585087 Elapsed time: 0 minutes
   10.585094 Number of molecules: 895
   10.585100 Number of complexes: 895
   10.585106 Total reaction matches: 8150711
   10.585116 Avg timestep duration: 2.952450e-04, iterations remaining: 69000, Time left: 0 minutes
   10.585123 Estimated end time: 2026-10-14 13:45:31
   10.585129 -------------------------------------------------
   10.895939 -------------------
   10.895970 End iteration: 32000, simulation time: 3.200000e-02 seconds.
   10.896075 System time: 2026-10-14 13:45:11
   10.896087 Elapsed time: 0 minutes
   10.896093 Number of molecules: 895
   10.896099 Number of complexes: 895
   10.896104 Total reaction matches: 8347036
   10.896113 Avg timestep duration: 3.107370e-04, iterations remaining: 68000, Time left: 0 minutes
   10.896120 Estimated end time: 2026-10-14 13:45:32
   10.896125 -------------------------------------------------
   11.214765 -------------------
   11.214797 End iteration: 33000, simulation time: 3.300000e-02 seconds.
   11.214905 System time: 2026-10-14 13:45:11
   11.214916 Elapsed time: 0 minutes
   11.214923 Number of molecules: 895
   11.214928 Number of complexes: 895
   11.214934 Total reaction matches: 8533979
   11.214943 Avg timestep duration: 3.185850e-04, iterations remaining: 67000, Time left: 0 minutes
   11.214949 Estimated end time: 2026-10-14 13:45:32
   11.214958 -------------------------------------------------
   11.558813 -------------------
   11.558848 End iteration: 34000, simulation time: 3.400000e-02 seconds.
   11.558954 System time: 2026-10-14 13:45:11
   11.558966 Elapsed time: 0 minutes
   11.558972 Number of molecules: 895
   11.558978 Number of complexes: 895
   11.558983 Total reaction matches: 8799413
   11.558993 Avg timestep duration: 3.438040e-04, iterations remaining: 66000, Time left: 0 minutes
   11.558999 Estimated end time: 2026-10-14 13:45:34
   11.559036 -------------------------------------------------
   11.890839 -------------------
   11.890872 End iteration: 35000, simulation time: 3.500000e-02 seconds.
   11.890960 System time: 2026-10-14 13:45:11
   11.890970 Elapsed time: 0 minutes
   11.890975 Number of molecules: 895
   11.890980 Number of complexes: 895
   11.890984 Total reaction matches: 9013853
   11.890992 Avg timestep duration: 3.317810e-04, iterations remaining: 65000, Time left: 0 minutes
   11.890997 Estimated end time: 2026-10-14 13:45:33
   11.891001 -------------------------------------------------
   12.179497 -------------------
   12.179525 End iteration: 36000, simulation time: 3.600000e-02 seconds.
   12.179618 System time: 2026-10-14 13:45:12
   12.179628 Elapsed time: 0 minutes
   12.179633 Number of molecules: 895
   12.179637 Number of complexes: 895
   12.179641 Total reaction matches: 9087678
   12.179650 Avg timestep duration: 2.884440e-04, iterations remaining: 64000, Time left: 0 minutes
   12.179655 Estimated end time: 2026-10-14 13:45:30
   12.179659 -------------------------------------------------
   12.475370 -------------------
   12.475402 End iteration: 37000, simulation time: 3.700000e-02 seconds.
   12.475493 System time: 2026-10-14 13:45:12
   12.475503 Elapsed time: 0 minutes
   12.475509 Number of molecules: 895
   12.475514 Number of complexes: 895
   12.475519 Total reaction matches: 9204810
   12.475528 Avg timestep duration: 2.956330e-04, iterations remaining: 63000, Time left: 0 minutes
   12.475534 Estimated end time: 2026-10-14 13:45:31
   12.475538 -------------------------------------------------
   12.771409 -------------------
   12.771437 End iteration: 38000, simulation time: 3.800000e-02 seconds.
   12.771528 System time: 2026-10-14 13:45:12
   12.771537 Elapsed time: 0 minutes
   12.771542 Number of molecules: 895
   12.771547 Number of complexes: 895
   12.771551 Total reaction matches: 9321331
   12.771560 Avg timestep duration: 2.958190e-04, iterations remaining: 62000, Time left: 0 minutes
   12.771565 Estimated end time: 2026-10-14 13:45:31
   12.771569 -------------------------------------------------
   13.071909 -------------------
   13.071942 End iteration: 39000, simulation time: 3.900000e-02 seconds.
   13.072033 System time: 2026-10-14 13:45:13
   13.072042 Elapsed time: 0 minutes
   13.072050 Number of molecules: 895
   13.072056 Number of complexes: 895
   13.072061 Total reaction matches: 9459259
   13.072069 Avg timestep duration: 3.002830e-04, iterations remaining: 61000, Time left: 0 minutes
   13.072074 Estimated end time: 2026-10-14 13:45:31
   13.072078 -------------------------------------------------
   13.370160 -------------------
   13.370209 End iteration: 40000, simulation time: 4.000000e-02 seconds.
   13.370350 System time: 2026-10-14 13:45:13
   13.370360 Elapsed time: 0 minutes
   13.370366 Number of molecules: 895
   13.370371 Number of complexes: 895
   13.370376 Total reaction matches: 9542581
   13.370384 Avg timestep duration: 2.979670e-04, iterations remaining: 60000, Time left: 0 minutes
   13.370390 Estimated end time: 2026-10-14 13:45:31
   13.370395 -------------------------------------------------
   13.726504 -------------------
   13.726537 End iteration: 41000, simulation time: 4.100000e-02 seconds.
   13.726611 System time: 2026-10-14 13:45:13
   13.726622 Elapsed time: 0 minutes
   13.726628 Number of molecules: 895
   13.726633 Number of complexes: 895
   13.726637 Total reaction matches: 9861714
   13.726647 Avg timestep duration: 3.560440e-04, iterations remaining: 59000, Time left: 0 minutes
   13.726652 Estimated end time: 2026-10-14 13:45:34
   13.726656 -------------------------------------------------
   14.083580 -------------------
   14.083611 End iteration: 42000, simulation time: 4.200000e-02 seconds.
   14.083709 System time: 2026-10-14 13:45:14
   14.083719 Elapsed time: 0 minutes
   14.083724 Number of molecules: 895
   14.083729 Number of complexes: 895
   14.083733 Total reaction matches: 10189925
   14.083754 Avg timestep duration: 3.568740e-04, iterations remaining: 58000, Time left: 0 minutes
   14.083759 Estimated end time: 2026-10-14 13:45:34
   14.083764 -------------------------------------------------
   14.388937 -------------------
   14.388991 End iteration: 43000, simulation time: 4.300000e-02 seconds.
   14.389076 System time: 2026-10-14 13:45:14
   14.389087 Elapsed time: 0 minutes
   14.389093 Number of molecules: 895
   14.389098 Number of complexes: 895
   14.389103 Total reaction matches: 10353000
   14.389112 Avg timestep duration: 3.050740e-04, iterations remaining: 57000, Time left: 0 minutes
   14.389118 Estimated end time: 2026-10-14 13:45:31
   14.389122 -------------------------------------------------
   14.687800 -------------------
   14.687833 End iteration: 44000, simulation time: 4.400000e-02 seconds.
   14.687935 System time: 2026-10-14 13:45:14
   14.687944 Elapsed time: 0 minutes
   14.687950 Number of molecules: 895
   14.687955 Number of complexes: 895
   14.687959 Total reaction matches: 10459023
   14.687968 Avg timestep duration: 2.986270e-04, iterations remaining: 56000, Time left: 0 minutes
   14.687973 Estimated end time: 2026-10-14 13:45:31
   14.687977 -------------------------------------------------
   15.003949 -------------------
   15.003977 End iteration: 45000, simulation time: 4.500000e-02 seconds.
   15.004073 System time: 2026-10-14 13:45:15
   15.004082 Elapsed time: 0 minutes
   15.004087 Number of molecules: 895
   15.004092 Number of complexes: 895
   15.004096 Total reaction matches: 10651232
   15.004104 Avg timestep duration: 3.159260e-04, iterations remaining: 55000, Time left: 0 minutes
   15.004109 Estimated end time: 2026-10-14 13:45:32
   15.004115 -------------------------------------------------
   15.334675 -------------------
   15.334726 End iteration: 46000, simulation time: 4.600000e-02 seconds.
   15.334864 System time: 2026-10-14 13:45:15
   15.334876 Elapsed time: 0 minutes
   15.334882 Number of molecules: 895
   15.334887 Number of complexes: 895
   15.334891 Total reaction matches: 10859226
   15.334900 Avg timestep duration: 3.303980e-04, iterations remaining: 54000, Time left: 0 minutes
   15.334906 Estimated end time: 2026-10-14 13:45:33
   15.334910 -------------------------------------------------
   15.650587 -------------------
   15.650629 End iteration: 47000, simulation time: 4.700000e-02 seconds.
   15.650753 System time: 2026-10-14 13:45:15
   15.650763 Elapsed time: 0 minutes
   15.650768 Number of molecules: 895
   15.650773 Number of complexes: 895
   15.650778 Total reaction matches: 10959047
   15.650787 Avg timestep duration: 3.155860e-04, iterations remaining: 53000, Time left: 0 minutes
   15.650793 Estimated end time: 2026-10-14 13:45:32
   15.650798 -------------------------------------------------
   15.991608 -------------------
   15.991635 End iteration: 48000, simulation time: 4.800000e-02 seconds.
   15.991750 System time: 2026-10-14 13:45:16
   15.991759 Elapsed time: 0 minutes
   15.991764 Number of molecules: 895
   15.991769 Number of complexes: 895
   15.991773 Total reaction matches: 11189128
   15.991781 Avg timestep duration: 3.407570e-04, iterations remaining: 52000, Time left: 0 minutes
   15.991786 Estimated end time: 2026-10-14 13:45:33
   15.991790 -------------------------------------------------
   16.365999 -------------------
   16.366045 End iteration: 49000, simulation time: 4.900000e-02 seconds.
   16.366130 System time: 2026-10-14 13:45:16
   16.366141 Elapsed time: 0 minutes
   16.366148 Number of molecules: 895
   16.366153 Number of complexes: 895
   16.366158 Total reaction matches: 11508179
   16.366167 Avg timestep duration: 3.740880e-04, iterations remaining: 51000, Time left: 0 minutes
   16.366176 Estimated end time: 2026-10-14 13:45:35
   16.366181 -------------------------------------------------
   16.735627 -------------------
   16.735663 End iteration: 50000, simulation time: 5.000000e-02 seconds.
   16.735800 System time: 2026-10-14 13:45:16
   16.735810 Elapsed time: 0 minutes
   16.735823 Number of molecules: 895
   16.735829 Number of complexes: 895
   16.735833 Total reaction matches: 11770511
   16.735841 Avg timestep duration: 3.693590e-04, iterations remaining: 50000, Time left: 0 minutes
   16.735847 Estimated end time: 2026-10-14 13:45:35
   16.735851 -------------------------------------------------
   17.109188 -------------------
   17.109215 End iteration: 51000, simulation time: 5.100000e-02 seconds.
   17.109312 System time: 2026-10-14 13:45:17
   17.109322 Elapsed time: 0 minutes
   17.109327 Number of molecules: 895
   17.109332 Number of complexes: 895
   17.109336 Total reaction matches: 12086723
   17.109344 Avg timestep duration: 3.732960e-04, iterations remaining: 49000, Time left: 0 minutes
   17.109349 Estimated end time: 2026-10-14 13:45:35
   17.109356 -------------------------------------------------
   17.471753 -------------------
   17.471786 End iteration: 52000, simulation time: 5.200000e-02 seconds.
   17.471864 System time: 2026-10-14 13:45:17
   17.471874 Elapsed time: 0 minutes
   17.471879 Number of molecules: 895
   17.471883 Number of complexes: 895
   17.471887 Total reaction matches: 12397713
   17.471895 Avg timestep duration: 3.623060e-04, iterations remaining: 48000, Time left: 0 minutes
   17.471901 Estimated end time: 2026-10-14 13:45:34
   17.471905 -------------------------------------------------
   17.820880 -------------------
   17.820905 End iteration: 53000, simulation time: 5.300000e-02 seconds.
   17.821004 System time: 2026-10-14 13:45:17
   17.821013 Elapsed time: 0 minutes
   17.821019 Number of molecules: 895
   17.821023 Number of complexes: 895
   17.821027 Total reaction matches: 12684423
   17.821035 Avg timestep duration: 3.489330e-04, iterations remaining: 47000, Time left: 0 minutes
   17.821040 Estimated end time: 2026-10-14 13:45:34
   17.821044 -------------------------------------------------
   18.154855 -------------------
   18.154888 End iteration: 54000, simulation time: 5.400000e-02 seconds.
   18.154979 System time: 2026-10-14 13:45:18
   18.154988 Elapsed time: 0 minutes
   18.154994 Number of molecules: 895
   18.154999 Number of complexes: 895
   18.155004 Total reaction matches: 12903525
   18.155012 Avg timestep duration: 3.337410e-04, iterations remaining: 46000, Time left: 0 minutes
   18.155018 Estimated end time: 2026-10-14 13:45:33
   18.155023 -------------------------------------------------
   18.494177 -------------------
   18.494217 End iteration: 55000, simulation time: 5.500000e-02 seconds.
   18.494288 System time: 2026-10-14 13:45:18
   18.494298 Elapsed time: 0 minutes
   18.494304 Number of molecules: 895
   18.494309 Number of complexes: 895
   18.494313 Total reaction matches: 13146370
   18.494321 Avg timestep duration: 3.390490e-04, iterations remaining: 45000, Time left: 0 minutes
   18.494326 Estimated end time: 2026-10-14 13:45:33
   18.494331 -------------------------------------------------
   18.822177 -------------------
   18.822202 End iteration: 56000, simulation time: 5.600000e-02 seconds.
   18.822314 System time: 2026-10-14 13:45:18
   18.822324 Elapsed time: 0 minutes
   18.822329 Number of molecules: 895
   18.822339 Number of complexes: 895
   18.822344 Total reaction matches: 13357327
   18.822353 Avg timestep duration: 3.277930e-04, iterations remaining: 44000, Time left: 0 minutes
   18.822358 Estimated end time: 2026-10-14 13:45:33
   18.822363 -------------------------------------------------
   19.137363 -------------------
   19.137393 End iteration: 57000, simulation time: 5.700000e-02 seconds.
   19.137488 System time: 2026-10-14 13:45:19
   19.137498 Elapsed time: 0 minutes
   19.137504 Number of molecules: 895
   19.137508 Number of complexes: 895
   19.137513 Total reaction matches: 13533890
   19.137523 Avg timestep duration: 3.149500e-04, iterations remaining: 43000, Time left: 0 minutes
   19.137528 Estimated end time: 2026-10-14 13:45:32
   19.137532 -------------------------------------------------
   19.450558 -------------------
   19.450604 End iteration: 58000, simulation time: 5.800000e-02 seconds.
   19.450691 System time: 2026-10-14 13:45:19
   19.450704 Elapsed time: 0 minutes
   19.450710 Number of molecules: 895
   19.450715 Number of complexes: 895
   19.450720 Total reaction matches: 13691978
   19.450731 Avg timestep duration: 3.129310e-04, iterations remaining: 42000, Time left: 0 minutes
   19.450738 Estimated end time: 2026-10-14 13:45:32
   19.450743 -------------------------------------------------
   19.786094 -------------------
   19.786125 End iteration: 59000, simulation time: 5.900000e-02 seconds.
   19.786219 System time: 2026-10-14 13:45:19
   19.786228 Elapsed time: 0 minutes
   19.786234 Number of molecules: 895
   19.786239 Number of complexes: 895
   19.786244 Total reaction matches: 13903223
   19.786252 Avg timestep duration: 3.352870e-04, iterations remaining: 41000, Time left: 0 minutes
   19.786258 Estimated end time: 2026-10-14 13:45:33
   19.786263 -------------------------------------------------
   20.191248 -------------------
   20.191278 End iteration: 60000, simulation time: 6.000000e-02 seconds.
   20.191401 System time: 2026-10-14 13:45:20
   20.191412 Elapsed time: 0 minutes
   20.191417 Number of molecules: 895
   20.191422 Number of complexes: 895
   20.191427 Total reaction matches: 14270743
   20.191435 Avg timestep duration: 4.049290e-04, iterations remaining: 40000, Time left: 0 minutes
   20.191446 Estimated end time: 2026-10-14 13:45:36
   20.191450 -------------------------------------------------
   20.501821 -------------------
   20.501851 End iteration: 61000, simulation time: 6.100000e-02 seconds.
   20.501951 System time: 2026-10-14 13:45:20
   20.501960 Elapsed time: 0 minutes
   20.501966 Number of molecules: 895
   20.501970 Number of complexes: 895
   20.501975 Total reaction matches: 14408703
   20.501983 Avg timestep duration: 3.103390e-04, iterations remaining: 39000, Time left: 0 minutes
   20.501988 Estimated end time: 2026-10-14 13:45:32
   20.501993 -------------------------------------------------
   20.797357 -------------------
   20.797385 End iteration: 62000, simulation time: 6.200000e-02 seconds.
   20.797479 System time: 2026-10-14 13:45:20
   20.797489 Elapsed time: 0 minutes
   20.797495 Number of molecules: 895
   20.797500 Number of complexes: 895
   20.797505 Total reaction matches: 14553068
   20.797513 Avg timestep duration: 2.952950e-04, iterations remaining: 38000, Time left: 0 minutes
   20.797519 Estimated end time: 2026-10-14 13:45:32
   20.797523 -------------------------------------------------
   21.074810 -------------------
   21.074839 End iteration: 63000, simulation time: 6.300000e-02 seconds.
   21.074939 System time: 2026-10-14 13:45:21
   21.074949 Elapsed time: 0 minutes
   21.074954 Number of molecules: 895
   21.074958 Number of complexes: 895
   21.074962 Total reaction matches: 14648576
   21.074975 Avg timestep duration: 2.772330e-04, iterations remaining: 37000, Time left: 0 minutes
   21.074981 Estimated end time: 2026-10-14 13:45:31
   21.074986 -------------------------------------------------
   21.397666 -------------------
   21.397709 End iteration: 64000, simulation time: 6.400000e-02 seconds.
   21.397822 System time: 2026-10-14 13:45:21
   21.397833 Elapsed time: 0 minutes
   21.397838 Number of molecules: 895
   21.397843 Number of complexes: 895
   21.397847 Total reaction matches: 14883668
   21.397856 Avg timestep duration: 3.225980e-04, iterations remaining: 36000, Time left: 0 minutes
   21.397862 Estimated end time: 2026-10-14 13:45:33
   21.397867 -------------------------------------------------
   21.736777 -------------------
   21.736821 End iteration: 65000, simulation time: 6.500000e-02 seconds.
   21.736906 System time: 2026-10-14 13:45:21
   21.736917 Elapsed time: 0 minutes
   21.736924 Number of molecules: 895
   21.736929 Number of complexes: 895
   21.736934 Total reaction matches: 15168249
   21.736943 Avg timestep duration: 3.388150e-04, iterations remaining: 35000, Time left: 0 minutes
   21.736949 Estimated end time: 2026-10-14 13:45:33
   21.736965 -------------------------------------------------
   22.063995 -------------------
   22.064024 End iteration: 66000, simulation time: 6.600000e-02 seconds.
   22.064117 System time: 2026-10-14 13:45:22
   22.064127 Elapsed time: 0 minutes
   22.064132 Number of molecules: 895
   22.064137 Number of complexes: 895
   22.064141 Total reaction matches: 15387965
   22.064150 Avg timestep duration: 3.269810e-04, iterations remaining: 34000, Time left: 0 minutes
   22.064155 Estimated end time: 2026-10-14 13:45:33
   22.064159 -------------------------------------------------
   22.403306 -------------------
   22.403334 End iteration: 67000, simulation time: 6.700000e-02 seconds.
   22.403436 System time: 2026-10-14 13:45:22
   22.403445 Elapsed time: 0 minutes
   22.403450 Number of molecules: 895
   22.403455 Number of complexes: 895
   22.403459 Total reaction matches: 15666452
   22.403468 Avg timestep duration: 3.391000e-04, iterations remaining: 33000, Time left: 0 minutes
   22.403473 Estimated end time: 2026-10-14 13:45:33
   22.403478 -------------------------------------------------
   22.760521 -------------------
   22.760557 End iteration: 68000, simulation time: 6.800000e-02 seconds.
   22.760643 System time: 2026-10-14 13:45:22
   22.760653 Elapsed time: 0 minutes
   22.760659 Number of molecules: 895
   22.760664 Number of complexes: 895
   22.760669 Total reaction matches: 15993355
   22.760679 Avg timestep duration: 3.569640e-04, iterations remaining: 32000, Time left: 0 minutes
   22.760685 Estimated end time: 2026-10-14 13:45:34
   22.760690 -------------------------------------------------
   23.119238 -------------------
   23.119270 End iteration: 69000, simulation time: 6.900000e-02 seconds.
   23.119371 System time: 2026-10-14 13:45:23
   23.119381 Elapsed time: 0 minutes
   23.119389 Number of molecules: 895
   23.119394 Number of complexes: 895
   23.119398 Total reaction matches: 16336209
   23.119406 Avg timestep duration: 3.584880e-04, iterations remaining: 31000, Time left: 0 minutes
   23.119412 Estimated end time: 2026-10-14 13:45:34
   23.119416 -------------------------------------------------
   23.473695 -------------------
   23.473727 End iteration: 70000, simulation time: 7.000000e-02 seconds.
   23.473850 System time: 2026-10-14 13:45:23
   23.473865 Elapsed time: 0 minutes
   23.473870 Number of molecules: 895
   23.473875 Number of complexes: 895
   23.473879 Total reaction matches: 16594924
   23.473887 Avg timestep duration: 3.542190e-04, iterations remaining: 30000, Time left: 0 minutes
   23.473893 Estimated end time: 2026-10-14 13:45:34
   23.473898 -------------------------------------------------
   23.817562 -------------------
   23.817604 End iteration: 71000, simulation time: 7.100000e-02 seconds.
   23.817725 System time: 2026-10-14 13:45:23
   23.817736 Elapsed time: 0 minutes
   23.817742 Number of molecules: 895
   23.817748 Number of complexes: 895
   23.817753 Total reaction matches: 16861581
   23.817762 Avg timestep duration: 3.435810e-04, iterations remaining: 29000, Time left: 0 minutes
   23.817769 Estimated end time: 2026-10-14 13:45:33
   23.817774 -------------------------------------------------
   24.204216 -------------------
   24.204249 End iteration: 72000, simulation time: 7.200000e-02 seconds.
   24.204339 System time: 2026-10-14 13:45:24
   24.204350 Elapsed time: 0 minutes
   24.204355 Number of molecules: 895
   24.204360 Number of complexes: 895
   24.204365 Total reaction matches: 17246280
   24.204374 Avg timestep duration: 3.863840e-04, iterations remaining: 28000, Time left: 0 minutes
   24.204379 Estimated end time: 2026-10-14 13:45:35
   24.204384 -------------------------------------------------
   24.575711 -------------------
   24.575740 End iteration: 73000, simulation time: 7.300000e-02 seconds.
   24.575842 System time: 2026-10-14 13:45:24
   24.575852 Elapsed time: 0 minutes
   24.575857 Number of molecules: 895
   24.575862 Number of complexes: 895
   24.575866 Total reaction matches: 17619589
   24.575874 Avg timestep duration: 3.712750e-04, iterations remaining: 27000, Time left: 0 minutes
   24.575891 Estimated end time: 2026-10-14 13:45:34
   24.575896 -------------------------------------------------
   24.891329 -------------------
   24.891360 End iteration: 74000, simulation time: 7.400000e-02 seconds.
   24.891467 System time: 2026-10-14 13:45:24
   24.891477 Elapsed time: 0 minutes
   24.891483 Number of molecules: 895
   24.891488 Number of complexes: 895
   24.891493 Total reaction matches: 17782739
   24.891502 Avg timestep duration: 3.153800e-04, iterations remaining: 26000, Time left: 0 minutes
   24.891507 Estimated end time: 2026-10-14 13:45:33
   24.891511 -------------------------------------------------
   25.185450 -------------------
   25.185480 End iteration: 75000, simulation time: 7.500000e-02 seconds.
   25.185592 System time: 2026-10-14 13:45:25
   25.185601 Elapsed time: 0 minutes
   25.185607 Number of molecules: 895
   25.185611 Number of complexes: 895
   25.185615 Total reaction matches: 17930097
   25.185624 Avg timestep duration: 2.938850e-04, iterations remaining: 25000, Time left: 0 minutes
   25.185629 Estimated end time: 2026-10-14 13:45:32
   25.185633 -------------------------------------------------
   25.495025 -------------------
   25.495075 End iteration: 76000, simulation time: 7.600000e-02 seconds.
   25.495154 System time: 2026-10-14 13:45:25
   25.495164 Elapsed time: 0 minutes
   25.495170 Number of molecules: 895
   25.495175 Number of complexes: 895
   25.495180 Total reaction matches: 18107507
   25.495190 Avg timestep duration: 3.092840e-04, iterations remaining: 24000, Time left: 0 minutes
   25.495197 Estimated end time: 2026-10-14 13:45:33
   25.495205 -------------------------------------------------
   25.844466 -------------------
   25.844505 End iteration: 77000, simulation time: 7.700000e-02 seconds.
   25.844591 System time: 2026-10-14 13:45:25
   25.844602 Elapsed time: 0 minutes
   25.844609 Number of molecules: 895
   25.844614 Number of complexes: 895
   25.844619 Total reaction matches: 18359935
   25.844629 Avg timestep duration: 3.491910e-04, iterations remaining: 23000, Time left: 0 minutes
   25.844635 Estimated end time: 2026-10-14 13:45:33
   25.844640 -------------------------------------------------
   26.193303 -------------------
   26.193338 End iteration: 78000, simulation time: 7.800000e-02 seconds.
   26.193521 System time: 2026-10-14 13:45:26
   26.193544 Elapsed time: 0 minutes
   26.193553 Number of molecules: 895
   26.193559 Number of complexes: 895
   26.193571 Total reaction matches: 18668940
   26.193578 Avg timestep duration: 3.485920e-04, iterations remaining: 22000, Time left: 0 minutes
   26.193583 Estimated end time: 2026-10-14 13:45:33
   26.193588 -------------------------------------------------
   26.504502 -------------------
   26.504544 End iteration: 79000, simulation time: 7.900000e-02 seconds.
   26.504686 System time: 2026-10-14 13:45:26
   26.504697 Elapsed time: 0 minutes
   26.504703 Number of molecules: 895
   26.504708 Number of complexes: 895
   26.504713 Total reaction matches: 18827147
   26.504724 Avg timestep duration: 3.108340e-04, iterations remaining: 21000, Time left: 0 minutes
   26.504731 Estimated end time: 2026-10-14 13:45:33
   26.504736 -------------------------------------------------
   26.876165 -------------------
   26.876196 End iteration: 80000, simulation time: 8.000000e-02 seconds.
   26.876326 System time: 2026-10-14 13:45:26
   26.876345 Elapsed time: 0 minutes
   26.876351 Number of molecules: 895
   26.876356 Number of complexes: 895
   26.876360 Total reaction matches: 19132556
   26.876371 Avg timestep duration: 3.713210e-04, iterations remaining: 20000, Time left: 0 minutes
   26.876376 Estimated end time: 2026-10-14 13:45:34
   26.876380 -------------------------------------------------
   27.235515 -------------------
   27.235547 End iteration: 81000, simulation time: 8.100000e-02 seconds.
   27.235647 System time: 2026-10-14 13:45:27
   27.235657 Elapsed time: 0 minutes
   27.235662 Number of molecules: 895
   27.235677 Number of complexes: 895
   27.235681 Total reaction matches: 19431150
   27.235690 Avg timestep duration: 3.591130e-04, iterations remaining: 19000, Time left: 0 minutes
   27.235695 Estimated end time: 2026-10-14 13:45:34
   27.235700 -------------------------------------------------
   27.627719 -------------------
   27.627753 End iteration: 82000, simulation time: 8.200000e-02 seconds.
   27.627841 System time: 2026-10-14 13:45:27
   27.627851 Elapsed time: 0 minutes
   27.627856 Number of molecules: 895
   27.627861 Number of complexes: 895
   27.627865 Total reaction matches: 19862055
   27.627874 Avg timestep duration: 3.919600e-04, iterations remaining: 18000, Time left: 0 minutes
   27.627880 Estimated end time: 2026-10-14 13:45:34
   27.627884 -------------------------------------------------
   27.987789 -------------------
   27.987817 End iteration: 83000, simulation time: 8.300000e-02 seconds.
   27.987911 System time: 2026-10-14 13:45:28
   27.987920 Elapsed time: 0 minutes
   27.987926 Number of molecules: 895
   27.987931 Number of complexes: 895
   27.987936 Total reaction matches: 20182476
   27.987945 Avg timestep duration: 3.598630e-04, iterations remaining: 17000, Time left: 0 minutes
   27.987950 Estimated end time: 2026-10-14 13:45:34
   27.987965 -------------------------------------------------
   28.329702 -------------------
   28.329749 End iteration: 84000, simulation time: 8.400000e-02 seconds.
   28.329839 System time: 2026-10-14 13:45:28
   28.329849 Elapsed time: 0 minutes
   28.329855 Number of molecules: 895
   28.329860 Number of complexes: 895
   28.329865 Total reaction matches: 20499425
   28.329876 Avg timestep duration: 3.416240e-04, iterations remaining: 16000, Time left: 0 minutes
   28.329882 Estimated end time: 2026-10-14 13:45:33
   28.329887 -------------------------------------------------
   28.680489 -------------------
   28.680522 End iteration: 85000, simulation time: 8.500000e-02 seconds.
   28.680627 System time: 2026-10-14 13:45:28
   28.680637 Elapsed time: 0 minutes
   28.680643 Number of molecules: 895
   28.680649 Number of complexes: 895
   28.680653 Total reaction matches: 20779145
   28.680662 Avg timestep duration: 3.505400e-04, iterations remaining: 15000, Time left: 0 minutes
   28.680668 Estimated end time: 2026-10-14 13:45:34
   28.680672 -------------------------------------------------
   29.030308 -------------------
   29.030334 End iteration: 86000, simulation time: 8.600000e-02 seconds.
   29.030432 System time: 2026-10-14 13:45:29
   29.030441 Elapsed time: 0 minutes
   29.030446 Number of molecules: 895
   29.030451 Number of complexes: 895
   29.030455 Total reaction matches: 21062094
   29.030464 Avg timestep duration: 3.495860e-04, iterations remaining: 14000, Time left: 0 minutes
   29.030469 Estimated end time: 2026-10-14 13:45:34
   29.030473 -------------------------------------------------
   29.351115 -------------------
   29.351187 End iteration: 87000, simulation time: 8.700000e-02 seconds.
   29.351263 System time: 2026-10-14 13:45:29
   29.351275 Elapsed time: 0 minutes
   29.351281 Number of molecules: 895
   29.351285 Number of complexes: 895
   29.351291 Total reaction matches: 21268290
   29.351300 Avg timestep duration: 3.204870e-04, iterations remaining: 13000, Time left: 0 minutes
   29.351306 Estimated end time: 2026-10-14 13:45:33
   29.351311 -------------------------------------------------
   29.634544 -------------------
   29.634571 End iteration: 88000, simulation time: 8.800000e-02 seconds.
   29.634656 System time: 2026-10-14 13:45:29
   29.634665 Elapsed time: 0 minutes
   29.634670 Number of molecules: 895
   29.634675 Number of complexes: 895
   29.634679 Total reaction matches: 21360437
   29.634688 Avg timestep duration: 2.831970e-04, iterations remaining: 12000, Time left: 0 minutes
   29.634694 Estimated end time: 2026-10-14 13:45:33
   29.634698 -------------------------------------------------
   29.918407 -------------------
   29.918437 End iteration: 89000, simulation time: 8.900000e-02 seconds.
   29.918508 System time: 2026-10-14 13:45:30
   29.918518 Elapsed time: 0 minutes
   29.918523 Number of molecules: 895
   29.918528 Number of complexes: 895
   29.918533 Total reaction matches: 21470201
   29.918541 Avg timestep duration: 2.836410e-04, iterations remaining: 11000, Time left: 0 minutes
   29.918546 Estimated end time: 2026-10-14 13:45:33
   29.918550 -------------------------------------------------
   30.244895 -------------------
   30.244923 End iteration: 90000, simulation time: 9.000000e-02 seconds.
   30.245049 System time: 2026-10-14 13:45:30
   30.245058 Elapsed time: 0 minutes
   30.245064 Number of molecules: 895
   30.245069 Number of complexes: 895
   30.245073 Total reaction matches: 21640833
   30.245082 Avg timestep duration: 3.262970e-04, iterations remaining: 10000, Time left: 0 minutes
   30.245087 Estimated end time: 2026-10-14 13:45:33
   30.245091 -------------------------------------------------
   30.601445 -------------------
   30.601479 End iteration: 91000, simulation time: 9.100000e-02 seconds.
   30.601604 System time: 2026-10-14 13:45:30
   30.601615 Elapsed time: 0 minutes
   30.601620 Number of molecules: 895
   30.601624 Number of complexes: 895
   30.601629 Total reaction matches: 21962989
   30.601641 Avg timestep duration: 3.562880e-04, iterations remaining: 9000, Time left: 0 minutes
   30.601646 Estimated end time: 2026-10-14 13:45:33
   30.601651 -------------------------------------------------
   30.913786 -------------------
   30.913820 End iteration: 92000, simulation time: 9.200000e-02 seconds.
   30.913906 System time: 2026-10-14 13:45:31
   30.913917 Elapsed time: 0 minutes
   30.913922 Number of molecules: 895
   30.913927 Number of complexes: 895
   30.913931 Total reaction matches: 22143942
   30.913947 Avg timestep duration: 3.121050e-04, iterations remaining: 8000, Time left: 0 minutes
   30.913952 Estimated end time: 2026-10-14 13:45:33
   30.913956 -------------------------------------------------
   31.221165 -------------------
   31.221200 End iteration: 93000, simulation time: 9.300000e-02 seconds.
   31.221299 System time: 2026-10-14 13:45:31
   31.221320 Elapsed time: 0 minutes
   31.221328 Number of molecules: 895
   31.221334 Number of complexes: 895
   31.221340 Total reaction matches: 22308880
   31.221345 Avg timestep duration: 3.071560e-04, iterations remaining: 7000, Time left: 0 minutes
   31.221357 Estimated end time: 2026-10-14 13:45:33
   31.221363 -------------------------------------------------
   31.559708 -------------------
   31.559741 End iteration: 94000, simulation time: 9.400000e-02 seconds.
   31.559831 System time: 2026-10-14 13:45:31
   31.559841 Elapsed time: 0 minutes
   31.559847 Number of molecules: 895
   31.559851 Number of complexes: 895
   31.559856 Total reaction matches: 22566652
   31.559865 Avg timestep duration: 3.383190e-04, iterations remaining: 6000, Time left: 0 minutes
   31.559870 Estimated end time: 2026-10-14 13:45:33
   31.559879 -------------------------------------------------
   31.852122 -------------------
   31.852147 End iteration: 95000, simulation time: 9.500000e-02 seconds.
   31.852237 System time: 2026-10-14 13:45:31
   31.852247 Elapsed time: 0 minutes
   31.852252 Number of molecules: 895
   31.852256 Number of complexes: 895
   31.852260 Total reaction matches: 22698126
   31.852269 Avg timestep duration: 2.921960e-04, iterations remaining: 5000, Time left: 0 minutes
   31.852274 Estimated end time: 2026-10-14 13:45:33
   31.852279 -------------------------------------------------
   32.156962 -------------------
   32.156991 End iteration: 96000, simulation time: 9.600000e-02 seconds.
   32.157090 System time: 2026-10-14 13:45:32
   32.157100 Elapsed time: 0 minutes
   32.157105 Number of molecules: 895
   32.157110 Number of complexes: 895
   32.157114 Total reaction matches: 22846686
   32.157123 Avg timestep duration: 3.046370e-04, iterations remaining: 4000, Time left: 0 minutes
   32.157127 Estimated end time: 2026-10-14 13:45:33
   32.157132 -------------------------------------------------
   32.534205 -------------------
   32.534235 End iteration: 97000, simulation time: 9.700000e-02 seconds.
   32.534414 System time: 2026-10-14 13:45:32
   32.534459 Elapsed time: 0 minutes
   32.534488 Number of molecules: 895
   32.534501 Number of complexes: 895
   32.534507 Total reaction matches: 23197245
   32.534511 Avg timestep duration: 3.770040e-04, iterations remaining: 3000, Time left: 0 minutes
   32.534517 Estimated end time: 2026-10-14 13:45:33
   32.534521 -------------------------------------------------
   32.843488 -------------------
   32.843515 End iteration: 98000, simulation time: 9.800000e-02 seconds.
   32.843596 System time: 2026-10-14 13:45:32
   32.843605 Elapsed time: 0 minutes
   32.843609 Number of molecules: 895
   32.843614 Number of complexes: 895
   32.843618 Total reaction matches: 23349669
   32.843626 Avg timestep duration: 3.090670e-04, iterations remaining: 2000, Time left: 0 minutes
   32.843631 Estimated end time: 2026-10-14 13:45:33
   32.843635 -------------------------------------------------
   33.159210 -------------------
   33.159260 End iteration: 99000, simulation time: 9.900000e-02 seconds.
   33.159462 System time: 2026-10-14 13:45:33
   33.159478 Elapsed time: 0 minutes
   33.159484 Number of molecules: 895
   33.159490 Number of complexes: 895
   33.159495 Total reaction matches: 23532933
   33.159505 Avg timestep duration: 3.155160e-04, iterations remaining: 1000, Time left: 0 minutes
   33.159511 Estimated end time: 2026-10-14 13:45:33
   33.159515 -------------------------------------------------
   33.502090 -------------------------------------------------
   33.502129 End simulation
   33.502137 End date: 2026-10-14 13:45:33
   33.502142 	Wall Time: 3.349526e+01 seconds
   33.503315 [exit status 0]
