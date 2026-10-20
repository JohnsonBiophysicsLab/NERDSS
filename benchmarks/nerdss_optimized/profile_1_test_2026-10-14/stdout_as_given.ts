# t0_monotonic 840663.028543
    0.008438 Command: /private/tmp/claude-501/-Users-yue-Workspace-NERDSS/795ae25f-3942-4661-85be-f52029fde836/scratchpad/nerdss_opt_g -f parms_0.inp
    0.008768 Start date: 2026-10-14 13:36:08
    0.008785 RNG Seed: 2190169512
    0.009015 
    0.009028 Parsing Input:
    0.009040 This is a new simulation with input file: parms_0.inp
    0.009484 Parsing simulation parameters...
    0.009509 Read in nItr: 100000 timeSteps
    0.009526 Read in timeStep: 1 us
    0.009541 Read in timeWrite: 1000 timeSteps
    0.009556 Read in trajWrite: -1 timeSteps
    0.009569 Read in restartWrite: 2000000 timeSteps
    0.009588 Read in checkPoint: -1 timeSteps
    0.009604 Read in pdbWrite: -1 timeSteps
    0.009618 Parsing simulation boundary conditions...
    0.009631 Read in waterBox: [125.387 nm, 125.387 nm, 5470.63 nm]
    0.009647 Parsing mol file information...
    0.009660 lipid.mol
    0.009755 Read in isLipid: true
    0.009767 Read in checkOverlap: false
    0.009888 Read in D: [1um^2s^-1, 1um^2s^-1, 0um^2s^-1]
    0.009924 Read in Dr: [0rad^2s^-1, 0rad^2s^-1, 0rad^2s^-1]
    0.009945 Coordinates:
    0.009965 COM: [0nm, 0nm, 0nm]
    0.010005 head: [0nm, 0nm, 1nm]
    0.010041 Bonds:
    0.010053 com-head,
    0.010063 States:
    0.010074 head has 2 state(s): U	P
    0.010084 radius calculated from the ifacesCoord: 1 nm
    0.010095 mass auto calculated from the radius: 1
    0.010107 
    0.010119 adaptor.mol
    0.010197 Read in checkOverlap: false
    0.010220 Read in D: [25um^2s^-1, 25um^2s^-1, 25um^2s^-1]
    0.010238 Read in Dr: [0.5rad^2s^-1, 0.5rad^2s^-1, 0.5rad^2s^-1]
    0.010262 Coordinates:
    0.010281 COM: [0nm, 0nm, 0nm]
    0.010297 lbs: [0nm, 0nm, -2nm]
    0.010313 ebs: [0nm, 1nm, -1nm]
    0.010344 Bonds:
    0.010361 com-lbs, com-ebs, Invalid atom in adaptor mol file. Ignoring bonds...
    0.010377 radius calculated from the ifacesCoord: 2 nm
    0.010390 mass auto calculated from the radius: 2
    0.010412 
    0.010440 kinase.mol
    0.010459 Read in checkOverlap: false
    0.010473 Read in D: [25um^2s^-1, 25um^2s^-1, 25um^2s^-1]
    0.010498 Read in Dr: [0.5rad^2s^-1, 0.5rad^2s^-1, 0.5rad^2s^-1]
    0.010514 Coordinates:
    0.010529 COM: [0nm, 0nm, 0nm]
    0.010548 abs: [0nm, 0nm, 0.5nm]
    0.010563 lbs: [0nm, 0nm, -0.5nm]
    0.010577 Bonds:
    0.010591 com-abs, com-lbs,
    0.010604 radius calculated from the ifacesCoord: 0.5 nm
    0.010618 mass auto calculated from the radius: 0.5
    0.010635 
    0.010649 phosphatase.mol
    0.010679 Read in checkOverlap: false
    0.010700 Read in D: [25um^2s^-1, 25um^2s^-1, 25um^2s^-1]
    0.010710 Read in Dr: [0.5rad^2s^-1, 0.5rad^2s^-1, 0.5rad^2s^-1]
    0.010720 Coordinates:
    0.010730 COM: [0nm, 0nm, 0nm]
    0.010740 abs: [0nm, 0nm, 0.5nm]
    0.010750 lbs: [0nm, 0nm, -0.5nm]
    0.010777 Bonds:
    0.010791 com-abs, com-lbs,
    0.010808 radius calculated from the ifacesCoord: 0.5 nm
    0.010821 mass auto calculated from the radius: 0.5
    0.010830 
    0.010839 Parsing reactions...
    0.010850 Parsing new reaction with form:
    0.010862 	kinase(lbs) + lipid(head~U) <-> kinase(lbs!1).lipid(head~U!1)
    0.010889 Reversible reaction.
    0.010902 Bimolecular reaction detected.
    0.010925 Parsing reactants...
    0.010937 Parsing products...
    0.010961 Read in value of onRate3DMacro: 0.00137276uM^-1s^-1
    0.010972 Read in value of offRateMacro: 3.31104s^-1
    0.010981 Read in value of kcat: 0.0894988
    0.010990 Read in value of norm1: [0i + 0j + 1k]
    0.011022 Read in value of norm2: [0i + 0j + 1k]
    0.011058 Read in value of assocAngles: [theta1: 1.5708, theta2: 1.5708, phi1: nan, phi2: nan, omega: 3.14159]
    0.011071 Read in value of sigma: 1nm
    0.011082 Read in value of length3Dto2D: 1nm
    0.011092 Read in value of coupledRxn: label, PhosphorylationByKinase
    0.011112 Finished reading reaction parameters, constructing reaction.
    0.011158 Checking for valid states for reactant kinase
    0.011168 Checking for valid states for reactant lipid
    0.011178 Found state, looking for state change...
    0.011201 Found corresponding interfaces: React: kinase(lbs) -- Prod: kinase(lbs)
    0.011218 Found corresponding interfaces: React: lipid(head) -- Prod: lipid(head)
    0.011247 Allow Michaelis-Menten for this reaction. Probability of catalysis: kcat/(kcat+offRatekb): 0.0263186 total rate: 3.40059
    0.011272 Parameter set complete for reaction [kinase(lbs)+lipid(head~U)<->kinase(lbs!1).lipid(head~U!1)].
    0.011282 
    0.011289 Parsing new reaction with form:
    0.011298 	lipid(head~U) -> lipid(head~P)
    0.011307 Irreversible reaction.
    0.011312 UniMolStateChange reaction detected
    0.011318 Parsing reactants...
    0.011323 Parsing products...
    0.011354 Read in value of rate: 0
    0.011367 Read in value of rxnLabel: PhosphorylationByKinase
    0.011376 Finished reading reaction parameters, constructing reaction.
    0.011384 Checking for valid states for reactant lipid
    0.011391 Found state, looking for state change...
    0.011401 State change found on species 0. Interface head changes state from U to P.
    0.011410 Parameter set complete for reaction [lipid(head~U)->lipid(head~P)].
    0.011417 
    0.011426 Parsing new reaction with form:
    0.011434 	phosphatase(lbs) + lipid(head~P) <-> phosphatase(lbs!1).lipid(head~P!1)
    0.011443 Reversible reaction.
    0.011452 Bimolecular reaction detected.
    0.011459 Parsing reactants...
    0.011466 Parsing products...
    0.011473 Read in value of onRate3DMacro: 0.0691851uM^-1s^-1
    0.011479 Read in value of offRateMacro: 0.00443464s^-1
    0.011485 Read in value of kcat: 0.0276077
    0.011492 Read in value of norm1: [0i + 0j + 1k]
    0.011516 Read in value of norm2: [0i + 0j + 1k]
    0.011527 Read in value of assocAngles: [theta1: 1.5708, theta2: 1.5708, phi1: nan, phi2: nan, omega: 3.14159]
    0.011539 Read in value of sigma: 1nm
    0.011550 Read in value of length3Dto2D: 1nm
    0.011559 Read in value of coupledRxn: label, DephosphorylationByPhosphatase
    0.011568 Finished reading reaction parameters, constructing reaction.
    0.011577 Checking for valid states for reactant phosphatase
    0.011586 Checking for valid states for reactant lipid
    0.011595 Found state, looking for state change...
    0.011604 Found corresponding interfaces: React: phosphatase(lbs) -- Prod: phosphatase(lbs)
    0.011615 Found corresponding interfaces: React: lipid(head) -- Prod: lipid(head)
    0.011625 Allow Michaelis-Menten for this reaction. Probability of catalysis: kcat/(kcat+offRatekb): 0.861516 total rate: 0.0320455
    0.011634 Parameter set complete for reaction [phosphatase(lbs)+lipid(head~P)<->phosphatase(lbs!1).lipid(head~P!1)].
    0.011643 
    0.011653 Parsing new reaction with form:
    0.011661 	lipid(head~P) -> lipid(head~U)
    0.011671 Irreversible reaction.
    0.011682 UniMolStateChange reaction detected
    0.011692 Parsing reactants...
    0.011707 Parsing products...
    0.011717 Read in value of rate: 0
    0.011727 Read in value of rxnLabel: DephosphorylationByPhosphatase
    0.011738 Finished reading reaction parameters, constructing reaction.
    0.011748 Checking for valid states for reactant lipid
    0.011757 Found state, looking for state change...
    0.011765 State change found on species 0. Interface head changes state from P to U.
    0.011775 Parameter set complete for reaction [lipid(head~P)->lipid(head~U)].
    0.011786 
    0.011795 Parsing new reaction with form:
    0.011804 	lipid(head~P) + adaptor(lbs) <-> lipid(head~P!1).adaptor(lbs!1)
    0.011813 Reversible reaction.
    0.011822 Bimolecular reaction detected.
    0.011831 Parsing reactants...
    0.011843 Parsing products...
    0.011853 Read in value of onRate3DMacro: 1.16403uM^-1s^-1
    0.011862 Read in value of offRateMacro: 0.0865906s^-1
    0.011874 Read in value of norm1: [0i + 0j + 1k]
    0.011884 Read in value of norm2: [0i + 0j + 1k]
    0.011893 Read in value of assocAngles: [theta1: 1.5708, theta2: 1.5708, phi1: nan, phi2: nan, omega: 3.14159]
    0.011903 Read in value of length3Dto2D: 1nm
    0.011912 Read in value of sigma: 1nm
    0.011920 Finished reading reaction parameters, constructing reaction.
    0.011953 Checking for valid states for reactant lipid
    0.011963 Found state, looking for state change...
    0.011974 Checking for valid states for reactant adaptor
    0.011984 Found corresponding interfaces: React: lipid(head) -- Prod: lipid(head)
    0.011993 Found corresponding interfaces: React: adaptor(lbs) -- Prod: adaptor(lbs)
    0.012002 Parameter set complete for reaction [lipid(head~P)+adaptor(lbs)<->lipid(head~P!1).adaptor(lbs!1)].
    0.012011 
    0.012020 Parsing new reaction with form:
    0.012029 	kinase(abs) + adaptor(ebs,lbs!*) <-> kinase(abs!1).adaptor(ebs!1,lbs!*)
    0.012038 Reversible reaction.
    0.012047 Bimolecular reaction detected.
    0.012056 Parsing reactants...
    0.012064 Parsing products...
    0.012074 Read in value of onRate3DMacro: 5.49406uM^-1s^-1
    0.012083 Read in value of offRateMacro: 0.00612925s^-1
    0.012094 Read in value of norm1: [0i + 0j + 1k]
    0.012104 Read in value of norm2: [0i + 0j + 1k]
    0.012114 Read in value of assocAngles: [theta1: 1.5708, theta2: 2.35619, phi1: nan, phi2: 3.14159, omega: 3.14159]
    0.012123 Read in value of length3Dto2D: 1nm
    0.012133 Read in value of sigma: 1nm
    0.012143 Finished reading reaction parameters, constructing reaction.
    0.012154 Checking for valid states for reactant kinase
    0.012167 Checking for valid states for reactant adaptor
    0.012177 Found corresponding interfaces: React: kinase(abs) -- Prod: kinase(abs)
    0.012186 Found corresponding interfaces: React: adaptor(ebs) -- Prod: adaptor(ebs)
    0.012195 Parameter set complete for reaction [kinase(abs)+adaptor(ebs,lbs!*)<->kinase(abs!1).adaptor(ebs!1,lbs!*)].
    0.012204 
    0.012215 Parsing new reaction with form:
    0.012224 	phosphatase(abs) + adaptor(ebs,lbs!*) <-> phosphatase(abs!1).adaptor(ebs!1,lbs!*)
    0.012236 Reversible reaction.
    0.012243 Bimolecular reaction detected.
    0.012252 Parsing reactants...
    0.012261 Parsing products...
    0.012272 Read in value of onRate3DMacro: 0.153275uM^-1s^-1
    0.012281 Read in value of offRateMacro: 0.00108378s^-1
    0.012291 Read in value of norm1: [0i + 0j + 1k]
    0.012302 Read in value of norm2: [0i + 0j + 1k]
    0.012311 Read in value of assocAngles: [theta1: 1.5708, theta2: 2.35619, phi1: nan, phi2: 3.14159, omega: 3.14159]
    0.012320 Read in value of length3Dto2D: 1nm
    0.012330 Read in value of sigma: 1nm
    0.012340 Checking for valid states for reactant phosphatase
    0.012350 Checking for valid states for reactant adaptor
    0.012366 Found corresponding interfaces: React: phosphatase(abs) -- Prod: phosphatase(abs)
    0.012378 Found corresponding interfaces: React: adaptor(ebs) -- Prod: adaptor(ebs)
    0.012387 Parameter set complete for reaction [phosphatase(abs)+adaptor(ebs,lbs!*)<->phosphatase(abs!1).adaptor(ebs!1,lbs!*)].
    0.012397 
    0.012409 
    0.012420 Input file parsing complete
    0.012430 NUMBER OF MOLECULE TYPES: 4NUMBER OF INTERFACES PLUS STATES, including PRODUCTS: 13
    0.012814 NUMBER OF MOLECULES IN GEN COORDS: 895
    0.012824 Finding and fixing overlapping proteins.
    0.014474 Iteration: 1
    0.014484 Overlapping Proteins (false): 0
    0.014500 No overlapping proteins found.
    0.033469 
    0.033486 Partitioning simulation box into sub-boxes...
    0.033494 Rmaxlimit: 55.2036
    0.033563 Simulation volume parameters:
    0.033573 Total sub-volumes: 396
    0.033581 	Dimensions: [2, 2, 99]
    0.033589 	Maximum sub-volume neighbors: 13
    0.033598 	Sub-volume size: [62.6935, 62.6935, 55.2589]
    0.042133  Total species calculated from init_speciesFile.cpp: 13
    0.042157 Initialized copyNumSpecies. Printing all values. total species: 13
    0.042167  absolute index: 0 Copy nums: 278
    0.042176  absolute index: 1 Copy nums: 0
    0.042185  absolute index: 2 Copy nums: 100
    0.042195  absolute index: 3 Copy nums: 100
    0.042204  absolute index: 4 Copy nums: 296
    0.042212  absolute index: 5 Copy nums: 296
    0.042220  absolute index: 6 Copy nums: 221
    0.042238  absolute index: 7 Copy nums: 221
    0.042266  absolute index: 8 Copy nums: 0
    0.042275  absolute index: 9 Copy nums: 0
    0.042283  absolute index: 10 Copy nums: 0
    0.042292  absolute index: 11 Copy nums: 0
    0.042302  absolute index: 12 Copy nums: 0
    0.042602 Pro pair:  p1: 0 p2: 1 index: 1 lipid,adaptor'	Pro pair:  p1: 0 p2: 2 index: 2 lipid,kinase'	Pro pair:  p1: 0 p2: 3 index: 3 lipid,phosphatase'	Pro pair:  p1: 1 p2: 2 index: 6 adaptor,kinase'	Pro pair:  p1: 1 p2: 3 index: 7 adaptor,phosphatase'	size of proPairlist: 5 Number of pairs: 16 numMolTypes: 4
    0.043028 Spherical shell index: inactive.
    0.043044 
    0.043055 Simulation Parameters
    0.043063 Number of iterations: 100000 timesteps
    0.043073 Timestep: 1 us
    0.043084 Timestep log interval (timeWrite): 1000 timesteps
    0.043102 Restart file write interval: 2000000 timesteps
    0.043114 Coordinate write interval (trajWrite): -1 timesteps
    0.043124 PDB Coordinate write interval: -1 timesteps
    0.043134 Checkpoint write interval: 10000 timesteps
    0.043145 overlapSepLimit: 0.1 nm
    0.043155 Transition matrix write interval: 10000 timesteps
    0.043166 ClusterOverlapCheck: false
    0.043176 RNGwrite: false
    0.043184 bondedComplexWrite: -1
    0.043195 Restart file format: JSON
    0.043202 Numerical settings:
    0.043213   2D table integration abs/rel error: 1e-07 / 1e-07
    0.043221   2D table integration fallback/tail cutoff: 1e-06 / 1e-10
    0.043230   normalization integration abs/rel error: 1e-06 / 1e-06
    0.043239   table rate abs/rel tolerance: 1e-08 / 0
    0.043249   table diffusion abs/rel tolerance: 0.0001 / 0
    0.043258   explicit/implicit lipid flat-diffusion cutoff: 1e-10 / 1e-15
    0.043269   association same-angle abs/rel tolerance: 0.0001 / 0
    0.043280   association rotation/endpoint-sign tolerance: 1e-08 / 1e-11
    0.043291   Vec3D coordinate precision: 10000
    0.043301 Molecule specific parameters:
    0.043310 Number of unique molecule types: 4
    0.043319 Total number of unique interfaces and states, including product states: 13
    0.043333 Total number of complexes in system at start: 895
    0.043342 Total number of units (molecules + interfaces) in system at start: 2407
    0.043352 Maximum allowed number of unique 2D reactions: 10000
    0.043363  isSphere? false
    0.043374  sphere Radius 0
    0.043383  hasCompartment? false
    0.043393  compartment Radius 0
    0.043401  BOX geometry, dimensions:
    0.043410 125.387 125.387 5470.63
    0.043419  hasImplicitLipid? false
    0.043428 
    0.043439 Molecule Information
    0.043449 Molecule template 0
    0.043459 Name: lipid
    0.043470 Copy number:278
    0.043480  Diffusion trans: 1 1 0
    0.043491  Diffusion Rot: 0 0 0
    0.043499 Is a lipid: true
    0.043508 Is a promoter: false
    0.043517 Is a rod: true
    0.043526 Is a point: false
    0.043535 Radius: 1
    0.043545 
    0.043555 Interfaces:
    0.043565 Name: head
    0.043574 Relative index: 0
    0.043582 States:
    0.043591 (identity: U, absolute index: 0)
    0.043600 
    0.043609 (identity: P, absolute index: 1)
    0.043619 
    0.043629 Molecule template 1
    0.043637 Name: adaptor
    0.043645 Copy number:100
    0.043654  Diffusion trans: 25 25 25
    0.043662  Diffusion Rot: 0.5 0.5 0.5
    0.043671 Is a promoter: false
    0.043680 Is a rod: false
    0.043689 Is a point: false
    0.043699 Radius: 2
    0.043712 
    0.043723 Interfaces:
    0.043733 Name: lbs
    0.043744 Relative index: 0
    0.043754 Absolute index: 2
    0.043762 Forward Reactions:  (4)
    0.043770 
    0.043779 Name: ebs
    0.043788 Relative index: 1
    0.043798 Absolute index: 3
    0.043807 Forward Reactions:  (5) (6)
    0.043816 
    0.043829 Molecule template 2
    0.043839 Name: kinase
    0.043849 Copy number:296
    0.043859  Diffusion trans: 25 25 25
    0.043869  Diffusion Rot: 0.5 0.5 0.5
    0.043878 Is a promoter: false
    0.043886 Is a rod: true
    0.043895 Is a point: false
    0.043903 Radius: 0.5
    0.043912 
    0.043921 Interfaces:
    0.043929 Name: abs
    0.043954 Relative index: 0
    0.043961 Absolute index: 4
    0.043967 Forward Reactions:  (5)
    0.043973 
    0.043979 Name: lbs
    0.043985 Relative index: 1
    0.043991 Absolute index: 5
    0.043996 Forward Reactions:  (0)
    0.044002 
    0.044011 Molecule template 3
    0.044017 Name: phosphatase
    0.044023 Copy number:221
    0.044029  Diffusion trans: 25 25 25
    0.044035  Diffusion Rot: 0.5 0.5 0.5
    0.044041 Is a promoter: false
    0.044046 Is a rod: true
    0.044052 Is a point: false
    0.044058 Radius: 0.5
    0.044064 
    0.044070 Interfaces:
    0.044077 Name: abs
    0.044083 Relative index: 0
    0.044090 Absolute index: 6
    0.044095 Forward Reactions:  (6)
    0.044101 
    0.044107 Name: lbs
    0.044113 Relative index: 1
    0.044119 Absolute index: 7
    0.044125 Forward Reactions:  (2)
    0.044131 
    0.044137 
    0.044143 Reactions
    0.044149 Forward Rxn: 0
    0.044155 Absolute index: 0
    0.044161 Type: Bimolecular association
    0.044168 Reactants:
    0.044173  [iface name: lbs, iface index: 5, required Free]
    0.044179  [iface name: head, iface index: 0, requires state U, required Free]
    0.044185 
    0.044191 Products:
    0.044197  [iface name: lbs, iface index: 8, requires interaction]
    0.044208  [iface name: head, iface index: 8, requires state U, requires interaction]
    0.044216 
    0.044223 Rate(s):
    0.044229 Rate 0: 0.00227955
    0.044235 Reactant 1 requires interfaces:
    0.044242 Reactant 2 requires interfaces:
    0.044248 Sigma: 1
    0.044254 Association angles:
    0.044260 Theta 1    1.5708
    0.044267 Theta 2    1.5708
    0.044273 Phi 1         nan
    0.044280 Phi 2         nan
    0.044286 Omega     3.14159
    0.044292 
    0.044297 label: nonspecified
    0.044303 Is On Membrane: false
    0.044309 bindRadSameCom 1.1
    0.044316 loopCoopFactor 1
    0.044322 length3Dto2D 1
    0.044328 area3Dto1D 12.5664
    0.044334 isCoupled? true
    0.044341  coupledRxn Number: 2 type: Unimolecular state change prob to perform coupled: 0.0263186
    0.044347 
    0.044353 Back Rxn: 0
    0.044359 Absolute index: 1
    0.044365 Type: Bimolecular association
    0.044371 Reactants:
    0.044377  [iface name: lbs, iface index: 8, requires interaction]
    0.044384  [iface name: head, iface index: 8, requires state U, requires interaction]
    0.044390 
    0.044396 Products:
    0.044402  [iface name: lbs, iface index: 5, required Free]
    0.044409  [iface name: head, iface index: 0, requires state U, required Free]
    0.044415 
    0.044421 
    0.044427 Rate(s):
    0.044432 Rate 0: 3.40059
    0.044438 Reactant 1 requires interfaces:
    0.044444 Reactant 2 requires interfaces:
    0.044451 On Membrane? false
    0.044456 Forward Rxn: 1
    0.044462 Absolute index: 2
    0.044469 Type: Unimolecular state change
    0.044479 State Change Reactant: [iface name: head, iface index: 0, requires state U, required Free]
    0.044486 State Change Product: [iface name: head, iface index: 1, requires state P, required Free]
    0.044492 Rate(s):
    0.044498 Rate 0: 0
    0.044504 label: PhosphorylationByKinase
    0.044510 Forward Rxn: 2
    0.044516 Absolute index: 3
    0.044522 Type: Bimolecular association
    0.044528 Reactants:
    0.044534  [iface name: lbs, iface index: 7, required Free]
    0.044540  [iface name: head, iface index: 1, requires state P, required Free]
    0.044546 
    0.044552 Products:
    0.044558  [iface name: lbs, iface index: 9, requires interaction]
    0.044564  [iface name: head, iface index: 9, requires state P, requires interaction]
    0.044570 
    0.044576 Rate(s):
    0.044582 Rate 0: 0.114966
    0.044588 Reactant 1 requires interfaces:
    0.044594 Reactant 2 requires interfaces:
    0.044599 Sigma: 1
    0.044605 Association angles:
    0.044611 Theta 1    1.5708
    0.044617 Theta 2    1.5708
    0.044623 Phi 1         nan
    0.044629 Phi 2         nan
    0.044635 Omega     3.14159
    0.044640 
    0.044646 label: nonspecified
    0.044652 Is On Membrane: false
    0.044658 bindRadSameCom 1.1
    0.044671 loopCoopFactor 1
    0.044677 length3Dto2D 1
    0.044683 area3Dto1D 12.5664
    0.044689 isCoupled? true
    0.044695  coupledRxn Number: 5 type: Unimolecular state change prob to perform coupled: 0.861516
    0.044701 
    0.044711 Back Rxn: 2
    0.044718 Absolute index: 4
    0.044724 Type: Bimolecular association
    0.044730 Reactants:
    0.044737  [iface name: lbs, iface index: 9, requires interaction]
    0.044743  [iface name: head, iface index: 9, requires state P, requires interaction]
    0.044749 
    0.044755 Products:
    0.044761  [iface name: lbs, iface index: 7, required Free]
    0.044767  [iface name: head, iface index: 1, requires state P, required Free]
    0.044773 
    0.044779 
    0.044785 Rate(s):
    0.044792 Rate 0: 0.0320455
    0.044798 Reactant 1 requires interfaces:
    0.044804 Reactant 2 requires interfaces:
    0.044810 On Membrane? false
    0.044816 Forward Rxn: 3
    0.044822 Absolute index: 5
    0.044828 Type: Unimolecular state change
    0.044834 State Change Reactant: [iface name: head, iface index: 1, requires state P, required Free]
    0.044840 State Change Product: [iface name: head, iface index: 0, requires state U, required Free]
    0.044846 Rate(s):
    0.044852 Rate 0: 0
    0.044858 label: DephosphorylationByPhosphatase
    0.044863 Forward Rxn: 4
    0.044870 Absolute index: 6
    0.044876 Type: Bimolecular association
    0.044882 Reactants:
    0.044887  [iface name: head, iface index: 1, requires state P, required Free]
    0.044893  [iface name: lbs, iface index: 2, required Free]
    0.044899 
    0.044905 Products:
    0.044915  [iface name: head, iface index: 10, requires state P, requires interaction]
    0.044922  [iface name: lbs, iface index: 10, requires interaction]
    0.044928 
    0.044934 Rate(s):
    0.044940 Rate 0: 1.95636
    0.044946 Reactant 1 requires interfaces:
    0.044952 Reactant 2 requires interfaces:
    0.044958 Sigma: 1
    0.044964 Association angles:
    0.044970 Theta 1    1.5708
    0.044977 Theta 2    1.5708
    0.044983 Phi 1         nan
    0.044989 Phi 2         nan
    0.044995 Omega     3.14159
    0.045002 
    0.045008 label: nonspecified
    0.045015 Is On Membrane: false
    0.045022 bindRadSameCom 1.1
    0.045028 loopCoopFactor 1
    0.045038 length3Dto2D 1
    0.045045 area3Dto1D 12.5664
    0.045052 isCoupled? false
    0.045059 
    0.045067 Back Rxn: 4
    0.045074 Absolute index: 7
    0.045081 Type: Bimolecular association
    0.045090 Reactants:
    0.045100  [iface name: head, iface index: 10, requires state P, requires interaction]
    0.045110  [iface name: lbs, iface index: 10, requires interaction]
    0.045120 
    0.045127 Products:
    0.045135  [iface name: head, iface index: 1, requires state P, required Free]
    0.045142  [iface name: lbs, iface index: 2, required Free]
    0.045152 
    0.045161 
    0.045168 Rate(s):
    0.045176 Rate 0: 0.0876411
    0.045186 Reactant 1 requires interfaces:
    0.045196 Reactant 2 requires interfaces:
    0.045204 On Membrane? false
    0.045213 Forward Rxn: 5
    0.045220 Absolute index: 8
    0.045230 Type: Bimolecular association
    0.045237 Reactants:
    0.045251  [iface name: abs, iface index: 4, required Free]
    0.045258  [iface name: ebs, iface index: 3, required Free]
    0.045265 
    0.045274 Products:
    0.045284  [iface name: abs, iface index: 11, requires interaction]
    0.045292  [iface name: ebs, iface index: 11, requires interaction]
    0.045300 
    0.045309 Rate(s):
    0.045318 Rate 0: 9.25751
    0.045325 Reactant 1 requires interfaces:
    0.045334 Reactant 2 requires interfaces:
    0.045342  [iface name: lbs, iface index: 2, requires interaction]
    0.045353 Sigma: 1
    0.045362 Association angles:
    0.045371 Theta 1    1.5708
    0.045381 Theta 2   2.35619
    0.045389 Phi 1         nan
    0.045398 Phi 2     3.14159
    0.045407 Omega     3.14159
    0.045415 
    0.045424 label: nonspecified
    0.045434 Is On Membrane: false
    0.045444 bindRadSameCom 1.1
    0.045453 loopCoopFactor 1
    0.045464 length3Dto2D 1
    0.045524 area3Dto1D 12.5664
    0.045533 isCoupled? false
    0.045541 
    0.045550 Back Rxn: 5
    0.045559 Absolute index: 9
    0.045567 Type: Bimolecular association
    0.045573 Reactants:
    0.045579  [iface name: abs, iface index: 11, requires interaction]
    0.045586  [iface name: ebs, iface index: 11, requires interaction]
    0.045592 
    0.045599 Products:
    0.045605  [iface name: abs, iface index: 4, required Free]
    0.045611  [iface name: ebs, iface index: 3, required Free]
    0.045618 
    0.045624 
    0.045630 Rate(s):
    0.045636 Rate 0: 0.00621957
    0.045643 Reactant 1 requires interfaces:
    0.045654 Reactant 2 requires interfaces:
    0.045663  [iface name: lbs, iface index: 2, requires interaction]
    0.045669 On Membrane? false
    0.045675 Forward Rxn: 6
    0.045682 Absolute index: 10
    0.045689 Type: Bimolecular association
    0.045695 Reactants:
    0.045702  [iface name: abs, iface index: 6, required Free]
    0.045709  [iface name: ebs, iface index: 3, required Free]
    0.045716 
    0.045723 Products:
    0.045731  [iface name: abs, iface index: 12, requires interaction]
    0.045737  [iface name: ebs, iface index: 12, requires interaction]
    0.045744 
    0.045751 Rate(s):
    0.045758 Rate 0: 0.254622
    0.045764 Reactant 1 requires interfaces:
    0.045771 Reactant 2 requires interfaces:
    0.045777  [iface name: lbs, iface index: 2, requires interaction]
    0.045784 Sigma: 1
    0.045790 Association angles:
    0.045796 Theta 1    1.5708
    0.045803 Theta 2   2.35619
    0.045809 Phi 1         nan
    0.045816 Phi 2     3.14159
    0.045822 Omega     3.14159
    0.045828 
    0.045835 label: nonspecified
    0.045841 Is On Membrane: false
    0.045847 bindRadSameCom 1.1
    0.045853 loopCoopFactor 1
    0.045860 length3Dto2D 1
    0.045866 area3Dto1D 12.5664
    0.045872 isCoupled? false
    0.045878 
    0.045884 Back Rxn: 6
    0.045891 Absolute index: 11
    0.045897 Type: Bimolecular association
    0.045903 Reactants:
    0.045909  [iface name: abs, iface index: 12, requires interaction]
    0.045916  [iface name: ebs, iface index: 12, requires interaction]
    0.045922 
    0.045928 Products:
    0.045939  [iface name: abs, iface index: 6, required Free]
    0.045946  [iface name: ebs, iface index: 3, required Free]
    0.046009 
    0.046020 
    0.046030 Rate(s):
    0.046050 Rate 0: 0.00108422
    0.046060 Reactant 1 requires interfaces:
    0.046067 Reactant 2 requires interfaces:
    0.046082  [iface name: lbs, iface index: 2, requires interaction]
    0.046104 On Membrane? false
    0.058742 *************** BEGIN SIMULATION ****************
    7.834207 -------------------
    7.835548 End iteration: 1000, simulation time: 1.000000e-03 seconds.
    7.835925 System time: 2026-10-14 13:36:16
    7.835944 Elapsed time: 0 minutes
    7.835961 Number of molecules: 895
    7.835970 Number of complexes: 895
    7.835977 Total reaction matches: 157092
    7.836082 Avg timestep duration: 7.775401e-03, iterations remaining: 99000, Time left: 12 minutes
    7.836096 Estimated end time: 2026-10-14 13:49:06
    7.836113 -------------------------------------------------
   13.617520 -------------------
   13.617685 End iteration: 2000, simulation time: 2.000000e-03 seconds.
   13.618106 System time: 2026-10-14 13:36:22
   13.618114 Elapsed time: 0 minutes
   13.618120 Number of molecules: 895
   13.618126 Number of complexes: 895
   13.618157 Total reaction matches: 421948
   13.618176 Avg timestep duration: 5.781311e-03, iterations remaining: 98000, Time left: 9 minutes
   13.618183 Estimated end time: 2026-10-14 13:45:49
   13.618189 -------------------------------------------------
   18.286841 -------------------
   18.286879 End iteration: 3000, simulation time: 3.000000e-03 seconds.
   18.286951 System time: 2026-10-14 13:36:27
   18.286962 Elapsed time: 0 minutes
   18.286967 Number of molecules: 895
   18.286972 Number of complexes: 895
   18.286978 Total reaction matches: 700125
   18.286987 Avg timestep duration: 4.668556e-03, iterations remaining: 97000, Time left: 7 minutes
   18.287009 Estimated end time: 2026-10-14 13:43:59
   18.287013 -------------------------------------------------
   22.979467 -------------------
   22.980763 End iteration: 4000, simulation time: 4.000000e-03 seconds.
   22.981206 System time: 2026-10-14 13:36:31
   22.981221 Elapsed time: 0 minutes
   22.981228 Number of molecules: 895
   22.981237 Number of complexes: 895
   22.981243 Total reaction matches: 1023484
   22.981259 Avg timestep duration: 4.692410e-03, iterations remaining: 96000, Time left: 7 minutes
   22.981266 Estimated end time: 2026-10-14 13:44:02
   22.981278 -------------------------------------------------
   27.746093 -------------------
   27.746155 End iteration: 5000, simulation time: 5.000000e-03 seconds.
   27.746177 System time: 2026-10-14 13:36:36
   27.746195 Elapsed time: 0 minutes
   27.746202 Number of molecules: 895
   27.746209 Number of complexes: 895
   27.746216 Total reaction matches: 1388515
   27.746221 Avg timestep duration: 4.764625e-03, iterations remaining: 95000, Time left: 7 minutes
   27.746226 Estimated end time: 2026-10-14 13:44:09
   27.746231 -------------------------------------------------
   32.482265 -------------------
   32.482306 End iteration: 6000, simulation time: 6.000000e-03 seconds.
   32.482851 System time: 2026-10-14 13:36:41
   32.482871 Elapsed time: 0 minutes
   32.482878 Number of molecules: 895
   32.482883 Number of complexes: 895
   32.482888 Total reaction matches: 1591285
   32.482987 Avg timestep duration: 4.736025e-03, iterations remaining: 94000, Time left: 7 minutes
   32.482997 Estimated end time: 2026-10-14 13:44:06
   32.483002 -------------------------------------------------
   37.134614 -------------------
   37.134644 End iteration: 7000, simulation time: 7.000000e-03 seconds.
   37.134709 System time: 2026-10-14 13:36:45
   37.134720 Elapsed time: 0 minutes
   37.134725 Number of molecules: 895
   37.134730 Number of complexes: 895
   37.134735 Total reaction matches: 1804703
   37.134744 Avg timestep duration: 4.651543e-03, iterations remaining: 93000, Time left: 7 minutes
   37.134750 Estimated end time: 2026-10-14 13:43:58
   37.134755 -------------------------------------------------
   43.404575 -------------------
   43.404600 End iteration: 8000, simulation time: 8.000000e-03 seconds.
   43.405043 System time: 2026-10-14 13:36:52
   43.405063 Elapsed time: 0 minutes
   43.405071 Number of molecules: 895
   43.405077 Number of complexes: 895
   43.405084 Total reaction matches: 2120792
   43.405153 Avg timestep duration: 6.269769e-03, iterations remaining: 92000, Time left: 9 minutes
   43.405162 Estimated end time: 2026-10-14 13:46:29
   43.405168 -------------------------------------------------
   51.501134 -------------------
   51.501164 End iteration: 9000, simulation time: 9.000000e-03 seconds.
   51.501631 System time: 2026-10-14 13:37:00
   51.501652 Elapsed time: 0 minutes
   51.501665 Number of molecules: 895
   51.501677 Number of complexes: 895
   51.501689 Total reaction matches: 2388709
   51.501744 Avg timestep duration: 8.095923e-03, iterations remaining: 91000, Time left: 12 minutes
   51.501761 Estimated end time: 2026-10-14 13:49:17
   51.501777 -------------------------------------------------
   59.782799 -------------------
   59.782827 End iteration: 10000, simulation time: 1.000000e-02 seconds.
   59.783274 System time: 2026-10-14 13:37:08
   59.783300 Elapsed time: 0 minutes
   59.783313 Number of molecules: 895
   59.783331 Number of complexes: 895
   59.783363 Total reaction matches: 2764706
   59.783406 Avg timestep duration: 8.280951e-03, iterations remaining: 90000, Time left: 12 minutes
   59.783426 Estimated end time: 2026-10-14 13:49:33
   59.783441 -------------------------------------------------
   68.162217 -------------------
   68.162253 End iteration: 11000, simulation time: 1.100000e-02 seconds.
   68.162828 System time: 2026-10-14 13:37:17
   68.162844 Elapsed time: 1 minutes
   68.162852 Number of molecules: 895
   68.162858 Number of complexes: 895
   68.162865 Total reaction matches: 3077891
   68.162947 Avg timestep duration: 8.378741e-03, iterations remaining: 89000, Time left: 12 minutes
   68.162954 Estimated end time: 2026-10-14 13:49:42
   68.162960 -------------------------------------------------
   76.280554 -------------------
   76.280760 End iteration: 12000, simulation time: 1.200000e-02 seconds.
   76.281113 System time: 2026-10-14 13:37:25
   76.281137 Elapsed time: 1 minutes
   76.281158 Number of molecules: 895
   76.281175 Number of complexes: 895
   76.281190 Total reaction matches: 3233680
   76.281214 Avg timestep duration: 8.117532e-03, iterations remaining: 88000, Time left: 11 minutes
   76.281231 Estimated end time: 2026-10-14 13:49:19
   76.281261 -------------------------------------------------
   84.438054 -------------------
   84.438083 End iteration: 13000, simulation time: 1.300000e-02 seconds.
   84.439751 System time: 2026-10-14 13:37:33
   84.439766 Elapsed time: 1 minutes
   84.439778 Number of molecules: 895
   84.439789 Number of complexes: 895
   84.439808 Total reaction matches: 3528293
   84.439825 Avg timestep duration: 8.156808e-03, iterations remaining: 87000, Time left: 11 minutes
   84.439838 Estimated end time: 2026-10-14 13:49:22
   84.439850 -------------------------------------------------
   92.729797 -------------------
   92.731012 End iteration: 14000, simulation time: 1.400000e-02 seconds.
   92.731345 System time: 2026-10-14 13:37:41
   92.731366 Elapsed time: 1 minutes
   92.731373 Number of molecules: 895
   92.731385 Number of complexes: 895
   92.731392 Total reaction matches: 3866336
   92.731406 Avg timestep duration: 8.289901e-03, iterations remaining: 86000, Time left: 11 minutes
   92.731413 Estimated end time: 2026-10-14 13:49:34
   92.731419 -------------------------------------------------
  100.747430 -------------------
  100.747611 End iteration: 15000, simulation time: 1.500000e-02 seconds.
  100.747968 System time: 2026-10-14 13:37:49
  100.747991 Elapsed time: 1 minutes
  100.748006 Number of molecules: 895
  100.748023 Number of complexes: 895
  100.748041 Total reaction matches: 4111294
  100.748063 Avg timestep duration: 8.015953e-03, iterations remaining: 85000, Time left: 11 minutes
  100.748075 Estimated end time: 2026-10-14 13:49:10
  100.748087 -------------------------------------------------
  108.771975 -------------------
  108.772168 End iteration: 16000, simulation time: 1.600000e-02 seconds.
  108.772518 System time: 2026-10-14 13:37:57
  108.772535 Elapsed time: 1 minutes
  108.772543 Number of molecules: 895
  108.772549 Number of complexes: 895
  108.772560 Total reaction matches: 4228316
  108.772585 Avg timestep duration: 8.023883e-03, iterations remaining: 84000, Time left: 11 minutes
  108.772594 Estimated end time: 2026-10-14 13:49:11
  108.772602 -------------------------------------------------
  116.725923 -------------------
  116.725951 End iteration: 17000, simulation time: 1.700000e-02 seconds.
  116.726020 System time: 2026-10-14 13:38:05
  116.726034 Elapsed time: 1 minutes
  116.726042 Number of molecules: 895
  116.726048 Number of complexes: 895
  116.726054 Total reaction matches: 4389989
  116.726066 Avg timestep duration: 7.953264e-03, iterations remaining: 83000, Time left: 11 minutes
  116.726073 Estimated end time: 2026-10-14 13:49:05
  116.726079 -------------------------------------------------
  124.990136 -------------------
  124.990329 End iteration: 18000, simulation time: 1.800000e-02 seconds.
  124.990660 System time: 2026-10-14 13:38:13
  124.990678 Elapsed time: 2 minutes
  124.990685 Number of molecules: 895
  124.990692 Number of complexes: 895
  124.990698 Total reaction matches: 4854619
  124.990734 Avg timestep duration: 8.264035e-03, iterations remaining: 82000, Time left: 11 minutes
  124.990742 Estimated end time: 2026-10-14 13:49:31
  124.990749 -------------------------------------------------
  133.292090 -------------------
  133.292268 End iteration: 19000, simulation time: 1.900000e-02 seconds.
  133.292621 System time: 2026-10-14 13:38:22
  133.292676 Elapsed time: 2 minutes
  133.292686 Number of molecules: 895
  133.292693 Number of complexes: 895
  133.292700 Total reaction matches: 5220348
  133.292723 Avg timestep duration: 8.301277e-03, iterations remaining: 81000, Time left: 11 minutes
  133.292730 Estimated end time: 2026-10-14 13:49:34
  133.292736 -------------------------------------------------
  141.359080 -------------------
  141.359103 End iteration: 20000, simulation time: 2.000000e-02 seconds.
  141.359323 System time: 2026-10-14 13:38:30
  141.359337 Elapsed time: 2 minutes
  141.359346 Number of molecules: 895
  141.359356 Number of complexes: 895
  141.359364 Total reaction matches: 5379417
  141.359378 Avg timestep duration: 8.066300e-03, iterations remaining: 80000, Time left: 10 minutes
  141.359387 Estimated end time: 2026-10-14 13:49:15
  141.359397 -------------------------------------------------
  149.469587 -------------------
  149.470514 End iteration: 21000, simulation time: 2.100000e-02 seconds.
  149.470838 System time: 2026-10-14 13:38:38
  149.470853 Elapsed time: 2 minutes
  149.470865 Number of molecules: 895
  149.470879 Number of complexes: 895
  149.470894 Total reaction matches: 5547560
  149.470909 Avg timestep duration: 8.110172e-03, iterations remaining: 79000, Time left: 10 minutes
  149.470917 Estimated end time: 2026-10-14 13:49:19
  149.470923 -------------------------------------------------
  157.470467 -------------------
  157.470492 End iteration: 22000, simulation time: 2.200000e-02 seconds.
  157.470637 System time: 2026-10-14 13:38:46
  157.470898 Elapsed time: 2 minutes
  157.470915 Number of molecules: 895
  157.470922 Number of complexes: 895
  157.470928 Total reaction matches: 5734094
  157.470950 Avg timestep duration: 7.999494e-03, iterations remaining: 78000, Time left: 10 minutes
  157.470957 Estimated end time: 2026-10-14 13:49:10
  157.470963 -------------------------------------------------
  163.623825 -------------------
  163.623940 End iteration: 23000, simulation time: 2.300000e-02 seconds.
  163.624367 System time: 2026-10-14 13:38:52
  163.624392 Elapsed time: 2 minutes
  163.624400 Number of molecules: 895
  163.624408 Number of complexes: 895
  163.624416 Total reaction matches: 5872138
  163.624430 Avg timestep duration: 6.152854e-03, iterations remaining: 77000, Time left: 7 minutes
  163.624437 Estimated end time: 2026-10-14 13:46:46
  163.624444 -------------------------------------------------
  168.307279 -------------------
  168.307315 End iteration: 24000, simulation time: 2.400000e-02 seconds.
  168.307361 System time: 2026-10-14 13:38:57
  168.307372 Elapsed time: 2 minutes
  168.307378 Number of molecules: 895
  168.307383 Number of complexes: 895
  168.307388 Total reaction matches: 6131902
  168.307397 Avg timestep duration: 4.682767e-03, iterations remaining: 76000, Time left: 5 minutes
  168.307403 Estimated end time: 2026-10-14 13:44:53
  168.307408 -------------------------------------------------
  172.937348 -------------------
  172.937586 End iteration: 25000, simulation time: 2.500000e-02 seconds.
  172.938008 System time: 2026-10-14 13:39:01
  172.938019 Elapsed time: 2 minutes
  172.938025 Number of molecules: 895
  172.938030 Number of complexes: 895
  172.938035 Total reaction matches: 6340431
  172.938096 Avg timestep duration: 4.629852e-03, iterations remaining: 75000, Time left: 5 minutes
  172.938104 Estimated end time: 2026-10-14 13:44:49
  172.938110 -------------------------------------------------
  177.706633 -------------------
  177.706662 End iteration: 26000, simulation time: 2.600000e-02 seconds.
  177.706726 System time: 2026-10-14 13:39:06
  177.706737 Elapsed time: 2 minutes
  177.706743 Number of molecules: 895
  177.706748 Number of complexes: 895
  177.706753 Total reaction matches: 6537930
  177.706762 Avg timestep duration: 4.768474e-03, iterations remaining: 74000, Time left: 5 minutes
  177.706768 Estimated end time: 2026-10-14 13:44:59
  177.706772 -------------------------------------------------
  182.462563 -------------------
  182.462781 End iteration: 27000, simulation time: 2.700000e-02 seconds.
  182.463193 System time: 2026-10-14 13:39:11
  182.463208 Elapsed time: 3 minutes
  182.463214 Number of molecules: 895
  182.463219 Number of complexes: 895
  182.463224 Total reaction matches: 6718849
  182.463274 Avg timestep duration: 4.755705e-03, iterations remaining: 73000, Time left: 5 minutes
  182.463282 Estimated end time: 2026-10-14 13:44:58
  182.463287 -------------------------------------------------
  187.266682 -------------------
  187.266710 End iteration: 28000, simulation time: 2.800000e-02 seconds.
  187.266772 System time: 2026-10-14 13:39:16
  187.266781 Elapsed time: 3 minutes
  187.266786 Number of molecules: 895
  187.266790 Number of complexes: 895
  187.266794 Total reaction matches: 6920102
  187.266802 Avg timestep duration: 4.803330e-03, iterations remaining: 72000, Time left: 5 minutes
  187.266807 Estimated end time: 2026-10-14 13:45:01
  187.266811 -------------------------------------------------
  193.132230 -------------------
  193.132449 End iteration: 29000, simulation time: 2.900000e-02 seconds.
  193.132781 System time: 2026-10-14 13:39:21
  193.132801 Elapsed time: 3 minutes
  193.132811 Number of molecules: 895
  193.132821 Number of complexes: 895
  193.132831 Total reaction matches: 7139267
  193.132890 Avg timestep duration: 5.865354e-03, iterations remaining: 71000, Time left: 6 minutes
  193.132903 Estimated end time: 2026-10-14 13:46:18
  193.132912 -------------------------------------------------
  200.801699 -------------------
  200.801727 End iteration: 30000, simulation time: 3.000000e-02 seconds.
  200.801985 System time: 2026-10-14 13:39:29
  200.802006 Elapsed time: 3 minutes
  200.802019 Number of molecules: 895
  200.802029 Number of complexes: 895
  200.802040 Total reaction matches: 7460408
  200.802058 Avg timestep duration: 7.668684e-03, iterations remaining: 70000, Time left: 8 minutes
  200.802071 Estimated end time: 2026-10-14 13:48:26
  200.802081 -------------------------------------------------
  205.911892 -------------------
  205.912141 End iteration: 31000, simulation time: 3.100000e-02 seconds.
  205.912610 System time: 2026-10-14 13:39:34
  205.912626 Elapsed time: 3 minutes
  205.912634 Number of molecules: 895
  205.912639 Number of complexes: 895
  205.912644 Total reaction matches: 7909921
  205.912689 Avg timestep duration: 5.109791e-03, iterations remaining: 69000, Time left: 5 minutes
  205.912705 Estimated end time: 2026-10-14 13:45:27
  205.912712 -------------------------------------------------
  210.616183 -------------------
  210.616230 End iteration: 32000, simulation time: 3.200000e-02 seconds.
  210.616272 System time: 2026-10-14 13:39:39
  210.616284 Elapsed time: 3 minutes
  210.616291 Number of molecules: 895
  210.616296 Number of complexes: 895
  210.616302 Total reaction matches: 8187105
  210.616312 Avg timestep duration: 4.703401e-03, iterations remaining: 68000, Time left: 5 minutes
  210.616319 Estimated end time: 2026-10-14 13:44:59
  210.616324 -------------------------------------------------
  215.386247 -------------------
  215.386519 End iteration: 33000, simulation time: 3.300000e-02 seconds.
  215.386936 System time: 2026-10-14 13:39:44
  215.386948 Elapsed time: 3 minutes
  215.386954 Number of molecules: 895
  215.386965 Number of complexes: 895
  215.386970 Total reaction matches: 8496237
  215.387069 Avg timestep duration: 4.769837e-03, iterations remaining: 67000, Time left: 5 minutes
  215.387078 Estimated end time: 2026-10-14 13:45:03
  215.387083 -------------------------------------------------
  220.053529 -------------------
  220.053572 End iteration: 34000, simulation time: 3.400000e-02 seconds.
  220.053623 System time: 2026-10-14 13:39:48
  220.053633 Elapsed time: 3 minutes
  220.053639 Number of molecules: 895
  220.053644 Number of complexes: 895
  220.053649 Total reaction matches: 8803570
  220.053661 Avg timestep duration: 4.666403e-03, iterations remaining: 66000, Time left: 5 minutes
  220.053667 Estimated end time: 2026-10-14 13:44:56
  220.053691 -------------------------------------------------
  224.664219 -------------------
  224.664460 End iteration: 35000, simulation time: 3.500000e-02 seconds.
  224.664756 System time: 2026-10-14 13:39:53
  224.664769 Elapsed time: 3 minutes
  224.664775 Number of molecules: 895
  224.664781 Number of complexes: 895
  224.664786 Total reaction matches: 9018350
  224.664828 Avg timestep duration: 4.610523e-03, iterations remaining: 65000, Time left: 4 minutes
  224.664836 Estimated end time: 2026-10-14 13:44:53
  224.664842 -------------------------------------------------
  229.295713 -------------------
  229.295738 End iteration: 36000, simulation time: 3.600000e-02 seconds.
  229.295810 System time: 2026-10-14 13:39:58
  229.295821 Elapsed time: 3 minutes
  229.295827 Number of molecules: 895
  229.295832 Number of complexes: 895
  229.295836 Total reaction matches: 9284145
  229.295846 Avg timestep duration: 4.630833e-03, iterations remaining: 64000, Time left: 4 minutes
  229.295851 Estimated end time: 2026-10-14 13:44:54
  229.295856 -------------------------------------------------
  233.973241 -------------------
  233.973464 End iteration: 37000, simulation time: 3.700000e-02 seconds.
  233.973762 System time: 2026-10-14 13:40:02
  233.973774 Elapsed time: 3 minutes
  233.973780 Number of molecules: 895
  233.973786 Number of complexes: 895
  233.973791 Total reaction matches: 9709087
  233.973868 Avg timestep duration: 4.677349e-03, iterations remaining: 63000, Time left: 4 minutes
  233.973876 Estimated end time: 2026-10-14 13:44:57
  233.973881 -------------------------------------------------
  238.582800 -------------------
  238.582827 End iteration: 38000, simulation time: 3.800000e-02 seconds.
  238.582889 System time: 2026-10-14 13:40:07
  238.582900 Elapsed time: 3 minutes
  238.582906 Number of molecules: 895
  238.582912 Number of complexes: 895
  238.582917 Total reaction matches: 9921366
  238.582926 Avg timestep duration: 4.608870e-03, iterations remaining: 62000, Time left: 4 minutes
  238.582932 Estimated end time: 2026-10-14 13:44:53
  238.582937 -------------------------------------------------
  243.369809 -------------------
  243.370018 End iteration: 39000, simulation time: 3.900000e-02 seconds.
  243.370369 System time: 2026-10-14 13:40:12
  243.370384 Elapsed time: 4 minutes
  243.370391 Number of molecules: 895
  243.370397 Number of complexes: 895
  243.370402 Total reaction matches: 10187257
  243.370472 Avg timestep duration: 4.786835e-03, iterations remaining: 61000, Time left: 4 minutes
  243.370480 Estimated end time: 2026-10-14 13:45:04
  243.370485 -------------------------------------------------
  248.010144 -------------------
  248.010171 End iteration: 40000, simulation time: 4.000000e-02 seconds.
  248.010268 System time: 2026-10-14 13:40:16
  248.010278 Elapsed time: 4 minutes
  248.010284 Number of molecules: 895
  248.010290 Number of complexes: 895
  248.010295 Total reaction matches: 10446902
  248.010304 Avg timestep duration: 4.639608e-03, iterations remaining: 60000, Time left: 4 minutes
  248.010310 Estimated end time: 2026-10-14 13:44:55
  248.010315 -------------------------------------------------
  252.720709 -------------------
  252.720936 End iteration: 41000, simulation time: 4.100000e-02 seconds.
  252.721189 System time: 2026-10-14 13:40:21
  252.721204 Elapsed time: 4 minutes
  252.721214 Number of molecules: 895
  252.721219 Number of complexes: 895
  252.721224 Total reaction matches: 10715376
  252.721297 Avg timestep duration: 4.710369e-03, iterations remaining: 59000, Time left: 4 minutes
  252.721305 Estimated end time: 2026-10-14 13:44:59
  252.721310 -------------------------------------------------
  257.318700 -------------------
  257.318725 End iteration: 42000, simulation time: 4.200000e-02 seconds.
  257.318793 System time: 2026-10-14 13:40:26
  257.318803 Elapsed time: 4 minutes
  257.318809 Number of molecules: 895
  257.318814 Number of complexes: 895
  257.318819 Total reaction matches: 10952960
  257.318857 Avg timestep duration: 4.597328e-03, iterations remaining: 58000, Time left: 4 minutes
  257.318862 Estimated end time: 2026-10-14 13:44:52
  257.318866 -------------------------------------------------
  261.942629 -------------------
  261.942833 End iteration: 43000, simulation time: 4.300000e-02 seconds.
  261.943297 System time: 2026-10-14 13:40:30
  261.943311 Elapsed time: 4 minutes
  261.943318 Number of molecules: 895
  261.943323 Number of complexes: 895
  261.943328 Total reaction matches: 11207594
  261.943371 Avg timestep duration: 4.623763e-03, iterations remaining: 57000, Time left: 4 minutes
  261.943382 Estimated end time: 2026-10-14 13:44:54
  261.943388 -------------------------------------------------
  266.690418 -------------------
  266.690443 End iteration: 44000, simulation time: 4.400000e-02 seconds.
  266.690508 System time: 2026-10-14 13:40:35
  266.690519 Elapsed time: 4 minutes
  266.690525 Number of molecules: 895
  266.690530 Number of complexes: 895
  266.690534 Total reaction matches: 11481081
  266.690544 Avg timestep duration: 4.746972e-03, iterations remaining: 56000, Time left: 4 minutes
  266.690549 Estimated end time: 2026-10-14 13:45:01
  266.690554 -------------------------------------------------
  271.369002 -------------------
  271.369190 End iteration: 45000, simulation time: 4.500000e-02 seconds.
  271.369487 System time: 2026-10-14 13:40:40
  271.369498 Elapsed time: 4 minutes
  271.369504 Number of molecules: 895
  271.369509 Number of complexes: 895
  271.369514 Total reaction matches: 11780117
  271.369558 Avg timestep duration: 4.678411e-03, iterations remaining: 55000, Time left: 4 minutes
  271.369567 Estimated end time: 2026-10-14 13:44:57
  271.369572 -------------------------------------------------
  276.060530 -------------------
  276.060554 End iteration: 46000, simulation time: 4.600000e-02 seconds.
  276.060624 System time: 2026-10-14 13:40:44
  276.060634 Elapsed time: 4 minutes
  276.060640 Number of molecules: 895
  276.060645 Number of complexes: 895
  276.060650 Total reaction matches: 12015448
  276.060659 Avg timestep duration: 4.690908e-03, iterations remaining: 54000, Time left: 4 minutes
  276.060664 Estimated end time: 2026-10-14 13:44:58
  276.060669 -------------------------------------------------
  280.714554 -------------------
  280.714754 End iteration: 47000, simulation time: 4.700000e-02 seconds.
  280.715063 System time: 2026-10-14 13:40:49
  280.715073 Elapsed time: 4 minutes
  280.715079 Number of molecules: 895
  280.715084 Number of complexes: 895
  280.715089 Total reaction matches: 12361747
  280.715128 Avg timestep duration: 4.653850e-03, iterations remaining: 53000, Time left: 4 minutes
  280.715136 Estimated end time: 2026-10-14 13:44:56
  280.715141 -------------------------------------------------
  285.384358 -------------------
  285.384386 End iteration: 48000, simulation time: 4.800000e-02 seconds.
  285.384585 System time: 2026-10-14 13:40:54
  285.384596 Elapsed time: 4 minutes
  285.384602 Number of molecules: 895
  285.384606 Number of complexes: 895
  285.384612 Total reaction matches: 12656629
  285.384616 Avg timestep duration: 4.669132e-03, iterations remaining: 52000, Time left: 4 minutes
  285.384621 Estimated end time: 2026-10-14 13:44:57
  285.384627 -------------------------------------------------
  290.020529 -------------------
  290.020712 End iteration: 49000, simulation time: 4.900000e-02 seconds.
  290.021036 System time: 2026-10-14 13:40:58
  290.021048 Elapsed time: 4 minutes
  290.021053 Number of molecules: 895
  290.021058 Number of complexes: 895
  290.021063 Total reaction matches: 12769121
  290.021109 Avg timestep duration: 4.635895e-03, iterations remaining: 51000, Time left: 3 minutes
  290.021117 Estimated end time: 2026-10-14 13:44:55
  290.021122 -------------------------------------------------
  294.710840 -------------------
  294.710866 End iteration: 50000, simulation time: 5.000000e-02 seconds.
  294.710985 System time: 2026-10-14 13:41:03
  294.711004 Elapsed time: 4 minutes
  294.711025 Number of molecules: 895
  294.711032 Number of complexes: 895
  294.711036 Total reaction matches: 12978474
  294.711040 Avg timestep duration: 4.689671e-03, iterations remaining: 50000, Time left: 3 minutes
  294.711052 Estimated end time: 2026-10-14 13:44:58
  294.711057 -------------------------------------------------
  299.413363 -------------------
  299.413584 End iteration: 51000, simulation time: 5.100000e-02 seconds.
  299.413882 System time: 2026-10-14 13:41:08
  299.413892 Elapsed time: 4 minutes
  299.413898 Number of molecules: 895
  299.413903 Number of complexes: 895
  299.413907 Total reaction matches: 13237541
  299.413951 Avg timestep duration: 4.702323e-03, iterations remaining: 49000, Time left: 3 minutes
  299.413959 Estimated end time: 2026-10-14 13:44:58
  299.413964 -------------------------------------------------
  304.274850 -------------------
  304.274881 End iteration: 52000, simulation time: 5.200000e-02 seconds.
  304.274947 System time: 2026-10-14 13:41:13
  304.274957 Elapsed time: 5 minutes
  304.274962 Number of molecules: 895
  304.274967 Number of complexes: 895
  304.274972 Total reaction matches: 13667069
  304.274981 Avg timestep duration: 4.860840e-03, iterations remaining: 48000, Time left: 3 minutes
  304.274986 Estimated end time: 2026-10-14 13:45:06
  304.274991 -------------------------------------------------
  308.928547 -------------------
  308.928732 End iteration: 53000, simulation time: 5.300000e-02 seconds.
  308.929054 System time: 2026-10-14 13:41:17
  308.929065 Elapsed time: 5 minutes
  308.929071 Number of molecules: 895
  308.929077 Number of complexes: 895
  308.929082 Total reaction matches: 13951664
  308.929130 Avg timestep duration: 4.653517e-03, iterations remaining: 47000, Time left: 3 minutes
  308.929147 Estimated end time: 2026-10-14 13:44:56
  308.929155 -------------------------------------------------
  313.582979 -------------------
  313.583007 End iteration: 54000, simulation time: 5.400000e-02 seconds.
  313.583069 System time: 2026-10-14 13:41:22
  313.583080 Elapsed time: 5 minutes
  313.583085 Number of molecules: 895
  313.583091 Number of complexes: 895
  313.583095 Total reaction matches: 14133506
  313.583104 Avg timestep duration: 4.653801e-03, iterations remaining: 46000, Time left: 3 minutes
  313.583110 Estimated end time: 2026-10-14 13:44:56
  313.583115 -------------------------------------------------
  318.248103 -------------------
  318.248290 End iteration: 55000, simulation time: 5.500000e-02 seconds.
  318.248631 System time: 2026-10-14 13:41:27
  318.248644 Elapsed time: 5 minutes
  318.248649 Number of molecules: 895
  318.248654 Number of complexes: 895
  318.248659 Total reaction matches: 14366437
  318.248703 Avg timestep duration: 4.664969e-03, iterations remaining: 45000, Time left: 3 minutes
  318.248711 Estimated end time: 2026-10-14 13:44:57
  318.248716 -------------------------------------------------
  322.812968 -------------------
  322.812992 End iteration: 56000, simulation time: 5.600000e-02 seconds.
  322.813066 System time: 2026-10-14 13:41:31
  322.813076 Elapsed time: 5 minutes
  322.813082 Number of molecules: 895
  322.813086 Number of complexes: 895
  322.813091 Total reaction matches: 14530544
  322.813100 Avg timestep duration: 4.564187e-03, iterations remaining: 44000, Time left: 3 minutes
  322.813106 Estimated end time: 2026-10-14 13:44:52
  322.813110 -------------------------------------------------
  327.452360 -------------------
  327.452546 End iteration: 57000, simulation time: 5.700000e-02 seconds.
  327.452875 System time: 2026-10-14 13:41:36
  327.452890 Elapsed time: 5 minutes
  327.452897 Number of molecules: 895
  327.452903 Number of complexes: 895
  327.452908 Total reaction matches: 14627649
  327.452977 Avg timestep duration: 4.639216e-03, iterations remaining: 43000, Time left: 3 minutes
  327.452987 Estimated end time: 2026-10-14 13:44:55
  327.452993 -------------------------------------------------
  332.230743 -------------------
  332.230783 End iteration: 58000, simulation time: 5.800000e-02 seconds.
  332.231068 System time: 2026-10-14 13:41:41
  332.231082 Elapsed time: 5 minutes
  332.231088 Number of molecules: 895
  332.231093 Number of complexes: 895
  332.231098 Total reaction matches: 14825965
  332.231107 Avg timestep duration: 4.777694e-03, iterations remaining: 42000, Time left: 3 minutes
  332.231112 Estimated end time: 2026-10-14 13:45:01
  332.231117 -------------------------------------------------
  337.005453 -------------------
  337.005568 End iteration: 59000, simulation time: 5.900000e-02 seconds.
  337.005887 System time: 2026-10-14 13:41:45
  337.005900 Elapsed time: 5 minutes
  337.005907 Number of molecules: 895
  337.005912 Number of complexes: 895
  337.005917 Total reaction matches: 15005203
  337.005992 Avg timestep duration: 4.774306e-03, iterations remaining: 41000, Time left: 3 minutes
  337.006000 Estimated end time: 2026-10-14 13:45:01
  337.006005 -------------------------------------------------
  341.628216 -------------------
  341.628243 End iteration: 60000, simulation time: 6.000000e-02 seconds.
  341.628325 System time: 2026-10-14 13:41:50
  341.628336 Elapsed time: 5 minutes
  341.628341 Number of molecules: 895
  341.628347 Number of complexes: 895
  341.628352 Total reaction matches: 15181068
  341.628361 Avg timestep duration: 4.622161e-03, iterations remaining: 40000, Time left: 3 minutes
  341.628367 Estimated end time: 2026-10-14 13:44:55
  341.628372 -------------------------------------------------
  346.163877 -------------------
  346.165164 End iteration: 61000, simulation time: 6.100000e-02 seconds.
  346.165475 System time: 2026-10-14 13:41:55
  346.165490 Elapsed time: 5 minutes
  346.165500 Number of molecules: 895
  346.165506 Number of complexes: 895
  346.165511 Total reaction matches: 15276881
  346.165564 Avg timestep duration: 4.535483e-03, iterations remaining: 39000, Time left: 2 minutes
  346.165573 Estimated end time: 2026-10-14 13:44:51
  346.165578 -------------------------------------------------
  350.750374 -------------------
  350.750479 End iteration: 62000, simulation time: 6.200000e-02 seconds.
  350.750570 System time: 2026-10-14 13:41:59
  350.750585 Elapsed time: 5 minutes
  350.750591 Number of molecules: 895
  350.750597 Number of complexes: 895
  350.750602 Total reaction matches: 15471661
  350.750612 Avg timestep duration: 4.584705e-03, iterations remaining: 38000, Time left: 2 minutes
  350.750618 Estimated end time: 2026-10-14 13:44:53
  350.750623 -------------------------------------------------
  355.383533 -------------------
  355.383731 End iteration: 63000, simulation time: 6.300000e-02 seconds.
  355.384141 System time: 2026-10-14 13:42:04
  355.384157 Elapsed time: 5 minutes
  355.384164 Number of molecules: 895
  355.384169 Number of complexes: 895
  355.384174 Total reaction matches: 15749301
  355.384214 Avg timestep duration: 4.632891e-03, iterations remaining: 37000, Time left: 2 minutes
  355.384222 Estimated end time: 2026-10-14 13:44:55
  355.384227 -------------------------------------------------
  360.057552 -------------------
  360.057577 End iteration: 64000, simulation time: 6.400000e-02 seconds.
  360.057639 System time: 2026-10-14 13:42:08
  360.057650 Elapsed time: 6 minutes
  360.057656 Number of molecules: 895
  360.057661 Number of complexes: 895
  360.057665 Total reaction matches: 15972348
  360.057678 Avg timestep duration: 4.673269e-03, iterations remaining: 36000, Time left: 2 minutes
  360.057683 Estimated end time: 2026-10-14 13:44:57
  360.057688 -------------------------------------------------
  364.732108 -------------------
  364.732340 End iteration: 65000, simulation time: 6.500000e-02 seconds.
  364.732638 System time: 2026-10-14 13:42:13
  364.732654 Elapsed time: 6 minutes
  364.732661 Number of molecules: 895
  364.732667 Number of complexes: 895
  364.732672 Total reaction matches: 16314281
  364.732750 Avg timestep duration: 4.674384e-03, iterations remaining: 35000, Time left: 2 minutes
  364.732758 Estimated end time: 2026-10-14 13:44:57
  364.732790 -------------------------------------------------
  369.419883 -------------------
  369.419907 End iteration: 66000, simulation time: 6.600000e-02 seconds.
  369.419978 System time: 2026-10-14 13:42:18
  369.419987 Elapsed time: 6 minutes
  369.419992 Number of molecules: 895
  369.419996 Number of complexes: 895
  369.420000 Total reaction matches: 16640761
  369.420009 Avg timestep duration: 4.687069e-03, iterations remaining: 34000, Time left: 2 minutes
  369.420013 Estimated end time: 2026-10-14 13:44:57
  369.420017 -------------------------------------------------
  374.105561 -------------------
  374.105751 End iteration: 67000, simulation time: 6.700000e-02 seconds.
  374.106036 System time: 2026-10-14 13:42:22
  374.106049 Elapsed time: 6 minutes
  374.106055 Number of molecules: 895
  374.106064 Number of complexes: 895
  374.106069 Total reaction matches: 17114974
  374.106122 Avg timestep duration: 4.685508e-03, iterations remaining: 33000, Time left: 2 minutes
  374.106139 Estimated end time: 2026-10-14 13:44:57
  374.106146 -------------------------------------------------
  378.832553 -------------------
  378.832599 End iteration: 68000, simulation time: 6.800000e-02 seconds.
  378.832637 System time: 2026-10-14 13:42:27
  378.832646 Elapsed time: 6 minutes
  378.832652 Number of molecules: 895
  378.832656 Number of complexes: 895
  378.832661 Total reaction matches: 17537196
  378.832671 Avg timestep duration: 4.726349e-03, iterations remaining: 32000, Time left: 2 minutes
  378.832677 Estimated end time: 2026-10-14 13:44:58
  378.832682 -------------------------------------------------
  383.571557 -------------------
  383.571592 End iteration: 69000, simulation time: 6.900000e-02 seconds.
  383.571948 System time: 2026-10-14 13:42:32
  383.571960 Elapsed time: 6 minutes
  383.571966 Number of molecules: 895
  383.571972 Number of complexes: 895
  383.571977 Total reaction matches: 17848778
  383.571986 Avg timestep duration: 4.738821e-03, iterations remaining: 31000, Time left: 2 minutes
  383.571992 Estimated end time: 2026-10-14 13:44:59
  383.571997 -------------------------------------------------
  388.241517 -------------------
  388.241672 End iteration: 70000, simulation time: 7.000000e-02 seconds.
  388.242027 System time: 2026-10-14 13:42:37
  388.242046 Elapsed time: 6 minutes
  388.242053 Number of molecules: 895
  388.242059 Number of complexes: 895
  388.242065 Total reaction matches: 18052373
  388.242154 Avg timestep duration: 4.669482e-03, iterations remaining: 30000, Time left: 2 minutes
  388.242191 Estimated end time: 2026-10-14 13:44:57
  388.242207 -------------------------------------------------
  393.005067 -------------------
  393.005121 End iteration: 71000, simulation time: 7.100000e-02 seconds.
  393.005456 System time: 2026-10-14 13:42:41
  393.005468 Elapsed time: 6 minutes
  393.005474 Number of molecules: 895
  393.005479 Number of complexes: 895
  393.005484 Total reaction matches: 18328102
  393.005494 Avg timestep duration: 4.762800e-03, iterations remaining: 29000, Time left: 2 minutes
  393.005500 Estimated end time: 2026-10-14 13:44:59
  393.005505 -------------------------------------------------
  397.669899 -------------------
  397.670223 End iteration: 72000, simulation time: 7.200000e-02 seconds.
  397.670520 System time: 2026-10-14 13:42:46
  397.670531 Elapsed time: 6 minutes
  397.670536 Number of molecules: 895
  397.670542 Number of complexes: 895
  397.670546 Total reaction matches: 18562662
  397.670663 Avg timestep duration: 4.664320e-03, iterations remaining: 28000, Time left: 2 minutes
  397.670672 Estimated end time: 2026-10-14 13:44:57
  397.670677 -------------------------------------------------
  402.316009 -------------------
  402.316042 End iteration: 73000, simulation time: 7.300000e-02 seconds.
  402.316391 System time: 2026-10-14 13:42:51
  402.316403 Elapsed time: 6 minutes
  402.316408 Number of molecules: 895
  402.316413 Number of complexes: 895
  402.316418 Total reaction matches: 18773925
  402.316427 Avg timestep duration: 4.645268e-03, iterations remaining: 27000, Time left: 2 minutes
  402.316635 Estimated end time: 2026-10-14 13:44:56
  402.316641 -------------------------------------------------
  406.963640 -------------------
  406.963875 End iteration: 74000, simulation time: 7.400000e-02 seconds.
  406.964196 System time: 2026-10-14 13:42:55
  406.964209 Elapsed time: 6 minutes
  406.964214 Number of molecules: 895
  406.964219 Number of complexes: 895
  406.964224 Total reaction matches: 18960311
  406.964261 Avg timestep duration: 4.647164e-03, iterations remaining: 26000, Time left: 2 minutes
  406.964269 Estimated end time: 2026-10-14 13:44:56
  406.964274 -------------------------------------------------
  411.691370 -------------------
  411.691391 End iteration: 75000, simulation time: 7.500000e-02 seconds.
  411.691630 System time: 2026-10-14 13:43:00
  411.691640 Elapsed time: 6 minutes
  411.691646 Number of molecules: 895
  411.691651 Number of complexes: 895
  411.691656 Total reaction matches: 19255614
  411.691661 Avg timestep duration: 4.727026e-03, iterations remaining: 25000, Time left: 1 minutes
  411.691672 Estimated end time: 2026-10-14 13:44:58
  411.691677 -------------------------------------------------
  416.356170 -------------------
  416.357579 End iteration: 76000, simulation time: 7.600000e-02 seconds.
  416.357910 System time: 2026-10-14 13:43:05
  416.357921 Elapsed time: 6 minutes
  416.357926 Number of molecules: 895
  416.357931 Number of complexes: 895
  416.357936 Total reaction matches: 19521111
  416.357991 Avg timestep duration: 4.664426e-03, iterations remaining: 24000, Time left: 1 minutes
  416.358008 Estimated end time: 2026-10-14 13:44:57
  416.358018 -------------------------------------------------
  420.964988 -------------------
  420.965037 End iteration: 77000, simulation time: 7.700000e-02 seconds.
  420.965372 System time: 2026-10-14 13:43:09
  420.965387 Elapsed time: 7 minutes
  420.965394 Number of molecules: 895
  420.965400 Number of complexes: 895
  420.965405 Total reaction matches: 19701396
  420.965416 Avg timestep duration: 4.606904e-03, iterations remaining: 23000, Time left: 1 minutes
  420.965422 Estimated end time: 2026-10-14 13:44:55
  420.965427 -------------------------------------------------
  425.772876 -------------------
  425.773165 End iteration: 78000, simulation time: 7.800000e-02 seconds.
  425.773542 System time: 2026-10-14 13:43:14
  425.773555 Elapsed time: 7 minutes
  425.773561 Number of molecules: 895
  425.773565 Number of complexes: 895
  425.773570 Total reaction matches: 19982264
  425.773685 Avg timestep duration: 4.807376e-03, iterations remaining: 22000, Time left: 1 minutes
  425.773694 Estimated end time: 2026-10-14 13:45:00
  425.773699 -------------------------------------------------
  430.548635 -------------------
  430.548659 End iteration: 79000, simulation time: 7.900000e-02 seconds.
  430.548927 System time: 2026-10-14 13:43:19
  430.548948 Elapsed time: 7 minutes
  430.548956 Number of molecules: 895
  430.548963 Number of complexes: 895
  430.548969 Total reaction matches: 20218250
  430.548976 Avg timestep duration: 4.774867e-03, iterations remaining: 21000, Time left: 1 minutes
  430.548983 Estimated end time: 2026-10-14 13:44:59
  430.548990 -------------------------------------------------
  435.278682 -------------------
  435.278805 End iteration: 80000, simulation time: 8.000000e-02 seconds.
  435.279116 System time: 2026-10-14 13:43:24
  435.279128 Elapsed time: 7 minutes
  435.279134 Number of molecules: 895
  435.279138 Number of complexes: 895
  435.279143 Total reaction matches: 20608382
  435.279196 Avg timestep duration: 4.729670e-03, iterations remaining: 20000, Time left: 1 minutes
  435.279204 Estimated end time: 2026-10-14 13:44:58
  435.279209 -------------------------------------------------
  439.976527 -------------------
  439.976561 End iteration: 81000, simulation time: 8.100000e-02 seconds.
  439.976881 System time: 2026-10-14 13:43:28
  439.976893 Elapsed time: 7 minutes
  439.977083 Number of molecules: 895
  439.977096 Number of complexes: 895
  439.977102 Total reaction matches: 21032023
  439.977118 Avg timestep duration: 4.697226e-03, iterations remaining: 19000, Time left: 1 minutes
  439.977124 Estimated end time: 2026-10-14 13:44:58
  439.977128 -------------------------------------------------
  444.752415 -------------------
  444.752741 End iteration: 82000, simulation time: 8.200000e-02 seconds.
  444.753047 System time: 2026-10-14 13:43:33
  444.753058 Elapsed time: 7 minutes
  444.753064 Number of molecules: 895
  444.753069 Number of complexes: 895
  444.753074 Total reaction matches: 21356442
  444.753182 Avg timestep duration: 4.775434e-03, iterations remaining: 18000, Time left: 1 minutes
  444.753190 Estimated end time: 2026-10-14 13:44:59
  444.753195 -------------------------------------------------
  449.423975 -------------------
  449.424026 End iteration: 83000, simulation time: 8.300000e-02 seconds.
  449.424176 System time: 2026-10-14 13:43:38
  449.424186 Elapsed time: 7 minutes
  449.424193 Number of molecules: 895
  449.424199 Number of complexes: 895
  449.424204 Total reaction matches: 21588656
  449.424213 Avg timestep duration: 4.670687e-03, iterations remaining: 17000, Time left: 1 minutes
  449.424222 Estimated end time: 2026-10-14 13:44:57
  449.424227 -------------------------------------------------
  454.543373 -------------------
  454.543596 End iteration: 84000, simulation time: 8.400000e-02 seconds.
  454.543955 System time: 2026-10-14 13:43:43
  454.543977 Elapsed time: 7 minutes
  454.543986 Number of molecules: 895
  454.543993 Number of complexes: 895
  454.544000 Total reaction matches: 21787031
  454.544065 Avg timestep duration: 5.119102e-03, iterations remaining: 16000, Time left: 1 minutes
  454.544081 Estimated end time: 2026-10-14 13:45:05
  454.544095 -------------------------------------------------
  460.431497 -------------------
  460.431542 End iteration: 85000, simulation time: 8.500000e-02 seconds.
  460.431584 System time: 2026-10-14 13:43:49
  460.431595 Elapsed time: 7 minutes
  460.431601 Number of molecules: 895
  460.431607 Number of complexes: 895
  460.431613 Total reaction matches: 22019538
  460.431623 Avg timestep duration: 5.887302e-03, iterations remaining: 15000, Time left: 1 minutes
  460.431629 Estimated end time: 2026-10-14 13:45:17
  460.431634 -------------------------------------------------
  465.079477 -------------------
  465.079740 End iteration: 86000, simulation time: 8.600000e-02 seconds.
  465.080182 System time: 2026-10-14 13:43:53
  465.080196 Elapsed time: 7 minutes
  465.080202 Number of molecules: 895
  465.080211 Number of complexes: 895
  465.080222 Total reaction matches: 22265585
  465.080248 Avg timestep duration: 4.647781e-03, iterations remaining: 14000, Time left: 1 minutes
  465.080256 Estimated end time: 2026-10-14 13:44:59
  465.080261 -------------------------------------------------
  469.735264 -------------------
  469.735313 End iteration: 87000, simulation time: 8.700000e-02 seconds.
  469.735334 System time: 2026-10-14 13:43:58
  469.735344 Elapsed time: 7 minutes
  469.735350 Number of molecules: 895
  469.735356 Number of complexes: 895
  469.735361 Total reaction matches: 22545244
  469.735369 Avg timestep duration: 4.654910e-03, iterations remaining: 13000, Time left: 1 minutes
  469.735376 Estimated end time: 2026-10-14 13:44:59
  469.735380 -------------------------------------------------
  474.391161 -------------------
  474.391211 End iteration: 88000, simulation time: 8.800000e-02 seconds.
  474.391249 System time: 2026-10-14 13:44:03
  474.391259 Elapsed time: 7 minutes
  474.391265 Number of molecules: 895
  474.391270 Number of complexes: 895
  474.391274 Total reaction matches: 22748595
  474.391284 Avg timestep duration: 4.655696e-03, iterations remaining: 12000, Time left: 0 minutes
  474.391290 Estimated end time: 2026-10-14 13:44:59
  474.391295 -------------------------------------------------
  479.139626 -------------------
  479.139873 End iteration: 89000, simulation time: 8.900000e-02 seconds.
  479.140342 System time: 2026-10-14 13:44:08
  479.140359 Elapsed time: 7 minutes
  479.140375 Number of molecules: 895
  479.140381 Number of complexes: 895
  479.140387 Total reaction matches: 22952628
  479.140419 Avg timestep duration: 4.748240e-03, iterations remaining: 11000, Time left: 0 minutes
  479.140427 Estimated end time: 2026-10-14 13:45:00
  479.140433 -------------------------------------------------
  484.001773 -------------------
  484.001834 End iteration: 90000, simulation time: 9.000000e-02 seconds.
  484.001932 System time: 2026-10-14 13:44:12
  484.001944 Elapsed time: 8 minutes
  484.001950 Number of molecules: 895
  484.001957 Number of complexes: 895
  484.001962 Total reaction matches: 23271479
  484.001972 Avg timestep duration: 4.861229e-03, iterations remaining: 10000, Time left: 0 minutes
  484.001978 Estimated end time: 2026-10-14 13:45:01
  484.001982 -------------------------------------------------
  488.728095 -------------------
  488.728130 End iteration: 91000, simulation time: 9.100000e-02 seconds.
  488.728190 System time: 2026-10-14 13:44:17
  488.728201 Elapsed time: 8 minutes
  488.728207 Number of molecules: 895
  488.728212 Number of complexes: 895
  488.728218 Total reaction matches: 23594198
  488.728231 Avg timestep duration: 4.726051e-03, iterations remaining: 9000, Time left: 0 minutes
  488.728237 Estimated end time: 2026-10-14 13:45:00
  488.728242 -------------------------------------------------
  493.398155 -------------------
  493.398424 End iteration: 92000, simulation time: 9.200000e-02 seconds.
  493.398832 System time: 2026-10-14 13:44:22
  493.398845 Elapsed time: 8 minutes
  493.398852 Number of molecules: 895
  493.398859 Number of complexes: 895
  493.398864 Total reaction matches: 23876563
  493.398932 Avg timestep duration: 4.669842e-03, iterations remaining: 8000, Time left: 0 minutes
  493.398941 Estimated end time: 2026-10-14 13:44:59
  493.398947 -------------------------------------------------
  498.089164 -------------------
  498.089202 End iteration: 93000, simulation time: 9.300000e-02 seconds.
  498.089251 System time: 2026-10-14 13:44:26
  498.089262 Elapsed time: 8 minutes
  498.089268 Number of molecules: 895
  498.089274 Number of complexes: 895
  498.089280 Total reaction matches: 24186479
  498.089290 Avg timestep duration: 4.690145e-03, iterations remaining: 7000, Time left: 0 minutes
  498.089296 Estimated end time: 2026-10-14 13:44:59
  498.089301 -------------------------------------------------
  502.888494 -------------------
  502.888538 End iteration: 94000, simulation time: 9.400000e-02 seconds.
  502.888576 System time: 2026-10-14 13:44:31
  502.888585 Elapsed time: 8 minutes
  502.888591 Number of molecules: 895
  502.888596 Number of complexes: 895
  502.888600 Total reaction matches: 24439455
  502.888610 Avg timestep duration: 4.799116e-03, iterations remaining: 6000, Time left: 0 minutes
  502.888616 Estimated end time: 2026-10-14 13:45:00
  502.888621 -------------------------------------------------
  507.531705 -------------------
  507.531950 End iteration: 95000, simulation time: 9.500000e-02 seconds.
  507.532430 System time: 2026-10-14 13:44:36
  507.532448 Elapsed time: 8 minutes
  507.532461 Number of molecules: 895
  507.532468 Number of complexes: 895
  507.532474 Total reaction matches: 24629194
  507.532488 Avg timestep duration: 4.642997e-03, iterations remaining: 5000, Time left: 0 minutes
  507.532495 Estimated end time: 2026-10-14 13:44:59
  507.532500 -------------------------------------------------
  512.482637 -------------------
  512.482692 End iteration: 96000, simulation time: 9.600000e-02 seconds.
  512.482766 System time: 2026-10-14 13:44:41
  512.482776 Elapsed time: 8 minutes
  512.482782 Number of molecules: 895
  512.482787 Number of complexes: 895
  512.482791 Total reaction matches: 24810079
  512.482801 Avg timestep duration: 4.949994e-03, iterations remaining: 4000, Time left: 0 minutes
  512.482806 Estimated end time: 2026-10-14 13:45:01
  512.482811 -------------------------------------------------
  517.111915 -------------------
  517.112197 End iteration: 97000, simulation time: 9.700000e-02 seconds.
  517.112644 System time: 2026-10-14 13:44:45
  517.112659 Elapsed time: 8 minutes
  517.112670 Number of molecules: 895
  517.112676 Number of complexes: 895
  517.112681 Total reaction matches: 24971714
  517.112693 Avg timestep duration: 4.629018e-03, iterations remaining: 3000, Time left: 0 minutes
  517.112702 Estimated end time: 2026-10-14 13:44:59
  517.112708 -------------------------------------------------
  521.802106 -------------------
  521.802142 End iteration: 98000, simulation time: 9.800000e-02 seconds.
  521.802202 System time: 2026-10-14 13:44:50
  521.802217 Elapsed time: 8 minutes
  521.802224 Number of molecules: 895
  521.802230 Number of complexes: 895
  521.802236 Total reaction matches: 25278618
  521.802248 Avg timestep duration: 4.689328e-03, iterations remaining: 2000, Time left: 0 minutes
  521.802255 Estimated end time: 2026-10-14 13:45:00
  521.802260 -------------------------------------------------
  526.489611 -------------------
  526.489659 End iteration: 99000, simulation time: 9.900000e-02 seconds.
  526.489700 System time: 2026-10-14 13:44:55
  526.489710 Elapsed time: 8 minutes
  526.489716 Number of molecules: 895
  526.489721 Number of complexes: 895
  526.489726 Total reaction matches: 25468519
  526.489735 Avg timestep duration: 4.687284e-03, iterations remaining: 1000, Time left: 0 minutes
  526.489740 Estimated end time: 2026-10-14 13:45:00
  526.489745 -------------------------------------------------
  531.190692 -------------------------------------------------
  531.190848 End simulation
  531.190865 End date: 2026-10-14 13:45:00
  531.190872 	Wall Time: 5.311853e+02 seconds
  531.191757 [exit status 0]
