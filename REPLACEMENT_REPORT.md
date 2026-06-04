# Repository Replacement Report

Generated: 2026-06-04 14:31:53 local time

## Scope

- Target Git repository: `/Users/yueying/Workspace/NERDSS`
- Source directory: `/Users/yueying/Downloads/nerdss_development-masterv2`
- Replacement command run:

```bash
rsync -a --checksum --delete --exclude='.git/' /Users/yueying/Downloads/nerdss_development-masterv2/ /Users/yueying/Workspace/NERDSS/
```

## Summary

- Added files planned by dry run: 443
- Deleted files planned by dry run: 641
- Content-modified files planned by dry run: 16
- Metadata-only file differences planned by dry run: 1544
- Unchanged source files identified by dry run: 1
- Directory metadata/delete entries in dry run: 155

## Skipped Files And Assumptions

- Skipped `.git/` intentionally, preserving the current repository history, remotes, branches, and index metadata.
- Did not skip `.gitignore`; it was replaced because the source intentionally contains a different `.gitignore`.
- Did not skip other source hidden paths. The source includes `.cache/clangd/...`, so those files were copied into the working tree.
- `REPLACEMENT_REPORT.md` was generated after the replacement as a local documentation artifact and is not part of the source directory.
- No commit was created.

## Git Verification Snapshot

This snapshot was captured after replacement and after creating this report file.

### git status --porcelain=v1 -uall

```
 M .gitignore
 M EXEs/nerdss.cpp
 M Makefile
 M Makefile_old
 M include/boundary_conditions/reflect_functions.hpp
 M include/classes/class_Molecule_Complex.hpp
 D sample_inputs/VALIDATE_SUITE/hetTrimer/PDB/0.pdb
 D sample_inputs/VALIDATE_SUITE/hetTrimer/RESTARTS/restart120000.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/RESTARTS/restart180000.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/RESTARTS/restart60000.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/assoc_dissoc_time.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/bound_pair_time.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/copy_numbers_time.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/event_counters_time.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/histogram_complexes_time.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/initial_crds.xyz
 D sample_inputs/VALIDATE_SUITE/hetTrimer/mono_dimer_time.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/observables_time.dat
 D sample_inputs/VALIDATE_SUITE/hetTrimer/system.psf
 D sample_inputs/VALIDATE_SUITE/hetTrimer/trajectory.xyz
 D sample_inputs/VALIDATE_SUITE/hetTrimer/transition_matrix_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/0.pdb
 D sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/1000000.pdb
 D sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/1999999.pdb
 D sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/599999.pdb
 D sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/9999.pdb
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1000000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart120000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1200000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1400000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1600000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart180000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1800000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart2000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart200000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart240000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart3000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart300000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart360000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart4000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart400000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart420000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart480000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart5000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart540000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart6000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart60000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart600000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart7000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart8000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart800000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart9000.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/assoc_dissoc_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/bound_pair_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/copy_numbers_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/event_counters_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/final_coords.xyz
 D sample_inputs/VALIDATE_SUITE/homoTrimer/histogram_complexes_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/initial_crds.xyz
 D sample_inputs/VALIDATE_SUITE/homoTrimer/mono_dimer_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/observables_time.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/restart.dat
 D sample_inputs/VALIDATE_SUITE/homoTrimer/system.psf
 D sample_inputs/VALIDATE_SUITE/homoTrimer/trajectory.xyz
 D sample_inputs/VALIDATE_SUITE/homoTrimer/transition_matrix_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/0.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/100000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1000000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1100000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1200000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1300000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1400000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1500000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1600000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1700000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1800000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/1900000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/200000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2000000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2100000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2200000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2300000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2400000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2500000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2600000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2700000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2800000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/2900000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/300000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3000000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3100000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3200000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3300000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3400000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3500000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3600000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3700000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3800000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/3900000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/400000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4000000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4100000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4200000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4300000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4400000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4500000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4600000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4700000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4800000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4900000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/4999999.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/500000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/600000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/700000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/800000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/PDB/900000.pdb
 D sample_inputs/VALIDATE_SUITE/trimer/assoc_dissoc_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/bound_pair_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/copy_numbers_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/event_counters_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/final_coords.xyz
 D sample_inputs/VALIDATE_SUITE/trimer/histogram_complexes_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/initial_crds.xyz
 D sample_inputs/VALIDATE_SUITE/trimer/mono_dimer_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/observables_time.dat
 D sample_inputs/VALIDATE_SUITE/trimer/restart.dat
 D sample_inputs/VALIDATE_SUITE/trimer/system.psf
 D sample_inputs/VALIDATE_SUITE/trimer/trajectory.xyz
 D sample_inputs/VALIDATE_SUITE/trimer/transition_matrix_time.dat
 M src/classes/class_Molecule_Complex.cpp
 M src/classes/class_Parameters.cpp
 M src/classes/class_Vector.cpp
 M src/io/read_restart.cpp
 M src/io/write_restart.cpp
 M src/reactions/associate_box.cpp
 M src/reactions/check_for_unimolecular_reactions_population.cpp
 M src/reactions/determine_1D_bimolecular_reaction_probability.cpp
 M src/reactions/pirr_pfree_ratio_psF_1D.cpp
 M src/trajectory_functions/calculate_update_position_interface.cpp
 M src/trajectory_functions/create_complex_propagation_vectors.cpp
 D system.psf
 M text_search.bsh
?? .cache/clangd/index/2D_reaction_table_functions.hpp.E01AB59A14764AC9.idx
?? .cache/clangd/index/DDpirr_pfree_ratio_ps.cpp.596F1BDD8B49E7DF.idx
?? .cache/clangd/index/DDpirr_pfree_ratio_ps.cpp.C4682E31B2684260.idx
?? .cache/clangd/index/Faddeeva.cpp.599BF93AF027009F.idx
?? .cache/clangd/index/Faddeeva.cpp.FEB6F97561210959.idx
?? .cache/clangd/index/Faddeeva.hpp.31C1D40B1BF1B9A0.idx
?? .cache/clangd/index/angleSignIsCorrect.cpp.346BCB05078F6C0C.idx
?? .cache/clangd/index/angleSignIsCorrect.cpp.9AD06B6E47C4053B.idx
?? .cache/clangd/index/areInVicinity.cpp.4B463C5EFAF9A5CA.idx
?? .cache/clangd/index/areInVicinity.cpp.4D8CC006E6E845BE.idx
?? .cache/clangd/index/areParallel.cpp.FBDA2F60160EC332.idx
?? .cache/clangd/index/areParallel.cpp.FCD65FAE99F26C0C.idx
?? .cache/clangd/index/areSameAngle.cpp.7BB77C900F78C285.idx
?? .cache/clangd/index/areSameAngle.cpp.D1C9C1D1812B119A.idx
?? .cache/clangd/index/areSameExceptState.cpp.688A3B589E0602E7.idx
?? .cache/clangd/index/areSameExceptState.cpp.D3D9FFD38DD5AEE7.idx
?? .cache/clangd/index/associate.cpp.88C8B5BDAE14774A.idx
?? .cache/clangd/index/associate.cpp.B2EB00E69DB9CAE7.idx
?? .cache/clangd/index/associate_ImplicitLipid.cpp.602E5A2F2AC34B7C.idx
?? .cache/clangd/index/associate_ImplicitLipid.cpp.994C548EC5BAEDF0.idx
?? .cache/clangd/index/associate_ImplicitLipid_box.cpp.1E1E65CA086E3A44.idx
?? .cache/clangd/index/associate_ImplicitLipid_box.cpp.AE174A59E0A21AB6.idx
?? .cache/clangd/index/associate_ImplicitLipid_sphere.cpp.48EEC1BA49BC613F.idx
?? .cache/clangd/index/associate_ImplicitLipid_sphere.cpp.F1D082D6F78F82E8.idx
?? .cache/clangd/index/associate_box.cpp.7C9F70CFE14BFF99.idx
?? .cache/clangd/index/associate_box.cpp.80218D6313AA983B.idx
?? .cache/clangd/index/associate_sphere.cpp.A1E17496A6B669BB.idx
?? .cache/clangd/index/associate_sphere.cpp.F0A5B19AAFC4B23E.idx
?? .cache/clangd/index/association.hpp.13A9340F987A3FAF.idx
?? .cache/clangd/index/bimolecular_reactions.hpp.A24B66D3C05A8948.idx
?? .cache/clangd/index/break_interaction.cpp.5212335BC8D27174.idx
?? .cache/clangd/index/break_interaction.cpp.9694D553B30D7489.idx
?? .cache/clangd/index/break_interaction_implicitlipid.cpp.4404381750380388.idx
?? .cache/clangd/index/break_interaction_implicitlipid.cpp.4D955D83E82D7C1B.idx
?? .cache/clangd/index/calc_angular_displacement.cpp.327AEAAA7C0A1090.idx
?? .cache/clangd/index/calc_angular_displacement.cpp.AEFA51D95ABB4532.idx
?? .cache/clangd/index/calc_one_angular_displacement.cpp.92AB01A168D1960B.idx
?? .cache/clangd/index/calc_one_angular_displacement.cpp.9B158482797B861E.idx
?? .cache/clangd/index/calc_pirr.cpp.3AF719EBB17BB3A4.idx
?? .cache/clangd/index/calc_pirr.cpp.C282BA08E48EEC99.idx
?? .cache/clangd/index/calculate_omega.cpp.26CA1EF283935390.idx
?? .cache/clangd/index/calculate_omega.cpp.720D33FAD6633F7C.idx
?? .cache/clangd/index/calculate_phi.cpp.BDADAD0C353D6443.idx
?? .cache/clangd/index/calculate_phi.cpp.E41D6C7D710850BC.idx
?? .cache/clangd/index/calculate_update_position_interface.cpp.0E1D62689CE6FDA1.idx
?? .cache/clangd/index/calculate_update_position_interface.cpp.3D1DD3C222C67D9B.idx
?? .cache/clangd/index/check_bases.cpp.705FEED0F752FED5.idx
?? .cache/clangd/index/check_bases.cpp.92B255123CB6A911.idx
?? .cache/clangd/index/check_bimolecular_reactions.cpp.2BAFDCC67DD415A8.idx
?? .cache/clangd/index/check_bimolecular_reactions.cpp.3D7D5567B8E471DB.idx
?? .cache/clangd/index/check_compartment_reaction.cpp.82A072B6658F1087.idx
?? .cache/clangd/index/check_compartment_reaction.cpp.DE4FCA58E2B79A6C.idx
?? .cache/clangd/index/check_dissociation.cpp.087600F372E7F5E3.idx
?? .cache/clangd/index/check_dissociation.cpp.52F4CD98B1355803.idx
?? .cache/clangd/index/check_dissociation_implicitlipid.cpp.7194057AA6673E4B.idx
?? .cache/clangd/index/check_dissociation_implicitlipid.cpp.9FAFF3FDFE4C886A.idx
?? .cache/clangd/index/check_for_state_change.cpp.0D2C35DEBE67D58A.idx
?? .cache/clangd/index/check_for_state_change.cpp.EC7CA378EA07A354.idx
?? .cache/clangd/index/check_for_structure_overlap.cpp.0EE5B4861890BCDC.idx
?? .cache/clangd/index/check_for_structure_overlap.cpp.605A7A3DE1FD173D.idx
?? .cache/clangd/index/check_for_structure_overlap_system.cpp.21C718F8CAA45044.idx
?? .cache/clangd/index/check_for_structure_overlap_system.cpp.ED27AC02F58D3B56.idx
?? .cache/clangd/index/check_for_unimolecular_reactions.cpp.1D95BB89B5AEF292.idx
?? .cache/clangd/index/check_for_unimolecular_reactions.cpp.E61F3EAE361A703A.idx
?? .cache/clangd/index/check_for_unimolecular_reactions_population.cpp.064A3C05887F547D.idx
?? .cache/clangd/index/check_for_unimolecular_reactions_population.cpp.90B110B8E372D347.idx
?? .cache/clangd/index/check_for_unimolstatechange_reactions.cpp.159A531B8BE6B6F4.idx
?? .cache/clangd/index/check_for_unimolstatechange_reactions.cpp.4B2B4CE87FDE1D7D.idx
?? .cache/clangd/index/check_for_valid_states.cpp.1C38B1412F4FE96D.idx
?? .cache/clangd/index/check_for_valid_states.cpp.37A9ACAD134AEFF5.idx
?? .cache/clangd/index/check_for_zeroth_order_creation.cpp.6F525395F7241967.idx
?? .cache/clangd/index/check_for_zeroth_order_creation.cpp.DBBFD5CFD238D079.idx
?? .cache/clangd/index/check_if_spans.cpp.83E8E1A2DB912F90.idx
?? .cache/clangd/index/check_if_spans.cpp.B25B0590DCB34027.idx
?? .cache/clangd/index/check_if_spans_box.cpp.5D284F79C405A37E.idx
?? .cache/clangd/index/check_if_spans_box.cpp.C5D8664F46DE09B0.idx
?? .cache/clangd/index/check_if_spans_sphere.cpp.E135304E3588AD48.idx
?? .cache/clangd/index/check_if_spans_sphere.cpp.E3CBCAB4A1C2E05C.idx
?? .cache/clangd/index/check_implicit_reactions.cpp.77FF45D85157C8AB.idx
?? .cache/clangd/index/check_implicit_reactions.cpp.C7A4AB4A551E5101.idx
?? .cache/clangd/index/class_Cluster.cpp.BD160A36A90F5798.idx
?? .cache/clangd/index/class_Cluster.cpp.D1477C0585549194.idx
?? .cache/clangd/index/class_Cluster.hpp.E89C938C2F55233D.idx
?? .cache/clangd/index/class_Coord.cpp.199DC27C605D7225.idx
?? .cache/clangd/index/class_Coord.cpp.A28391A9D1673528.idx
?? .cache/clangd/index/class_Coord.hpp.3AB0AB56D94C613E.idx
?? .cache/clangd/index/class_MDTimer.cpp.29402414A915C06F.idx
?? .cache/clangd/index/class_MDTimer.cpp.E55D4F42BFFEE3F6.idx
?? .cache/clangd/index/class_MDTimer.hpp.6E2DDE8C259E6173.idx
?? .cache/clangd/index/class_Membrane.hpp.D1CB8CD4E39E5BE1.idx
?? .cache/clangd/index/class_MolTemplate.cpp.0047A1AC1C5A577C.idx
?? .cache/clangd/index/class_MolTemplate.cpp.3DBF9D6A5F5AC21E.idx
?? .cache/clangd/index/class_MolTemplate.hpp.D2862D54454C5B3D.idx
?? .cache/clangd/index/class_Molecule_Complex.cpp.06DC17412EC80C5B.idx
?? .cache/clangd/index/class_Molecule_Complex.cpp.DE1C3BB70572DEDF.idx
?? .cache/clangd/index/class_Molecule_Complex.hpp.A5431CAD3364BAD6.idx
?? .cache/clangd/index/class_Observable.cpp.3DD6007F75742128.idx
?? .cache/clangd/index/class_Observable.cpp.F9748ACB8339B771.idx
?? .cache/clangd/index/class_Observable.hpp.B06702238E12CEA8.idx
?? .cache/clangd/index/class_Parameters.cpp.6530F1C58E6F75D9.idx
?? .cache/clangd/index/class_Parameters.cpp.9E55625509C6FE93.idx
?? .cache/clangd/index/class_Parameters.hpp.824301A6ACE9B9C7.idx
?? .cache/clangd/index/class_Quat.cpp.E352F5EDC32944E7.idx
?? .cache/clangd/index/class_Quat.cpp.E8318091DA02AC57.idx
?? .cache/clangd/index/class_Quat.hpp.76D2E15A147C9793.idx
?? .cache/clangd/index/class_Rxns.cpp.66F6DFCFFC272689.idx
?? .cache/clangd/index/class_Rxns.cpp.CA479397F29C8CC0.idx
?? .cache/clangd/index/class_Rxns.hpp.37F41A5206ACC8CA.idx
?? .cache/clangd/index/class_SimulVolume.cpp.6346D705697E3012.idx
?? .cache/clangd/index/class_SimulVolume.cpp.FAA20F64BD9149F7.idx
?? .cache/clangd/index/class_SimulVolume.hpp.35E9A0BA21F5EA27.idx
?? .cache/clangd/index/class_Vector.cpp.01B952DD7BFA8EB3.idx
?? .cache/clangd/index/class_Vector.cpp.AFC2EE1F913CE571.idx
?? .cache/clangd/index/class_Vector.hpp.D2667BB5B7554218.idx
?? .cache/clangd/index/class_bngl_parser.hpp.9331DCCCB8184F52.idx
?? .cache/clangd/index/class_bngl_parser_functions.cpp.55F7A9D1332DC984.idx
?? .cache/clangd/index/class_bngl_parser_functions.cpp.576239292709E75E.idx
?? .cache/clangd/index/class_copyCounters.hpp.09CB36E1C86E15FA.idx
?? .cache/clangd/index/clear_reweight_vecs.cpp.442DD75210EA6839.idx
?? .cache/clangd/index/clear_reweight_vecs.cpp.C3D0747F4E7A31DC.idx
?? .cache/clangd/index/com_of_two_tmp_complexes.cpp.A5E45E0CA2588523.idx
?? .cache/clangd/index/com_of_two_tmp_complexes.cpp.CC0FDF914757913A.idx
?? .cache/clangd/index/conservedMags.cpp.C940BA855CBAA293.idx
?? .cache/clangd/index/conservedMags.cpp.E0DF51220127458E.idx
?? .cache/clangd/index/conservedRigid.cpp.235A441583D8610D.idx
?? .cache/clangd/index/conservedRigid.cpp.86323E869171B346.idx
?? .cache/clangd/index/constants.hpp.FDF13DEA60CE30B6.idx
?? .cache/clangd/index/create_DDMatrices.cpp.B38F0B34E290FACE.idx
?? .cache/clangd/index/create_DDMatrices.cpp.E8CAA1FA5A352C85.idx
?? .cache/clangd/index/create_arbitrary_vector.cpp.44CC616EED7BAB7D.idx
?? .cache/clangd/index/create_arbitrary_vector.cpp.4A85B44D1CA088D6.idx
?? .cache/clangd/index/create_complex_propagation_vectors.cpp.567D0926945E3DF7.idx
?? .cache/clangd/index/create_complex_propagation_vectors.cpp.767D759168397DF7.idx
?? .cache/clangd/index/create_complex_propagation_vectors_on_sphere.cpp.81DF5D883F10092D.idx
?? .cache/clangd/index/create_complex_propagation_vectors_on_sphere.cpp.C8F9B73F6D69E81B.idx
?? .cache/clangd/index/create_conjugate_reaction_itrs.cpp.0F50AC44BCE0F615.idx
?? .cache/clangd/index/create_conjugate_reaction_itrs.cpp.B6929A3863FC2A7C.idx
?? .cache/clangd/index/create_molecule_and_complex_from_rxn.cpp.AA2ED144BE580052.idx
?? .cache/clangd/index/create_molecule_and_complex_from_rxn.cpp.F88B63867BEA535F.idx
?? .cache/clangd/index/create_molecule_and_complex_from_transmission_rxn.cpp.6A453F4C7557C781.idx
?? .cache/clangd/index/create_molecule_and_complex_from_transmission_rxn.cpp.B4BE3EE0144F450A.idx
?? .cache/clangd/index/create_normMatrix.cpp.4B4823C0237F61D7.idx
?? .cache/clangd/index/create_normMatrix.cpp.98BF6B26BF50FF30.idx
?? .cache/clangd/index/create_pirMatrix.cpp.3B2EAF45E9D23109.idx
?? .cache/clangd/index/create_pirMatrix.cpp.9801504D1C08EB78.idx
?? .cache/clangd/index/create_survMatrix.cpp.88418B8D81960E00.idx
?? .cache/clangd/index/create_survMatrix.cpp.9434F5360EA93885.idx
?? .cache/clangd/index/determine_1D_bimolecular_reaction_probability.cpp.3124089C1A564B84.idx
?? .cache/clangd/index/determine_2D_bimolecular_reaction_probability.cpp.51ADD49D63731880.idx
?? .cache/clangd/index/determine_2D_bimolecular_reaction_probability.cpp.62809B98425D6EF3.idx
?? .cache/clangd/index/determine_2D_implicitlipid_reaction_probability.cpp.4DB82439422FF127.idx
?? .cache/clangd/index/determine_2D_implicitlipid_reaction_probability.cpp.D63E9DC35F63DAF1.idx
?? .cache/clangd/index/determine_3D_bimolecular_reaction_probability.cpp.55EDD4001C5898F4.idx
?? .cache/clangd/index/determine_3D_bimolecular_reaction_probability.cpp.E52AB0E0C63B2289.idx
?? .cache/clangd/index/determine_3D_implicitlipid_reaction_probability.cpp.93F32C41AB2FA0C0.idx
?? .cache/clangd/index/determine_3D_implicitlipid_reaction_probability.cpp.BD52A187A620821A.idx
?? .cache/clangd/index/determine_bound_iface_index.cpp.57893AD913278BBE.idx
?? .cache/clangd/index/determine_bound_iface_index.cpp.5D4FDF29AC0B5CAC.idx
?? .cache/clangd/index/determine_entering_compartment_probability.cpp.0584F27799043DCC.idx
?? .cache/clangd/index/determine_entering_compartment_probability.cpp.A6057F141E49B096.idx
?? .cache/clangd/index/determine_exiting_compartment_probability.cpp.0FCF25CCDEED9A14.idx
?? .cache/clangd/index/determine_exiting_compartment_probability.cpp.1A292F4787C47097.idx
?? .cache/clangd/index/determine_if_reaction_occurs.cpp.1166F7BD9ACCF1CA.idx
?? .cache/clangd/index/determine_if_reaction_occurs.cpp.BDDD131A11CA9C74.idx
?? .cache/clangd/index/determine_iface_indices.cpp.6598A92C6837D473.idx
?? .cache/clangd/index/determine_iface_indices.cpp.EF0D6C4D7BC35D4E.idx
?? .cache/clangd/index/determine_normal.cpp.7D409ABAAD08B489.idx
?? .cache/clangd/index/determine_normal.cpp.E5111DA7D4CD5CEB.idx
?? .cache/clangd/index/determine_parent_complex.cpp.2A89C8EF52F98AFD.idx
?? .cache/clangd/index/determine_parent_complex.cpp.494AE6A4CFD159DD.idx
?? .cache/clangd/index/determine_parent_complex_IL.cpp.2CB3A606FC30A567.idx
?? .cache/clangd/index/determine_parent_complex_IL.cpp.C5A36D94B6BEBD5F.idx
?? .cache/clangd/index/determine_rotation_angles.cpp.1F92E765FE04DD04.idx
?? .cache/clangd/index/determine_rotation_angles.cpp.200E7DE777341D93.idx
?? .cache/clangd/index/determine_shape_molecule.cpp.0A060223FBC79DFD.idx
?? .cache/clangd/index/determine_shape_molecule.cpp.1F7C4351A0C5AB7F.idx
?? .cache/clangd/index/display_all_MolTemplates.cpp.22E69D0CC546E283.idx
?? .cache/clangd/index/display_all_MolTemplates.cpp.AC2BA8CFCB08679F.idx
?? .cache/clangd/index/display_all_reactions.cpp.B45C56E1E0A39E93.idx
?? .cache/clangd/index/display_all_reactions.cpp.D690C19E94611240.idx
?? .cache/clangd/index/error_handling.hpp.00DB20F9513BFBB4.idx
?? .cache/clangd/index/evaluate_binding_within_complex.cpp.8F7EEBF22E084286.idx
?? .cache/clangd/index/evaluate_binding_within_complex.cpp.C39FACEBBB326A95.idx
?? .cache/clangd/index/eye_candy.cpp.9E7D79E279A61098.idx
?? .cache/clangd/index/eye_candy.cpp.A6F4FA8C8D902CD8.idx
?? .cache/clangd/index/find_reaction_rate_state.cpp.2783325BB3D9E0AE.idx
?? .cache/clangd/index/find_reaction_rate_state.cpp.DD748A03243AC975.idx
?? .cache/clangd/index/find_which_reaction.cpp.19A721BAEC7E453A.idx
?? .cache/clangd/index/find_which_reaction.cpp.42C0A8EE96B6859C.idx
?? .cache/clangd/index/find_which_state_change_reaction.cpp.6823B9A51CDD9CA2.idx
?? .cache/clangd/index/find_which_state_change_reaction.cpp.DAD941B7FA0A20E6.idx
?? .cache/clangd/index/functions_for_spherical_system.cpp.13BE19BE8FC5AE2C.idx
?? .cache/clangd/index/functions_for_spherical_system.cpp.BAF403A10FD19FF3.idx
?? .cache/clangd/index/functions_for_spherical_system.hpp.44C914A5BACE35CE.idx
?? .cache/clangd/index/functions_implicitlipid.cpp.801D27A088B90481.idx
?? .cache/clangd/index/functions_implicitlipid.cpp.95F8F6317607B0EA.idx
?? .cache/clangd/index/generate_coordinates.cpp.16B2A7A3181F9215.idx
?? .cache/clangd/index/generate_coordinates.cpp.D69DCC233D5F38E8.idx
?? .cache/clangd/index/generate_coordinates_for_restart.cpp.3D1771392B23CC5D.idx
?? .cache/clangd/index/generate_coordinates_for_restart.cpp.B16A2647FD8B39F6.idx
?? .cache/clangd/index/get_distance.cpp.13BC466368CA022D.idx
?? .cache/clangd/index/get_distance.cpp.1E4410EE36BCFAB3.idx
?? .cache/clangd/index/get_distance_to_surface.cpp.9CE97970184AD8D8.idx
?? .cache/clangd/index/get_distance_to_surface.cpp.C6FA80CB76043D41.idx
?? .cache/clangd/index/get_geodesic_distance.cpp.D1B19D2ADB17334F.idx
?? .cache/clangd/index/get_geodesic_distance.cpp.F573EE9A7B78961E.idx
?? .cache/clangd/index/get_prevNorm.cpp.721FD5AFF7C1D1E8.idx
?? .cache/clangd/index/get_prevNorm.cpp.9F3A592FDD47A064.idx
?? .cache/clangd/index/get_prevSurv.cpp.7E87C920C709D9E2.idx
?? .cache/clangd/index/get_prevSurv.cpp.B8F8ABDB13308443.idx
?? .cache/clangd/index/hasIntangibles.cpp.750DA381622252E8.idx
?? .cache/clangd/index/hasIntangibles.cpp.C403F7E13B382785.idx
?? .cache/clangd/index/implicitlipid_reactions.hpp.BE41CF893F3E3202.idx
?? .cache/clangd/index/init_NboundPairs.cpp.160EDCA1D2DDACD1.idx
?? .cache/clangd/index/init_NboundPairs.cpp.E64DA4E6B33DB78F.idx
?? .cache/clangd/index/init_association_events.cpp.656FD99258AD8593.idx
?? .cache/clangd/index/init_association_events.cpp.78E68E54C27C808B.idx
?? .cache/clangd/index/init_counterCopyNums.cpp.169C9D87C9522E8D.idx
?? .cache/clangd/index/init_counterCopyNums.cpp.369A76E22ABF6B4E.idx
?? .cache/clangd/index/init_print_dimers.cpp.0E5397223F7F851A.idx
?? .cache/clangd/index/init_print_dimers.cpp.7A0BB8D3B7ACF8C7.idx
?? .cache/clangd/index/init_speciesFile.cpp.0EE2A3920D1A1EAF.idx
?? .cache/clangd/index/init_speciesFile.cpp.CE2403885BC44A99.idx
?? .cache/clangd/index/initialize_complex.cpp.58914BDD0A132399.idx
?? .cache/clangd/index/initialize_complex.cpp.B12D093B13832147.idx
?? .cache/clangd/index/initialize_molecule.cpp.3C498FF5B7FA505E.idx
?? .cache/clangd/index/initialize_molecule.cpp.95EE1954C436E306.idx
?? .cache/clangd/index/initialize_molecule_after_reaction.cpp.B6549DDE11D54A33.idx
?? .cache/clangd/index/initialize_molecule_after_reaction.cpp.CC85C6F5EBB17AC5.idx
?? .cache/clangd/index/initialize_molecule_after_transmission_reaction.cpp.B414F109F5F2F89B.idx
?? .cache/clangd/index/initialize_molecule_after_transmission_reaction.cpp.E430319412D81C54.idx
?? .cache/clangd/index/initialize_parameters_for_implicitlipid_model.cpp.22FBDEAA95954A37.idx
?? .cache/clangd/index/initialize_parameters_for_implicitlipid_model.cpp.33812D28E82C97CB.idx
?? .cache/clangd/index/initialize_states.cpp.5C441F6474FB6FFC.idx
?? .cache/clangd/index/initialize_states.cpp.6435AFD291B3653F.idx
?? .cache/clangd/index/integrator.cpp.821750084BDEC9E4.idx
?? .cache/clangd/index/integrator.cpp.DF0F7B03A7E90F05.idx
?? .cache/clangd/index/io.hpp.8FE3FE838B7920FC.idx
?? .cache/clangd/index/isReactant.cpp.3161D866610F7246.idx
?? .cache/clangd/index/isReactant.cpp.568AE496327CB29F.idx
?? .cache/clangd/index/math_functions.cpp.90EFF53D13A8D1A8.idx
?? .cache/clangd/index/math_functions.cpp.DE2E8C9B8E4DC3CD.idx
?? .cache/clangd/index/math_functions.hpp.1B4F6843444C0DD4.idx
?? .cache/clangd/index/matrix.hpp.2CBD37937179E515.idx
?? .cache/clangd/index/matrix_functions.cpp.96DE76A2160DF10B.idx
?? .cache/clangd/index/matrix_functions.cpp.FAC2AAE0B7737A30.idx
?? .cache/clangd/index/measure_complex_displacement.cpp.0BB1F795AF1055EB.idx
?? .cache/clangd/index/measure_complex_displacement.cpp.29EEDE0A9A6A3454.idx
?? .cache/clangd/index/measure_overlap_free_protein_interfaces.cpp.3746C0BE56149FC5.idx
?? .cache/clangd/index/measure_overlap_free_protein_interfaces.cpp.F84E6B2A389C4918.idx
?? .cache/clangd/index/measure_overlap_protein_interfaces.cpp.96FEED36AEB5BF51.idx
?? .cache/clangd/index/measure_overlap_protein_interfaces.cpp.AC627E953452A87D.idx
?? .cache/clangd/index/nerdss.cpp.80B63EE2100376EF.idx
?? .cache/clangd/index/nerdss.cpp.98F4FD2E6C3E5F90.idx
?? .cache/clangd/index/norm_function.cpp.672A0510FB9DB49C.idx
?? .cache/clangd/index/norm_function.cpp.CC096EA8EE978C74.idx
?? .cache/clangd/index/omega_rotation.cpp.0D4288A3494E3C5F.idx
?? .cache/clangd/index/omega_rotation.cpp.20880CC9A02FB952.idx
?? .cache/clangd/index/orient_crds_to_template.cpp.41AF77913A823B96.idx
?? .cache/clangd/index/orient_crds_to_template.cpp.994CA3097D6421EF.idx
?? .cache/clangd/index/parse_command.cpp.25146AACF03FAAA3.idx
?? .cache/clangd/index/parse_command.cpp.3E2D6E42F9FAD93D.idx
?? .cache/clangd/index/parse_input.cpp.B00710BD2C81E1F1.idx
?? .cache/clangd/index/parse_input.cpp.BD0494BEFFE0F147.idx
?? .cache/clangd/index/parse_input_array.cpp.112811B41C24C6BD.idx
?? .cache/clangd/index/parse_input_array.cpp.87A69827A61A915D.idx
?? .cache/clangd/index/parse_molFile.cpp.2EA7734D8B936620.idx
?? .cache/clangd/index/parse_molFile.cpp.8695447D0075BDCE.idx
?? .cache/clangd/index/parse_molecule_bngl.cpp.78E8E3AC5CB39E26.idx
?? .cache/clangd/index/parse_molecule_bngl.cpp.9BD78A04AADDC2A4.idx
?? .cache/clangd/index/parse_observable.cpp.3015F4BE8974906F.idx
?? .cache/clangd/index/parse_observable.cpp.69451BB6793785EB.idx
?? .cache/clangd/index/parse_reaction.cpp.01B52812B8D37DB4.idx
?? .cache/clangd/index/parse_reaction.cpp.341EB28634FFCD91.idx
?? .cache/clangd/index/parse_states.cpp.5A62E9149BB859DF.idx
?? .cache/clangd/index/parse_states.cpp.B9CB12992CF4A0A5.idx
?? .cache/clangd/index/parser_functions.hpp.C6AC7DDAB8B9AD5F.idx
?? .cache/clangd/index/passocF.cpp.5B71984DFDBBFDBF.idx
?? .cache/clangd/index/passocF.cpp.AAFBFD0C5EBDE062.idx
?? .cache/clangd/index/passocF_1D.cpp.77B4214CF419885C.idx
?? .cache/clangd/index/perform_bimolecular_state_change.cpp.333EFC1BDC163D63.idx
?? .cache/clangd/index/perform_bimolecular_state_change.cpp.41C1D83B440B0C8A.idx
?? .cache/clangd/index/perform_bimolecular_state_change_box.cpp.CE3CEDF3129DFF1D.idx
?? .cache/clangd/index/perform_bimolecular_state_change_box.cpp.F0E282660FA3685C.idx
?? .cache/clangd/index/perform_bimolecular_state_change_sphere.cpp.AF9E179F35BC8C0A.idx
?? .cache/clangd/index/perform_bimolecular_state_change_sphere.cpp.D7FCA77491AFB616.idx
?? .cache/clangd/index/perform_implicitlipid_state_change.cpp.309F94221A9D606D.idx
?? .cache/clangd/index/perform_implicitlipid_state_change.cpp.3898CBB324762504.idx
?? .cache/clangd/index/perform_implicitlipid_state_change_box.cpp.21ABAEBCFC800E01.idx
?? .cache/clangd/index/perform_implicitlipid_state_change_box.cpp.C8077F00E11A2386.idx
?? .cache/clangd/index/perform_implicitlipid_state_change_sphere.cpp.4C9D562D02BB33E1.idx
?? .cache/clangd/index/perform_implicitlipid_state_change_sphere.cpp.B5FBA8045B7A81C3.idx
?? .cache/clangd/index/perform_transmission_reaction.cpp.8E79A0B3C8027EB3.idx
?? .cache/clangd/index/perform_transmission_reaction.cpp.D92311C4B8220054.idx
?? .cache/clangd/index/phi_rotation.cpp.DC11EAC0CA07F981.idx
?? .cache/clangd/index/phi_rotation.cpp.FD90D2E6A30FC2A1.idx
?? .cache/clangd/index/pir_function.cpp.91674D199853ADDB.idx
?? .cache/clangd/index/pir_function.cpp.E4F055D38B805296.idx
?? .cache/clangd/index/pirr_pfree_ratio_psF.cpp.4EC38F841FF03218.idx
?? .cache/clangd/index/pirr_pfree_ratio_psF.cpp.A691C31BAB345CE5.idx
?? .cache/clangd/index/pirr_pfree_ratio_psF_1D.cpp.A68F9C0908B44A34.idx
?? .cache/clangd/index/populate_reaction_lists.cpp.304E85BC5473A9A8.idx
?? .cache/clangd/index/populate_reaction_lists.cpp.5D9D75E42416741B.idx
?? .cache/clangd/index/print_association_events.cpp.7A5E4282ACA3BEB9.idx
?? .cache/clangd/index/print_association_events.cpp.DF9D99C5A6413683.idx
?? .cache/clangd/index/print_complex_hist.cpp.B48E54801E501868.idx
?? .cache/clangd/index/print_complex_hist.cpp.C58EEE4813E17764.idx
?? .cache/clangd/index/print_dimers.cpp.0914723F53EEE3E4.idx
?? .cache/clangd/index/print_dimers.cpp.A2D1B639045737BB.idx
?? .cache/clangd/index/print_system_information.cpp.1232EF4EE6F75EC5.idx
?? .cache/clangd/index/print_system_information.cpp.958D24918492AB2B.idx
?? .cache/clangd/index/rand_gsl.cpp.42664C848C74200A.idx
?? .cache/clangd/index/rand_gsl.cpp.80B80CBE7885345A.idx
?? .cache/clangd/index/rand_gsl.hpp.ACEB269AA9B0E2AB.idx
?? .cache/clangd/index/read_bonds.cpp.2DEAF93439058A98.idx
?? .cache/clangd/index/read_bonds.cpp.5956AAD7DA6DD01D.idx
?? .cache/clangd/index/read_boolean.cpp.3EE6F4D4190E4E3D.idx
?? .cache/clangd/index/read_boolean.cpp.84083415C1F5909C.idx
?? .cache/clangd/index/read_internal_coordinates.cpp.1E850579D8F6DFF2.idx
?? .cache/clangd/index/read_internal_coordinates.cpp.8F6468D3BEDDFF0D.idx
?? .cache/clangd/index/read_restart.cpp.1EF5DBBFD7224790.idx
?? .cache/clangd/index/read_restart.cpp.EF13FDB3E22CAD27.idx
?? .cache/clangd/index/reflect_complex_compartment.cpp.1224F15970C63878.idx
?? .cache/clangd/index/reflect_complex_compartment.cpp.8E39336CAD854B7A.idx
?? .cache/clangd/index/reflect_complex_rad_rot.cpp.83470A1201D5DF43.idx
?? .cache/clangd/index/reflect_complex_rad_rot.cpp.ADDFECB32A8D664A.idx
?? .cache/clangd/index/reflect_complex_rad_rot_box.cpp.66DEFA766CE9D1F2.idx
?? .cache/clangd/index/reflect_complex_rad_rot_box.cpp.EB1E345E25642D52.idx
?? .cache/clangd/index/reflect_complex_rad_rot_sphere.cpp.CCDAA21419BDEB14.idx
?? .cache/clangd/index/reflect_complex_rad_rot_sphere.cpp.DB0B097E1535342A.idx
?? .cache/clangd/index/reflect_functions.hpp.F43662F0B421E913.idx
?? .cache/clangd/index/reflect_traj_check_span.cpp.2CB1A01A393740C1.idx
?? .cache/clangd/index/reflect_traj_check_span.cpp.4E7126A71A38F442.idx
?? .cache/clangd/index/reflect_traj_check_span_box.cpp.A62E3DFA1A2699F8.idx
?? .cache/clangd/index/reflect_traj_check_span_box.cpp.F75D584F03D368D1.idx
?? .cache/clangd/index/reflect_traj_check_span_sphere.cpp.3E51FA003DA5CBDF.idx
?? .cache/clangd/index/reflect_traj_check_span_sphere.cpp.DC74740111A9922B.idx
?? .cache/clangd/index/reflect_traj_complex_compartment.cpp.CA24B2A6CD199AAB.idx
?? .cache/clangd/index/reflect_traj_complex_compartment.cpp.E191C793A1DE239C.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot.cpp.9B45E05045508237.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot.cpp.EA66AD9BBB3D5099.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_box.cpp.167D430E49C19648.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_box.cpp.C6CD9CC8DBD39373.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_nocheck.cpp.1D91C7BF93825375.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_nocheck.cpp.B26AD61B263755ED.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_box.cpp.5E56F80FE9B6BFE0.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_box.cpp.B8FB5D899D7207E6.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_sphere.cpp.346A7DF421B7A0D0.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_sphere.cpp.579412256539EA4F.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_sphere.cpp.8C26C2DB6679076A.idx
?? .cache/clangd/index/reflect_traj_complex_rad_rot_sphere.cpp.A999161725CAD16F.idx
?? .cache/clangd/index/reflect_traj_tmp_crds.cpp.300F4C7090381D20.idx
?? .cache/clangd/index/reflect_traj_tmp_crds.cpp.FC4AB1972B08A76C.idx
?? .cache/clangd/index/reflect_traj_tmp_crds_box.cpp.877F4FD732DC927A.idx
?? .cache/clangd/index/reflect_traj_tmp_crds_box.cpp.C5087CF17B06EF9D.idx
?? .cache/clangd/index/reflect_traj_tmp_crds_compartment.cpp.5F4D7E7C3A8C31C6.idx
?? .cache/clangd/index/reflect_traj_tmp_crds_compartment.cpp.DC9A160C71AF3260.idx
?? .cache/clangd/index/reflect_traj_tmp_crds_sphere.cpp.15407B68AAFDA41A.idx
?? .cache/clangd/index/reflect_traj_tmp_crds_sphere.cpp.A0593DC7A5EB3B6E.idx
?? .cache/clangd/index/remove_comment.cpp.52FDE4B10B1B7D5F.idx
?? .cache/clangd/index/remove_comment.cpp.80E4A2E2672773D8.idx
?? .cache/clangd/index/requiresSignFlip.cpp.12FA43DABADF6829.idx
?? .cache/clangd/index/requiresSignFlip.cpp.BFA627775DD27381.idx
?? .cache/clangd/index/reverse_rotation.cpp.2357217350D1937F.idx
?? .cache/clangd/index/reverse_rotation.cpp.FF87CFB3C78B8A6F.idx
?? .cache/clangd/index/rotate.cpp.93616F7AAD2ED0DF.idx
?? .cache/clangd/index/rotate.cpp.FFB632E877B19079.idx
?? .cache/clangd/index/save_mem_orientation.cpp.4318E508B6D3CCFF.idx
?? .cache/clangd/index/save_mem_orientation.cpp.EF1D0CA7F40CF579.idx
?? .cache/clangd/index/set_rMaxLimit.cpp.6817D3F55C8D917B.idx
?? .cache/clangd/index/set_rMaxLimit.cpp.EF935BE6CE9C368D.idx
?? .cache/clangd/index/shared_reaction_functions.hpp.78F5129E8A21E9AE.idx
?? .cache/clangd/index/size_lookup.cpp.3D203DD11EA0F57F.idx
?? .cache/clangd/index/size_lookup.cpp.F5E59B74A1B92311.idx
?? .cache/clangd/index/subtract_off_com_position.cpp.51D15F5904BDEBCA.idx
?? .cache/clangd/index/subtract_off_com_position.cpp.74131E55141860B3.idx
?? .cache/clangd/index/survival_function.cpp.137F96F31D4F39B5.idx
?? .cache/clangd/index/survival_function.cpp.61AD740FDA38C84E.idx
?? .cache/clangd/index/sweep_separation_complex_rot.cpp.138C679ECC97BD67.idx
?? .cache/clangd/index/sweep_separation_complex_rot.cpp.2FDE65E038AA1F8D.idx
?? .cache/clangd/index/sweep_separation_complex_rot_box.cpp.AB5601CE2E1F47F6.idx
?? .cache/clangd/index/sweep_separation_complex_rot_box.cpp.FB720D7F93035D3B.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest.cpp.2AC6B4BD97E8D867.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest.cpp.8DD995669F88DE3A.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_box.cpp.67AEB26631344EF9.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_box.cpp.84BB62414E7A9EB1.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_cluster.cpp.6AA1F8C1C98421CB.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_cluster.cpp.FF8159554073BAEC.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_box.cpp.D34CD4E3BFF1A907.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_box.cpp.D91799B387F49E2B.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_sphere.cpp.7CB8285DB34C9BC3.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_sphere.cpp.C096DD9B75FEF2E6.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_sphere.cpp.3544A93E7245AA23.idx
?? .cache/clangd/index/sweep_separation_complex_rot_memtest_sphere.cpp.7533B7D5FF46A997.idx
?? .cache/clangd/index/sweep_separation_complex_rot_sphere.cpp.23E5A7C15BAD5C20.idx
?? .cache/clangd/index/sweep_separation_complex_rot_sphere.cpp.E96C225D7E451FB1.idx
?? .cache/clangd/index/system_setup.hpp.975D395A733EAA9B.idx
?? .cache/clangd/index/theta_rotation.cpp.725C02B950F89CBB.idx
?? .cache/clangd/index/theta_rotation.cpp.E7F8FCBAAF5EDD25.idx
?? .cache/clangd/index/tracing.hpp.1E86E9BD1E155318.idx
?? .cache/clangd/index/track_association_events.cpp.5EC1BB0518EF4C88.idx
?? .cache/clangd/index/track_association_events.cpp.BC9C3A5899771857.idx
?? .cache/clangd/index/trajectory_functions.hpp.AC79A5E4DB6676EE.idx
?? .cache/clangd/index/transform.cpp.0799BA799CB62539.idx
?? .cache/clangd/index/transform.cpp.A294C2E06FEAA8EC.idx
?? .cache/clangd/index/unimolecular_reactions.hpp.E68FB70E56D09C26.idx
?? .cache/clangd/index/update_Nboundpairs.cpp.2F87FB400C0CBE76.idx
?? .cache/clangd/index/update_Nboundpairs.cpp.60C03811F69BE3CC.idx
?? .cache/clangd/index/update_complex_tmp_com_crds.cpp.0019E67F368F77F2.idx
?? .cache/clangd/index/update_complex_tmp_com_crds.cpp.B455D89D8F568DB2.idx
?? .cache/clangd/index/write_NboundPairs.cpp.24ECA0CF0D20459B.idx
?? .cache/clangd/index/write_NboundPairs.cpp.AC875BA9E909A18F.idx
?? .cache/clangd/index/write_all_species.cpp.29EB2818CDE699E7.idx
?? .cache/clangd/index/write_all_species.cpp.5CEA5725D26E7FE7.idx
?? .cache/clangd/index/write_complex_components.cpp.7C6E3BF548909F1B.idx
?? .cache/clangd/index/write_complex_components.cpp.E10F59063DA2209C.idx
?? .cache/clangd/index/write_complex_crds.cpp.5D2BB1B5982AB2E4.idx
?? .cache/clangd/index/write_complex_crds.cpp.F956FAA0CC551A0D.idx
?? .cache/clangd/index/write_crds.cpp.87B153880311ECEC.idx
?? .cache/clangd/index/write_crds.cpp.8EC8485290ED5644.idx
?? .cache/clangd/index/write_mol_iface.cpp.A132C69B74F5FBA4.idx
?? .cache/clangd/index/write_mol_iface.cpp.A70EADBF89E29E06.idx
?? .cache/clangd/index/write_observables.cpp.215B012073FF5ADC.idx
?? .cache/clangd/index/write_observables.cpp.DFEFE8F8472E1B98.idx
?? .cache/clangd/index/write_pdb.cpp.7AE30E5B95F74C2D.idx
?? .cache/clangd/index/write_pdb.cpp.E8E9139B4AFC19F2.idx
?? .cache/clangd/index/write_psf.cpp.24855AE09A3509C7.idx
?? .cache/clangd/index/write_psf.cpp.A63F54587DF31E90.idx
?? .cache/clangd/index/write_restart.cpp.80E381FC88A4FACC.idx
?? .cache/clangd/index/write_restart.cpp.D1EF0F3CE1EE6D2A.idx
?? .cache/clangd/index/write_timestep_information.cpp.217117A0DC33022D.idx
?? .cache/clangd/index/write_timestep_information.cpp.31443B470577CB15.idx
?? .cache/clangd/index/write_traj.cpp.0CE11FFD455FCD3A.idx
?? .cache/clangd/index/write_traj.cpp.260F4AC3A3565830.idx
?? .cache/clangd/index/write_transition.cpp.1BC5F68470F3D951.idx
?? .cache/clangd/index/write_transition.cpp.1BF6ED8E56F75945.idx
?? .cache/clangd/index/write_xyz.cpp.2C4B77EE4DC45A9D.idx
?? .cache/clangd/index/write_xyz.cpp.AF1F327DB8F0F77A.idx
?? .cache/clangd/index/write_xyz_assoc.cpp.0F06761FFD82FEB5.idx
?? .cache/clangd/index/write_xyz_assoc.cpp.0FEBE136BE2B9E57.idx
?? .cache/clangd/index/write_xyz_assoc_cout.cpp.2435F8AAB24E9634.idx
?? .cache/clangd/index/write_xyz_assoc_cout.cpp.2C6D704A7AF33BA3.idx
?? .github/workflows/c-cpp.yml
?? REPLACEMENT_REPORT.md
```

### git diff --stat

```
 .gitignore                                         |   26 +-
 EXEs/nerdss.cpp                                    |    2 +-
 Makefile                                           |    2 +-
 Makefile_old                                       |    0
 include/boundary_conditions/reflect_functions.hpp  |    2 +-
 include/classes/class_Molecule_Complex.hpp         |    2 +-
 sample_inputs/VALIDATE_SUITE/hetTrimer/PDB/0.pdb   | 3002 -------
 .../hetTrimer/RESTARTS/restart120000.dat           |  Bin 1765436 -> 0 bytes
 .../hetTrimer/RESTARTS/restart180000.dat           |  Bin 2306993 -> 0 bytes
 .../hetTrimer/RESTARTS/restart60000.dat            |  Bin 1220009 -> 0 bytes
 .../VALIDATE_SUITE/hetTrimer/assoc_dissoc_time.dat |    0
 .../VALIDATE_SUITE/hetTrimer/bound_pair_time.dat   |   46 -
 .../VALIDATE_SUITE/hetTrimer/copy_numbers_time.dat |   46 -
 .../hetTrimer/event_counters_time.dat              |  131 -
 .../hetTrimer/histogram_complexes_time.dat         |  266 -
 .../VALIDATE_SUITE/hetTrimer/initial_crds.xyz      | 2999 -------
 .../VALIDATE_SUITE/hetTrimer/mono_dimer_time.dat   |   46 -
 .../VALIDATE_SUITE/hetTrimer/observables_time.dat  |    0
 sample_inputs/VALIDATE_SUITE/hetTrimer/system.psf  | 3508 --------
 .../VALIDATE_SUITE/hetTrimer/trajectory.xyz        | 2999 -------
 .../hetTrimer/transition_matrix_time.dat           | 8505 ------------------
 sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/0.pdb  | 3003 -------
 .../VALIDATE_SUITE/homoTrimer/PDB/1000000.pdb      | 3003 -------
 .../VALIDATE_SUITE/homoTrimer/PDB/1999999.pdb      | 3003 -------
 .../VALIDATE_SUITE/homoTrimer/PDB/599999.pdb       | 3003 -------
 .../VALIDATE_SUITE/homoTrimer/PDB/9999.pdb         | 3003 -------
 .../homoTrimer/RESTARTS/restart1000.dat            |  Bin 605789 -> 0 bytes
 .../homoTrimer/RESTARTS/restart1000000.dat         |  Bin 666685 -> 0 bytes
 .../homoTrimer/RESTARTS/restart120000.dat          |  Bin 657238 -> 0 bytes
 .../homoTrimer/RESTARTS/restart1200000.dat         |  Bin 651798 -> 0 bytes
 .../homoTrimer/RESTARTS/restart1400000.dat         |  Bin 684642 -> 0 bytes
 .../homoTrimer/RESTARTS/restart1600000.dat         |  Bin 666670 -> 0 bytes
 .../homoTrimer/RESTARTS/restart180000.dat          |  Bin 651659 -> 0 bytes
 .../homoTrimer/RESTARTS/restart1800000.dat         |  Bin 677870 -> 0 bytes
 .../homoTrimer/RESTARTS/restart2000.dat            |  Bin 590939 -> 0 bytes
 .../homoTrimer/RESTARTS/restart200000.dat          |  Bin 660997 -> 0 bytes
 .../homoTrimer/RESTARTS/restart240000.dat          |  Bin 651538 -> 0 bytes
 .../homoTrimer/RESTARTS/restart3000.dat            |  Bin 588434 -> 0 bytes
 .../homoTrimer/RESTARTS/restart300000.dat          |  Bin 668739 -> 0 bytes
 .../homoTrimer/RESTARTS/restart360000.dat          |  Bin 660725 -> 0 bytes
 .../homoTrimer/RESTARTS/restart4000.dat            |  Bin 590327 -> 0 bytes
 .../homoTrimer/RESTARTS/restart400000.dat          |  Bin 671970 -> 0 bytes
 .../homoTrimer/RESTARTS/restart420000.dat          |  Bin 684440 -> 0 bytes
 .../homoTrimer/RESTARTS/restart480000.dat          |  Bin 664606 -> 0 bytes
 .../homoTrimer/RESTARTS/restart5000.dat            |  Bin 587655 -> 0 bytes
 .../homoTrimer/RESTARTS/restart540000.dat          |  Bin 666147 -> 0 bytes
 .../homoTrimer/RESTARTS/restart6000.dat            |  Bin 592663 -> 0 bytes
 .../homoTrimer/RESTARTS/restart60000.dat           |  Bin 675259 -> 0 bytes
 .../homoTrimer/RESTARTS/restart600000.dat          |  Bin 675517 -> 0 bytes
 .../homoTrimer/RESTARTS/restart7000.dat            |  Bin 591875 -> 0 bytes
 .../homoTrimer/RESTARTS/restart8000.dat            |  Bin 585524 -> 0 bytes
 .../homoTrimer/RESTARTS/restart800000.dat          |  Bin 664237 -> 0 bytes
 .../homoTrimer/RESTARTS/restart9000.dat            |  Bin 592661 -> 0 bytes
 .../homoTrimer/assoc_dissoc_time.dat               |    0
 .../VALIDATE_SUITE/homoTrimer/bound_pair_time.dat  |  402 -
 .../homoTrimer/copy_numbers_time.dat               |  402 -
 .../homoTrimer/event_counters_time.dat             | 1201 ---
 .../VALIDATE_SUITE/homoTrimer/final_coords.xyz     | 3002 -------
 .../homoTrimer/histogram_complexes_time.dat        | 1602 ----
 .../VALIDATE_SUITE/homoTrimer/initial_crds.xyz     | 3002 -------
 .../VALIDATE_SUITE/homoTrimer/mono_dimer_time.dat  |  402 -
 .../VALIDATE_SUITE/homoTrimer/observables_time.dat |    0
 .../VALIDATE_SUITE/homoTrimer/restart.dat          |  Bin 668183 -> 0 bytes
 sample_inputs/VALIDATE_SUITE/homoTrimer/system.psf | 3011 -------
 .../VALIDATE_SUITE/homoTrimer/trajectory.xyz       | 9006 --------------------
 .../homoTrimer/transition_matrix_time.dat          | 1203 ---
 sample_inputs/VALIDATE_SUITE/trimer/PDB/0.pdb      |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/100000.pdb |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1000000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1100000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1200000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1300000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1400000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1500000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1600000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1700000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1800000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/1900000.pdb          |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/200000.pdb |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2000000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2100000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2200000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2300000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2400000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2500000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2600000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2700000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2800000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/2900000.pdb          |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/300000.pdb |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3000000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3100000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3200000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3300000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3400000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3500000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3600000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3700000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3800000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/3900000.pdb          |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/400000.pdb |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4000000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4100000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4200000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4300000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4400000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4500000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4600000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4700000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4800000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4900000.pdb          |  905 --
 .../VALIDATE_SUITE/trimer/PDB/4999999.pdb          |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/500000.pdb |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/600000.pdb |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/700000.pdb |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/800000.pdb |  905 --
 sample_inputs/VALIDATE_SUITE/trimer/PDB/900000.pdb |  905 --
 .../VALIDATE_SUITE/trimer/assoc_dissoc_time.dat    |    0
 .../VALIDATE_SUITE/trimer/bound_pair_time.dat      |   24 -
 .../VALIDATE_SUITE/trimer/copy_numbers_time.dat    |   24 -
 .../VALIDATE_SUITE/trimer/event_counters_time.dat  |   63 -
 .../VALIDATE_SUITE/trimer/final_coords.xyz         |  902 --
 .../trimer/histogram_complexes_time.dat            |  178 -
 .../VALIDATE_SUITE/trimer/initial_crds.xyz         |  902 --
 .../VALIDATE_SUITE/trimer/mono_dimer_time.dat      |   24 -
 .../VALIDATE_SUITE/trimer/observables_time.dat     |    0
 sample_inputs/VALIDATE_SUITE/trimer/restart.dat    |  Bin 1322863 -> 0 bytes
 sample_inputs/VALIDATE_SUITE/trimer/system.psf     | 1061 ---
 sample_inputs/VALIDATE_SUITE/trimer/trajectory.xyz |  902 --
 .../trimer/transition_matrix_time.dat              |  189 -
 src/classes/class_Molecule_Complex.cpp             |   97 +-
 src/classes/class_Parameters.cpp                   |   27 +-
 src/classes/class_Vector.cpp                       |    2 +-
 src/io/read_restart.cpp                            |    4 +
 src/io/write_restart.cpp                           |    1 +
 src/reactions/associate_box.cpp                    |   31 +-
 ...check_for_unimolecular_reactions_population.cpp |    6 +-
 ...termine_1D_bimolecular_reaction_probability.cpp |   14 +-
 src/reactions/pirr_pfree_ratio_psF_1D.cpp          |    2 +-
 .../calculate_update_position_interface.cpp        |    2 +-
 .../create_complex_propagation_vectors.cpp         |   20 +-
 system.psf                                         |   15 -
 text_search.bsh                                    |    0
 143 files changed, 154 insertions(+), 110321 deletions(-)
```

## Added Files

- `.cache/clangd/index/2D_reaction_table_functions.hpp.E01AB59A14764AC9.idx`
- `.cache/clangd/index/DDpirr_pfree_ratio_ps.cpp.596F1BDD8B49E7DF.idx`
- `.cache/clangd/index/DDpirr_pfree_ratio_ps.cpp.C4682E31B2684260.idx`
- `.cache/clangd/index/Faddeeva.cpp.599BF93AF027009F.idx`
- `.cache/clangd/index/Faddeeva.cpp.FEB6F97561210959.idx`
- `.cache/clangd/index/Faddeeva.hpp.31C1D40B1BF1B9A0.idx`
- `.cache/clangd/index/angleSignIsCorrect.cpp.346BCB05078F6C0C.idx`
- `.cache/clangd/index/angleSignIsCorrect.cpp.9AD06B6E47C4053B.idx`
- `.cache/clangd/index/areInVicinity.cpp.4B463C5EFAF9A5CA.idx`
- `.cache/clangd/index/areInVicinity.cpp.4D8CC006E6E845BE.idx`
- `.cache/clangd/index/areParallel.cpp.FBDA2F60160EC332.idx`
- `.cache/clangd/index/areParallel.cpp.FCD65FAE99F26C0C.idx`
- `.cache/clangd/index/areSameAngle.cpp.7BB77C900F78C285.idx`
- `.cache/clangd/index/areSameAngle.cpp.D1C9C1D1812B119A.idx`
- `.cache/clangd/index/areSameExceptState.cpp.688A3B589E0602E7.idx`
- `.cache/clangd/index/areSameExceptState.cpp.D3D9FFD38DD5AEE7.idx`
- `.cache/clangd/index/associate.cpp.88C8B5BDAE14774A.idx`
- `.cache/clangd/index/associate.cpp.B2EB00E69DB9CAE7.idx`
- `.cache/clangd/index/associate_ImplicitLipid.cpp.602E5A2F2AC34B7C.idx`
- `.cache/clangd/index/associate_ImplicitLipid.cpp.994C548EC5BAEDF0.idx`
- `.cache/clangd/index/associate_ImplicitLipid_box.cpp.1E1E65CA086E3A44.idx`
- `.cache/clangd/index/associate_ImplicitLipid_box.cpp.AE174A59E0A21AB6.idx`
- `.cache/clangd/index/associate_ImplicitLipid_sphere.cpp.48EEC1BA49BC613F.idx`
- `.cache/clangd/index/associate_ImplicitLipid_sphere.cpp.F1D082D6F78F82E8.idx`
- `.cache/clangd/index/associate_box.cpp.7C9F70CFE14BFF99.idx`
- `.cache/clangd/index/associate_box.cpp.80218D6313AA983B.idx`
- `.cache/clangd/index/associate_sphere.cpp.A1E17496A6B669BB.idx`
- `.cache/clangd/index/associate_sphere.cpp.F0A5B19AAFC4B23E.idx`
- `.cache/clangd/index/association.hpp.13A9340F987A3FAF.idx`
- `.cache/clangd/index/bimolecular_reactions.hpp.A24B66D3C05A8948.idx`
- `.cache/clangd/index/break_interaction.cpp.5212335BC8D27174.idx`
- `.cache/clangd/index/break_interaction.cpp.9694D553B30D7489.idx`
- `.cache/clangd/index/break_interaction_implicitlipid.cpp.4404381750380388.idx`
- `.cache/clangd/index/break_interaction_implicitlipid.cpp.4D955D83E82D7C1B.idx`
- `.cache/clangd/index/calc_angular_displacement.cpp.327AEAAA7C0A1090.idx`
- `.cache/clangd/index/calc_angular_displacement.cpp.AEFA51D95ABB4532.idx`
- `.cache/clangd/index/calc_one_angular_displacement.cpp.92AB01A168D1960B.idx`
- `.cache/clangd/index/calc_one_angular_displacement.cpp.9B158482797B861E.idx`
- `.cache/clangd/index/calc_pirr.cpp.3AF719EBB17BB3A4.idx`
- `.cache/clangd/index/calc_pirr.cpp.C282BA08E48EEC99.idx`
- `.cache/clangd/index/calculate_omega.cpp.26CA1EF283935390.idx`
- `.cache/clangd/index/calculate_omega.cpp.720D33FAD6633F7C.idx`
- `.cache/clangd/index/calculate_phi.cpp.BDADAD0C353D6443.idx`
- `.cache/clangd/index/calculate_phi.cpp.E41D6C7D710850BC.idx`
- `.cache/clangd/index/calculate_update_position_interface.cpp.0E1D62689CE6FDA1.idx`
- `.cache/clangd/index/calculate_update_position_interface.cpp.3D1DD3C222C67D9B.idx`
- `.cache/clangd/index/check_bases.cpp.705FEED0F752FED5.idx`
- `.cache/clangd/index/check_bases.cpp.92B255123CB6A911.idx`
- `.cache/clangd/index/check_bimolecular_reactions.cpp.2BAFDCC67DD415A8.idx`
- `.cache/clangd/index/check_bimolecular_reactions.cpp.3D7D5567B8E471DB.idx`
- `.cache/clangd/index/check_compartment_reaction.cpp.82A072B6658F1087.idx`
- `.cache/clangd/index/check_compartment_reaction.cpp.DE4FCA58E2B79A6C.idx`
- `.cache/clangd/index/check_dissociation.cpp.087600F372E7F5E3.idx`
- `.cache/clangd/index/check_dissociation.cpp.52F4CD98B1355803.idx`
- `.cache/clangd/index/check_dissociation_implicitlipid.cpp.7194057AA6673E4B.idx`
- `.cache/clangd/index/check_dissociation_implicitlipid.cpp.9FAFF3FDFE4C886A.idx`
- `.cache/clangd/index/check_for_state_change.cpp.0D2C35DEBE67D58A.idx`
- `.cache/clangd/index/check_for_state_change.cpp.EC7CA378EA07A354.idx`
- `.cache/clangd/index/check_for_structure_overlap.cpp.0EE5B4861890BCDC.idx`
- `.cache/clangd/index/check_for_structure_overlap.cpp.605A7A3DE1FD173D.idx`
- `.cache/clangd/index/check_for_structure_overlap_system.cpp.21C718F8CAA45044.idx`
- `.cache/clangd/index/check_for_structure_overlap_system.cpp.ED27AC02F58D3B56.idx`
- `.cache/clangd/index/check_for_unimolecular_reactions.cpp.1D95BB89B5AEF292.idx`
- `.cache/clangd/index/check_for_unimolecular_reactions.cpp.E61F3EAE361A703A.idx`
- `.cache/clangd/index/check_for_unimolecular_reactions_population.cpp.064A3C05887F547D.idx`
- `.cache/clangd/index/check_for_unimolecular_reactions_population.cpp.90B110B8E372D347.idx`
- `.cache/clangd/index/check_for_unimolstatechange_reactions.cpp.159A531B8BE6B6F4.idx`
- `.cache/clangd/index/check_for_unimolstatechange_reactions.cpp.4B2B4CE87FDE1D7D.idx`
- `.cache/clangd/index/check_for_valid_states.cpp.1C38B1412F4FE96D.idx`
- `.cache/clangd/index/check_for_valid_states.cpp.37A9ACAD134AEFF5.idx`
- `.cache/clangd/index/check_for_zeroth_order_creation.cpp.6F525395F7241967.idx`
- `.cache/clangd/index/check_for_zeroth_order_creation.cpp.DBBFD5CFD238D079.idx`
- `.cache/clangd/index/check_if_spans.cpp.83E8E1A2DB912F90.idx`
- `.cache/clangd/index/check_if_spans.cpp.B25B0590DCB34027.idx`
- `.cache/clangd/index/check_if_spans_box.cpp.5D284F79C405A37E.idx`
- `.cache/clangd/index/check_if_spans_box.cpp.C5D8664F46DE09B0.idx`
- `.cache/clangd/index/check_if_spans_sphere.cpp.E135304E3588AD48.idx`
- `.cache/clangd/index/check_if_spans_sphere.cpp.E3CBCAB4A1C2E05C.idx`
- `.cache/clangd/index/check_implicit_reactions.cpp.77FF45D85157C8AB.idx`
- `.cache/clangd/index/check_implicit_reactions.cpp.C7A4AB4A551E5101.idx`
- `.cache/clangd/index/class_Cluster.cpp.BD160A36A90F5798.idx`
- `.cache/clangd/index/class_Cluster.cpp.D1477C0585549194.idx`
- `.cache/clangd/index/class_Cluster.hpp.E89C938C2F55233D.idx`
- `.cache/clangd/index/class_Coord.cpp.199DC27C605D7225.idx`
- `.cache/clangd/index/class_Coord.cpp.A28391A9D1673528.idx`
- `.cache/clangd/index/class_Coord.hpp.3AB0AB56D94C613E.idx`
- `.cache/clangd/index/class_MDTimer.cpp.29402414A915C06F.idx`
- `.cache/clangd/index/class_MDTimer.cpp.E55D4F42BFFEE3F6.idx`
- `.cache/clangd/index/class_MDTimer.hpp.6E2DDE8C259E6173.idx`
- `.cache/clangd/index/class_Membrane.hpp.D1CB8CD4E39E5BE1.idx`
- `.cache/clangd/index/class_MolTemplate.cpp.0047A1AC1C5A577C.idx`
- `.cache/clangd/index/class_MolTemplate.cpp.3DBF9D6A5F5AC21E.idx`
- `.cache/clangd/index/class_MolTemplate.hpp.D2862D54454C5B3D.idx`
- `.cache/clangd/index/class_Molecule_Complex.cpp.06DC17412EC80C5B.idx`
- `.cache/clangd/index/class_Molecule_Complex.cpp.DE1C3BB70572DEDF.idx`
- `.cache/clangd/index/class_Molecule_Complex.hpp.A5431CAD3364BAD6.idx`
- `.cache/clangd/index/class_Observable.cpp.3DD6007F75742128.idx`
- `.cache/clangd/index/class_Observable.cpp.F9748ACB8339B771.idx`
- `.cache/clangd/index/class_Observable.hpp.B06702238E12CEA8.idx`
- `.cache/clangd/index/class_Parameters.cpp.6530F1C58E6F75D9.idx`
- `.cache/clangd/index/class_Parameters.cpp.9E55625509C6FE93.idx`
- `.cache/clangd/index/class_Parameters.hpp.824301A6ACE9B9C7.idx`
- `.cache/clangd/index/class_Quat.cpp.E352F5EDC32944E7.idx`
- `.cache/clangd/index/class_Quat.cpp.E8318091DA02AC57.idx`
- `.cache/clangd/index/class_Quat.hpp.76D2E15A147C9793.idx`
- `.cache/clangd/index/class_Rxns.cpp.66F6DFCFFC272689.idx`
- `.cache/clangd/index/class_Rxns.cpp.CA479397F29C8CC0.idx`
- `.cache/clangd/index/class_Rxns.hpp.37F41A5206ACC8CA.idx`
- `.cache/clangd/index/class_SimulVolume.cpp.6346D705697E3012.idx`
- `.cache/clangd/index/class_SimulVolume.cpp.FAA20F64BD9149F7.idx`
- `.cache/clangd/index/class_SimulVolume.hpp.35E9A0BA21F5EA27.idx`
- `.cache/clangd/index/class_Vector.cpp.01B952DD7BFA8EB3.idx`
- `.cache/clangd/index/class_Vector.cpp.AFC2EE1F913CE571.idx`
- `.cache/clangd/index/class_Vector.hpp.D2667BB5B7554218.idx`
- `.cache/clangd/index/class_bngl_parser.hpp.9331DCCCB8184F52.idx`
- `.cache/clangd/index/class_bngl_parser_functions.cpp.55F7A9D1332DC984.idx`
- `.cache/clangd/index/class_bngl_parser_functions.cpp.576239292709E75E.idx`
- `.cache/clangd/index/class_copyCounters.hpp.09CB36E1C86E15FA.idx`
- `.cache/clangd/index/clear_reweight_vecs.cpp.442DD75210EA6839.idx`
- `.cache/clangd/index/clear_reweight_vecs.cpp.C3D0747F4E7A31DC.idx`
- `.cache/clangd/index/com_of_two_tmp_complexes.cpp.A5E45E0CA2588523.idx`
- `.cache/clangd/index/com_of_two_tmp_complexes.cpp.CC0FDF914757913A.idx`
- `.cache/clangd/index/conservedMags.cpp.C940BA855CBAA293.idx`
- `.cache/clangd/index/conservedMags.cpp.E0DF51220127458E.idx`
- `.cache/clangd/index/conservedRigid.cpp.235A441583D8610D.idx`
- `.cache/clangd/index/conservedRigid.cpp.86323E869171B346.idx`
- `.cache/clangd/index/constants.hpp.FDF13DEA60CE30B6.idx`
- `.cache/clangd/index/create_DDMatrices.cpp.B38F0B34E290FACE.idx`
- `.cache/clangd/index/create_DDMatrices.cpp.E8CAA1FA5A352C85.idx`
- `.cache/clangd/index/create_arbitrary_vector.cpp.44CC616EED7BAB7D.idx`
- `.cache/clangd/index/create_arbitrary_vector.cpp.4A85B44D1CA088D6.idx`
- `.cache/clangd/index/create_complex_propagation_vectors.cpp.567D0926945E3DF7.idx`
- `.cache/clangd/index/create_complex_propagation_vectors.cpp.767D759168397DF7.idx`
- `.cache/clangd/index/create_complex_propagation_vectors_on_sphere.cpp.81DF5D883F10092D.idx`
- `.cache/clangd/index/create_complex_propagation_vectors_on_sphere.cpp.C8F9B73F6D69E81B.idx`
- `.cache/clangd/index/create_conjugate_reaction_itrs.cpp.0F50AC44BCE0F615.idx`
- `.cache/clangd/index/create_conjugate_reaction_itrs.cpp.B6929A3863FC2A7C.idx`
- `.cache/clangd/index/create_molecule_and_complex_from_rxn.cpp.AA2ED144BE580052.idx`
- `.cache/clangd/index/create_molecule_and_complex_from_rxn.cpp.F88B63867BEA535F.idx`
- `.cache/clangd/index/create_molecule_and_complex_from_transmission_rxn.cpp.6A453F4C7557C781.idx`
- `.cache/clangd/index/create_molecule_and_complex_from_transmission_rxn.cpp.B4BE3EE0144F450A.idx`
- `.cache/clangd/index/create_normMatrix.cpp.4B4823C0237F61D7.idx`
- `.cache/clangd/index/create_normMatrix.cpp.98BF6B26BF50FF30.idx`
- `.cache/clangd/index/create_pirMatrix.cpp.3B2EAF45E9D23109.idx`
- `.cache/clangd/index/create_pirMatrix.cpp.9801504D1C08EB78.idx`
- `.cache/clangd/index/create_survMatrix.cpp.88418B8D81960E00.idx`
- `.cache/clangd/index/create_survMatrix.cpp.9434F5360EA93885.idx`
- `.cache/clangd/index/determine_1D_bimolecular_reaction_probability.cpp.3124089C1A564B84.idx`
- `.cache/clangd/index/determine_2D_bimolecular_reaction_probability.cpp.51ADD49D63731880.idx`
- `.cache/clangd/index/determine_2D_bimolecular_reaction_probability.cpp.62809B98425D6EF3.idx`
- `.cache/clangd/index/determine_2D_implicitlipid_reaction_probability.cpp.4DB82439422FF127.idx`
- `.cache/clangd/index/determine_2D_implicitlipid_reaction_probability.cpp.D63E9DC35F63DAF1.idx`
- `.cache/clangd/index/determine_3D_bimolecular_reaction_probability.cpp.55EDD4001C5898F4.idx`
- `.cache/clangd/index/determine_3D_bimolecular_reaction_probability.cpp.E52AB0E0C63B2289.idx`
- `.cache/clangd/index/determine_3D_implicitlipid_reaction_probability.cpp.93F32C41AB2FA0C0.idx`
- `.cache/clangd/index/determine_3D_implicitlipid_reaction_probability.cpp.BD52A187A620821A.idx`
- `.cache/clangd/index/determine_bound_iface_index.cpp.57893AD913278BBE.idx`
- `.cache/clangd/index/determine_bound_iface_index.cpp.5D4FDF29AC0B5CAC.idx`
- `.cache/clangd/index/determine_entering_compartment_probability.cpp.0584F27799043DCC.idx`
- `.cache/clangd/index/determine_entering_compartment_probability.cpp.A6057F141E49B096.idx`
- `.cache/clangd/index/determine_exiting_compartment_probability.cpp.0FCF25CCDEED9A14.idx`
- `.cache/clangd/index/determine_exiting_compartment_probability.cpp.1A292F4787C47097.idx`
- `.cache/clangd/index/determine_if_reaction_occurs.cpp.1166F7BD9ACCF1CA.idx`
- `.cache/clangd/index/determine_if_reaction_occurs.cpp.BDDD131A11CA9C74.idx`
- `.cache/clangd/index/determine_iface_indices.cpp.6598A92C6837D473.idx`
- `.cache/clangd/index/determine_iface_indices.cpp.EF0D6C4D7BC35D4E.idx`
- `.cache/clangd/index/determine_normal.cpp.7D409ABAAD08B489.idx`
- `.cache/clangd/index/determine_normal.cpp.E5111DA7D4CD5CEB.idx`
- `.cache/clangd/index/determine_parent_complex.cpp.2A89C8EF52F98AFD.idx`
- `.cache/clangd/index/determine_parent_complex.cpp.494AE6A4CFD159DD.idx`
- `.cache/clangd/index/determine_parent_complex_IL.cpp.2CB3A606FC30A567.idx`
- `.cache/clangd/index/determine_parent_complex_IL.cpp.C5A36D94B6BEBD5F.idx`
- `.cache/clangd/index/determine_rotation_angles.cpp.1F92E765FE04DD04.idx`
- `.cache/clangd/index/determine_rotation_angles.cpp.200E7DE777341D93.idx`
- `.cache/clangd/index/determine_shape_molecule.cpp.0A060223FBC79DFD.idx`
- `.cache/clangd/index/determine_shape_molecule.cpp.1F7C4351A0C5AB7F.idx`
- `.cache/clangd/index/display_all_MolTemplates.cpp.22E69D0CC546E283.idx`
- `.cache/clangd/index/display_all_MolTemplates.cpp.AC2BA8CFCB08679F.idx`
- `.cache/clangd/index/display_all_reactions.cpp.B45C56E1E0A39E93.idx`
- `.cache/clangd/index/display_all_reactions.cpp.D690C19E94611240.idx`
- `.cache/clangd/index/error_handling.hpp.00DB20F9513BFBB4.idx`
- `.cache/clangd/index/evaluate_binding_within_complex.cpp.8F7EEBF22E084286.idx`
- `.cache/clangd/index/evaluate_binding_within_complex.cpp.C39FACEBBB326A95.idx`
- `.cache/clangd/index/eye_candy.cpp.9E7D79E279A61098.idx`
- `.cache/clangd/index/eye_candy.cpp.A6F4FA8C8D902CD8.idx`
- `.cache/clangd/index/find_reaction_rate_state.cpp.2783325BB3D9E0AE.idx`
- `.cache/clangd/index/find_reaction_rate_state.cpp.DD748A03243AC975.idx`
- `.cache/clangd/index/find_which_reaction.cpp.19A721BAEC7E453A.idx`
- `.cache/clangd/index/find_which_reaction.cpp.42C0A8EE96B6859C.idx`
- `.cache/clangd/index/find_which_state_change_reaction.cpp.6823B9A51CDD9CA2.idx`
- `.cache/clangd/index/find_which_state_change_reaction.cpp.DAD941B7FA0A20E6.idx`
- `.cache/clangd/index/functions_for_spherical_system.cpp.13BE19BE8FC5AE2C.idx`
- `.cache/clangd/index/functions_for_spherical_system.cpp.BAF403A10FD19FF3.idx`
- `.cache/clangd/index/functions_for_spherical_system.hpp.44C914A5BACE35CE.idx`
- `.cache/clangd/index/functions_implicitlipid.cpp.801D27A088B90481.idx`
- `.cache/clangd/index/functions_implicitlipid.cpp.95F8F6317607B0EA.idx`
- `.cache/clangd/index/generate_coordinates.cpp.16B2A7A3181F9215.idx`
- `.cache/clangd/index/generate_coordinates.cpp.D69DCC233D5F38E8.idx`
- `.cache/clangd/index/generate_coordinates_for_restart.cpp.3D1771392B23CC5D.idx`
- `.cache/clangd/index/generate_coordinates_for_restart.cpp.B16A2647FD8B39F6.idx`
- `.cache/clangd/index/get_distance.cpp.13BC466368CA022D.idx`
- `.cache/clangd/index/get_distance.cpp.1E4410EE36BCFAB3.idx`
- `.cache/clangd/index/get_distance_to_surface.cpp.9CE97970184AD8D8.idx`
- `.cache/clangd/index/get_distance_to_surface.cpp.C6FA80CB76043D41.idx`
- `.cache/clangd/index/get_geodesic_distance.cpp.D1B19D2ADB17334F.idx`
- `.cache/clangd/index/get_geodesic_distance.cpp.F573EE9A7B78961E.idx`
- `.cache/clangd/index/get_prevNorm.cpp.721FD5AFF7C1D1E8.idx`
- `.cache/clangd/index/get_prevNorm.cpp.9F3A592FDD47A064.idx`
- `.cache/clangd/index/get_prevSurv.cpp.7E87C920C709D9E2.idx`
- `.cache/clangd/index/get_prevSurv.cpp.B8F8ABDB13308443.idx`
- `.cache/clangd/index/hasIntangibles.cpp.750DA381622252E8.idx`
- `.cache/clangd/index/hasIntangibles.cpp.C403F7E13B382785.idx`
- `.cache/clangd/index/implicitlipid_reactions.hpp.BE41CF893F3E3202.idx`
- `.cache/clangd/index/init_NboundPairs.cpp.160EDCA1D2DDACD1.idx`
- `.cache/clangd/index/init_NboundPairs.cpp.E64DA4E6B33DB78F.idx`
- `.cache/clangd/index/init_association_events.cpp.656FD99258AD8593.idx`
- `.cache/clangd/index/init_association_events.cpp.78E68E54C27C808B.idx`
- `.cache/clangd/index/init_counterCopyNums.cpp.169C9D87C9522E8D.idx`
- `.cache/clangd/index/init_counterCopyNums.cpp.369A76E22ABF6B4E.idx`
- `.cache/clangd/index/init_print_dimers.cpp.0E5397223F7F851A.idx`
- `.cache/clangd/index/init_print_dimers.cpp.7A0BB8D3B7ACF8C7.idx`
- `.cache/clangd/index/init_speciesFile.cpp.0EE2A3920D1A1EAF.idx`
- `.cache/clangd/index/init_speciesFile.cpp.CE2403885BC44A99.idx`
- `.cache/clangd/index/initialize_complex.cpp.58914BDD0A132399.idx`
- `.cache/clangd/index/initialize_complex.cpp.B12D093B13832147.idx`
- `.cache/clangd/index/initialize_molecule.cpp.3C498FF5B7FA505E.idx`
- `.cache/clangd/index/initialize_molecule.cpp.95EE1954C436E306.idx`
- `.cache/clangd/index/initialize_molecule_after_reaction.cpp.B6549DDE11D54A33.idx`
- `.cache/clangd/index/initialize_molecule_after_reaction.cpp.CC85C6F5EBB17AC5.idx`
- `.cache/clangd/index/initialize_molecule_after_transmission_reaction.cpp.B414F109F5F2F89B.idx`
- `.cache/clangd/index/initialize_molecule_after_transmission_reaction.cpp.E430319412D81C54.idx`
- `.cache/clangd/index/initialize_parameters_for_implicitlipid_model.cpp.22FBDEAA95954A37.idx`
- `.cache/clangd/index/initialize_parameters_for_implicitlipid_model.cpp.33812D28E82C97CB.idx`
- `.cache/clangd/index/initialize_states.cpp.5C441F6474FB6FFC.idx`
- `.cache/clangd/index/initialize_states.cpp.6435AFD291B3653F.idx`
- `.cache/clangd/index/integrator.cpp.821750084BDEC9E4.idx`
- `.cache/clangd/index/integrator.cpp.DF0F7B03A7E90F05.idx`
- `.cache/clangd/index/io.hpp.8FE3FE838B7920FC.idx`
- `.cache/clangd/index/isReactant.cpp.3161D866610F7246.idx`
- `.cache/clangd/index/isReactant.cpp.568AE496327CB29F.idx`
- `.cache/clangd/index/math_functions.cpp.90EFF53D13A8D1A8.idx`
- `.cache/clangd/index/math_functions.cpp.DE2E8C9B8E4DC3CD.idx`
- `.cache/clangd/index/math_functions.hpp.1B4F6843444C0DD4.idx`
- `.cache/clangd/index/matrix.hpp.2CBD37937179E515.idx`
- `.cache/clangd/index/matrix_functions.cpp.96DE76A2160DF10B.idx`
- `.cache/clangd/index/matrix_functions.cpp.FAC2AAE0B7737A30.idx`
- `.cache/clangd/index/measure_complex_displacement.cpp.0BB1F795AF1055EB.idx`
- `.cache/clangd/index/measure_complex_displacement.cpp.29EEDE0A9A6A3454.idx`
- `.cache/clangd/index/measure_overlap_free_protein_interfaces.cpp.3746C0BE56149FC5.idx`
- `.cache/clangd/index/measure_overlap_free_protein_interfaces.cpp.F84E6B2A389C4918.idx`
- `.cache/clangd/index/measure_overlap_protein_interfaces.cpp.96FEED36AEB5BF51.idx`
- `.cache/clangd/index/measure_overlap_protein_interfaces.cpp.AC627E953452A87D.idx`
- `.cache/clangd/index/nerdss.cpp.80B63EE2100376EF.idx`
- `.cache/clangd/index/nerdss.cpp.98F4FD2E6C3E5F90.idx`
- `.cache/clangd/index/norm_function.cpp.672A0510FB9DB49C.idx`
- `.cache/clangd/index/norm_function.cpp.CC096EA8EE978C74.idx`
- `.cache/clangd/index/omega_rotation.cpp.0D4288A3494E3C5F.idx`
- `.cache/clangd/index/omega_rotation.cpp.20880CC9A02FB952.idx`
- `.cache/clangd/index/orient_crds_to_template.cpp.41AF77913A823B96.idx`
- `.cache/clangd/index/orient_crds_to_template.cpp.994CA3097D6421EF.idx`
- `.cache/clangd/index/parse_command.cpp.25146AACF03FAAA3.idx`
- `.cache/clangd/index/parse_command.cpp.3E2D6E42F9FAD93D.idx`
- `.cache/clangd/index/parse_input.cpp.B00710BD2C81E1F1.idx`
- `.cache/clangd/index/parse_input.cpp.BD0494BEFFE0F147.idx`
- `.cache/clangd/index/parse_input_array.cpp.112811B41C24C6BD.idx`
- `.cache/clangd/index/parse_input_array.cpp.87A69827A61A915D.idx`
- `.cache/clangd/index/parse_molFile.cpp.2EA7734D8B936620.idx`
- `.cache/clangd/index/parse_molFile.cpp.8695447D0075BDCE.idx`
- `.cache/clangd/index/parse_molecule_bngl.cpp.78E8E3AC5CB39E26.idx`
- `.cache/clangd/index/parse_molecule_bngl.cpp.9BD78A04AADDC2A4.idx`
- `.cache/clangd/index/parse_observable.cpp.3015F4BE8974906F.idx`
- `.cache/clangd/index/parse_observable.cpp.69451BB6793785EB.idx`
- `.cache/clangd/index/parse_reaction.cpp.01B52812B8D37DB4.idx`
- `.cache/clangd/index/parse_reaction.cpp.341EB28634FFCD91.idx`
- `.cache/clangd/index/parse_states.cpp.5A62E9149BB859DF.idx`
- `.cache/clangd/index/parse_states.cpp.B9CB12992CF4A0A5.idx`
- `.cache/clangd/index/parser_functions.hpp.C6AC7DDAB8B9AD5F.idx`
- `.cache/clangd/index/passocF.cpp.5B71984DFDBBFDBF.idx`
- `.cache/clangd/index/passocF.cpp.AAFBFD0C5EBDE062.idx`
- `.cache/clangd/index/passocF_1D.cpp.77B4214CF419885C.idx`
- `.cache/clangd/index/perform_bimolecular_state_change.cpp.333EFC1BDC163D63.idx`
- `.cache/clangd/index/perform_bimolecular_state_change.cpp.41C1D83B440B0C8A.idx`
- `.cache/clangd/index/perform_bimolecular_state_change_box.cpp.CE3CEDF3129DFF1D.idx`
- `.cache/clangd/index/perform_bimolecular_state_change_box.cpp.F0E282660FA3685C.idx`
- `.cache/clangd/index/perform_bimolecular_state_change_sphere.cpp.AF9E179F35BC8C0A.idx`
- `.cache/clangd/index/perform_bimolecular_state_change_sphere.cpp.D7FCA77491AFB616.idx`
- `.cache/clangd/index/perform_implicitlipid_state_change.cpp.309F94221A9D606D.idx`
- `.cache/clangd/index/perform_implicitlipid_state_change.cpp.3898CBB324762504.idx`
- `.cache/clangd/index/perform_implicitlipid_state_change_box.cpp.21ABAEBCFC800E01.idx`
- `.cache/clangd/index/perform_implicitlipid_state_change_box.cpp.C8077F00E11A2386.idx`
- `.cache/clangd/index/perform_implicitlipid_state_change_sphere.cpp.4C9D562D02BB33E1.idx`
- `.cache/clangd/index/perform_implicitlipid_state_change_sphere.cpp.B5FBA8045B7A81C3.idx`
- `.cache/clangd/index/perform_transmission_reaction.cpp.8E79A0B3C8027EB3.idx`
- `.cache/clangd/index/perform_transmission_reaction.cpp.D92311C4B8220054.idx`
- `.cache/clangd/index/phi_rotation.cpp.DC11EAC0CA07F981.idx`
- `.cache/clangd/index/phi_rotation.cpp.FD90D2E6A30FC2A1.idx`
- `.cache/clangd/index/pir_function.cpp.91674D199853ADDB.idx`
- `.cache/clangd/index/pir_function.cpp.E4F055D38B805296.idx`
- `.cache/clangd/index/pirr_pfree_ratio_psF.cpp.4EC38F841FF03218.idx`
- `.cache/clangd/index/pirr_pfree_ratio_psF.cpp.A691C31BAB345CE5.idx`
- `.cache/clangd/index/pirr_pfree_ratio_psF_1D.cpp.A68F9C0908B44A34.idx`
- `.cache/clangd/index/populate_reaction_lists.cpp.304E85BC5473A9A8.idx`
- `.cache/clangd/index/populate_reaction_lists.cpp.5D9D75E42416741B.idx`
- `.cache/clangd/index/print_association_events.cpp.7A5E4282ACA3BEB9.idx`
- `.cache/clangd/index/print_association_events.cpp.DF9D99C5A6413683.idx`
- `.cache/clangd/index/print_complex_hist.cpp.B48E54801E501868.idx`
- `.cache/clangd/index/print_complex_hist.cpp.C58EEE4813E17764.idx`
- `.cache/clangd/index/print_dimers.cpp.0914723F53EEE3E4.idx`
- `.cache/clangd/index/print_dimers.cpp.A2D1B639045737BB.idx`
- `.cache/clangd/index/print_system_information.cpp.1232EF4EE6F75EC5.idx`
- `.cache/clangd/index/print_system_information.cpp.958D24918492AB2B.idx`
- `.cache/clangd/index/rand_gsl.cpp.42664C848C74200A.idx`
- `.cache/clangd/index/rand_gsl.cpp.80B80CBE7885345A.idx`
- `.cache/clangd/index/rand_gsl.hpp.ACEB269AA9B0E2AB.idx`
- `.cache/clangd/index/read_bonds.cpp.2DEAF93439058A98.idx`
- `.cache/clangd/index/read_bonds.cpp.5956AAD7DA6DD01D.idx`
- `.cache/clangd/index/read_boolean.cpp.3EE6F4D4190E4E3D.idx`
- `.cache/clangd/index/read_boolean.cpp.84083415C1F5909C.idx`
- `.cache/clangd/index/read_internal_coordinates.cpp.1E850579D8F6DFF2.idx`
- `.cache/clangd/index/read_internal_coordinates.cpp.8F6468D3BEDDFF0D.idx`
- `.cache/clangd/index/read_restart.cpp.1EF5DBBFD7224790.idx`
- `.cache/clangd/index/read_restart.cpp.EF13FDB3E22CAD27.idx`
- `.cache/clangd/index/reflect_complex_compartment.cpp.1224F15970C63878.idx`
- `.cache/clangd/index/reflect_complex_compartment.cpp.8E39336CAD854B7A.idx`
- `.cache/clangd/index/reflect_complex_rad_rot.cpp.83470A1201D5DF43.idx`
- `.cache/clangd/index/reflect_complex_rad_rot.cpp.ADDFECB32A8D664A.idx`
- `.cache/clangd/index/reflect_complex_rad_rot_box.cpp.66DEFA766CE9D1F2.idx`
- `.cache/clangd/index/reflect_complex_rad_rot_box.cpp.EB1E345E25642D52.idx`
- `.cache/clangd/index/reflect_complex_rad_rot_sphere.cpp.CCDAA21419BDEB14.idx`
- `.cache/clangd/index/reflect_complex_rad_rot_sphere.cpp.DB0B097E1535342A.idx`
- `.cache/clangd/index/reflect_functions.hpp.F43662F0B421E913.idx`
- `.cache/clangd/index/reflect_traj_check_span.cpp.2CB1A01A393740C1.idx`
- `.cache/clangd/index/reflect_traj_check_span.cpp.4E7126A71A38F442.idx`
- `.cache/clangd/index/reflect_traj_check_span_box.cpp.A62E3DFA1A2699F8.idx`
- `.cache/clangd/index/reflect_traj_check_span_box.cpp.F75D584F03D368D1.idx`
- `.cache/clangd/index/reflect_traj_check_span_sphere.cpp.3E51FA003DA5CBDF.idx`
- `.cache/clangd/index/reflect_traj_check_span_sphere.cpp.DC74740111A9922B.idx`
- `.cache/clangd/index/reflect_traj_complex_compartment.cpp.CA24B2A6CD199AAB.idx`
- `.cache/clangd/index/reflect_traj_complex_compartment.cpp.E191C793A1DE239C.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot.cpp.9B45E05045508237.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot.cpp.EA66AD9BBB3D5099.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_box.cpp.167D430E49C19648.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_box.cpp.C6CD9CC8DBD39373.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_nocheck.cpp.1D91C7BF93825375.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_nocheck.cpp.B26AD61B263755ED.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_box.cpp.5E56F80FE9B6BFE0.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_box.cpp.B8FB5D899D7207E6.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_sphere.cpp.346A7DF421B7A0D0.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_nocheck_sphere.cpp.579412256539EA4F.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_sphere.cpp.8C26C2DB6679076A.idx`
- `.cache/clangd/index/reflect_traj_complex_rad_rot_sphere.cpp.A999161725CAD16F.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds.cpp.300F4C7090381D20.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds.cpp.FC4AB1972B08A76C.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds_box.cpp.877F4FD732DC927A.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds_box.cpp.C5087CF17B06EF9D.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds_compartment.cpp.5F4D7E7C3A8C31C6.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds_compartment.cpp.DC9A160C71AF3260.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds_sphere.cpp.15407B68AAFDA41A.idx`
- `.cache/clangd/index/reflect_traj_tmp_crds_sphere.cpp.A0593DC7A5EB3B6E.idx`
- `.cache/clangd/index/remove_comment.cpp.52FDE4B10B1B7D5F.idx`
- `.cache/clangd/index/remove_comment.cpp.80E4A2E2672773D8.idx`
- `.cache/clangd/index/requiresSignFlip.cpp.12FA43DABADF6829.idx`
- `.cache/clangd/index/requiresSignFlip.cpp.BFA627775DD27381.idx`
- `.cache/clangd/index/reverse_rotation.cpp.2357217350D1937F.idx`
- `.cache/clangd/index/reverse_rotation.cpp.FF87CFB3C78B8A6F.idx`
- `.cache/clangd/index/rotate.cpp.93616F7AAD2ED0DF.idx`
- `.cache/clangd/index/rotate.cpp.FFB632E877B19079.idx`
- `.cache/clangd/index/save_mem_orientation.cpp.4318E508B6D3CCFF.idx`
- `.cache/clangd/index/save_mem_orientation.cpp.EF1D0CA7F40CF579.idx`
- `.cache/clangd/index/set_rMaxLimit.cpp.6817D3F55C8D917B.idx`
- `.cache/clangd/index/set_rMaxLimit.cpp.EF935BE6CE9C368D.idx`
- `.cache/clangd/index/shared_reaction_functions.hpp.78F5129E8A21E9AE.idx`
- `.cache/clangd/index/size_lookup.cpp.3D203DD11EA0F57F.idx`
- `.cache/clangd/index/size_lookup.cpp.F5E59B74A1B92311.idx`
- `.cache/clangd/index/subtract_off_com_position.cpp.51D15F5904BDEBCA.idx`
- `.cache/clangd/index/subtract_off_com_position.cpp.74131E55141860B3.idx`
- `.cache/clangd/index/survival_function.cpp.137F96F31D4F39B5.idx`
- `.cache/clangd/index/survival_function.cpp.61AD740FDA38C84E.idx`
- `.cache/clangd/index/sweep_separation_complex_rot.cpp.138C679ECC97BD67.idx`
- `.cache/clangd/index/sweep_separation_complex_rot.cpp.2FDE65E038AA1F8D.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_box.cpp.AB5601CE2E1F47F6.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_box.cpp.FB720D7F93035D3B.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest.cpp.2AC6B4BD97E8D867.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest.cpp.8DD995669F88DE3A.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_box.cpp.67AEB26631344EF9.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_box.cpp.84BB62414E7A9EB1.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_cluster.cpp.6AA1F8C1C98421CB.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_cluster.cpp.FF8159554073BAEC.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_box.cpp.D34CD4E3BFF1A907.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_box.cpp.D91799B387F49E2B.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_sphere.cpp.7CB8285DB34C9BC3.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_cluster_sphere.cpp.C096DD9B75FEF2E6.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_sphere.cpp.3544A93E7245AA23.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_memtest_sphere.cpp.7533B7D5FF46A997.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_sphere.cpp.23E5A7C15BAD5C20.idx`
- `.cache/clangd/index/sweep_separation_complex_rot_sphere.cpp.E96C225D7E451FB1.idx`
- `.cache/clangd/index/system_setup.hpp.975D395A733EAA9B.idx`
- `.cache/clangd/index/theta_rotation.cpp.725C02B950F89CBB.idx`
- `.cache/clangd/index/theta_rotation.cpp.E7F8FCBAAF5EDD25.idx`
- `.cache/clangd/index/tracing.hpp.1E86E9BD1E155318.idx`
- `.cache/clangd/index/track_association_events.cpp.5EC1BB0518EF4C88.idx`
- `.cache/clangd/index/track_association_events.cpp.BC9C3A5899771857.idx`
- `.cache/clangd/index/trajectory_functions.hpp.AC79A5E4DB6676EE.idx`
- `.cache/clangd/index/transform.cpp.0799BA799CB62539.idx`
- `.cache/clangd/index/transform.cpp.A294C2E06FEAA8EC.idx`
- `.cache/clangd/index/unimolecular_reactions.hpp.E68FB70E56D09C26.idx`
- `.cache/clangd/index/update_Nboundpairs.cpp.2F87FB400C0CBE76.idx`
- `.cache/clangd/index/update_Nboundpairs.cpp.60C03811F69BE3CC.idx`
- `.cache/clangd/index/update_complex_tmp_com_crds.cpp.0019E67F368F77F2.idx`
- `.cache/clangd/index/update_complex_tmp_com_crds.cpp.B455D89D8F568DB2.idx`
- `.cache/clangd/index/write_NboundPairs.cpp.24ECA0CF0D20459B.idx`
- `.cache/clangd/index/write_NboundPairs.cpp.AC875BA9E909A18F.idx`
- `.cache/clangd/index/write_all_species.cpp.29EB2818CDE699E7.idx`
- `.cache/clangd/index/write_all_species.cpp.5CEA5725D26E7FE7.idx`
- `.cache/clangd/index/write_complex_components.cpp.7C6E3BF548909F1B.idx`
- `.cache/clangd/index/write_complex_components.cpp.E10F59063DA2209C.idx`
- `.cache/clangd/index/write_complex_crds.cpp.5D2BB1B5982AB2E4.idx`
- `.cache/clangd/index/write_complex_crds.cpp.F956FAA0CC551A0D.idx`
- `.cache/clangd/index/write_crds.cpp.87B153880311ECEC.idx`
- `.cache/clangd/index/write_crds.cpp.8EC8485290ED5644.idx`
- `.cache/clangd/index/write_mol_iface.cpp.A132C69B74F5FBA4.idx`
- `.cache/clangd/index/write_mol_iface.cpp.A70EADBF89E29E06.idx`
- `.cache/clangd/index/write_observables.cpp.215B012073FF5ADC.idx`
- `.cache/clangd/index/write_observables.cpp.DFEFE8F8472E1B98.idx`
- `.cache/clangd/index/write_pdb.cpp.7AE30E5B95F74C2D.idx`
- `.cache/clangd/index/write_pdb.cpp.E8E9139B4AFC19F2.idx`
- `.cache/clangd/index/write_psf.cpp.24855AE09A3509C7.idx`
- `.cache/clangd/index/write_psf.cpp.A63F54587DF31E90.idx`
- `.cache/clangd/index/write_restart.cpp.80E381FC88A4FACC.idx`
- `.cache/clangd/index/write_restart.cpp.D1EF0F3CE1EE6D2A.idx`
- `.cache/clangd/index/write_timestep_information.cpp.217117A0DC33022D.idx`
- `.cache/clangd/index/write_timestep_information.cpp.31443B470577CB15.idx`
- `.cache/clangd/index/write_traj.cpp.0CE11FFD455FCD3A.idx`
- `.cache/clangd/index/write_traj.cpp.260F4AC3A3565830.idx`
- `.cache/clangd/index/write_transition.cpp.1BC5F68470F3D951.idx`
- `.cache/clangd/index/write_transition.cpp.1BF6ED8E56F75945.idx`
- `.cache/clangd/index/write_xyz.cpp.2C4B77EE4DC45A9D.idx`
- `.cache/clangd/index/write_xyz.cpp.AF1F327DB8F0F77A.idx`
- `.cache/clangd/index/write_xyz_assoc.cpp.0F06761FFD82FEB5.idx`
- `.cache/clangd/index/write_xyz_assoc.cpp.0FEBE136BE2B9E57.idx`
- `.cache/clangd/index/write_xyz_assoc_cout.cpp.2435F8AAB24E9634.idx`
- `.cache/clangd/index/write_xyz_assoc_cout.cpp.2C6D704A7AF33BA3.idx`
- `.github/workflows/c-cpp.yml`

## Deleted Files

- `obj/trajectory_functions/sweep_separation_complex_rot_sphere.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_memtest_sphere.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_memtest_cluster_sphere.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_memtest_cluster_box.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_memtest_cluster.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_memtest_box.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_memtest.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_fiber_box.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_fiber.o`
- `obj/trajectory_functions/sweep_separation_complex_rot_box.o`
- `obj/trajectory_functions/sweep_separation_complex_rot.o`
- `obj/trajectory_functions/create_complex_propagation_vectors_on_sphere.o`
- `obj/trajectory_functions/create_complex_propagation_vectors.o`
- `obj/trajectory_functions/clear_reweight_vecs.o`
- `obj/trajectory_functions/calculate_update_position_interface.o`
- `obj/system_setup/set_rMaxLimit.o`
- `obj/system_setup/set_exclude_volume_bound.o`
- `obj/system_setup/initialize_states.o`
- `obj/system_setup/initialize_parameters_for_implicitlipid_and_compartment_model.o`
- `obj/system_setup/initialize_molecule.o`
- `obj/system_setup/initialize_complex.o`
- `obj/system_setup/generate_coordinates_for_restart.o`
- `obj/system_setup/generate_coordinates.o`
- `obj/system_setup/determine_shape_molecule.o`
- `obj/system_setup/are_molecules_in_vicinity.o`
- `obj/system_setup/areInVicinity.o`
- `obj/reactions/update_complex_tmp_com_crds.o`
- `obj/reactions/update_Nboundpairs.o`
- `obj/reactions/transform.o`
- `obj/reactions/track_association_events.o`
- `obj/reactions/theta_rotation.o`
- `obj/reactions/survival_function.o`
- `obj/reactions/subtract_off_com_position.o`
- `obj/reactions/size_lookup.o`
- `obj/reactions/save_mem_orientation.o`
- `obj/reactions/rotate.o`
- `obj/reactions/reverse_rotation.o`
- `obj/reactions/requiresSignFlip.o`
- `obj/reactions/remove_empty_slots.o`
- `obj/reactions/print_association_events.o`
- `obj/reactions/pirr_pfree_ratio_psF_1D.o`
- `obj/reactions/pirr_pfree_ratio_psF.o`
- `obj/reactions/pir_function.o`
- `obj/reactions/phi_rotation.o`
- `obj/reactions/perform_transmission_reaction.o`
- `obj/reactions/perform_implicitlipid_state_change_sphere.o`
- `obj/reactions/perform_implicitlipid_state_change_box.o`
- `obj/reactions/perform_implicitlipid_state_change.o`
- `obj/reactions/perform_bimolecular_state_change_sphere.o`
- `obj/reactions/perform_bimolecular_state_change_box.o`
- `obj/reactions/perform_bimolecular_state_change.o`
- `obj/reactions/perform_bimolecular_reactions.o`
- `obj/reactions/passocF_1D.o`
- `obj/reactions/passocF.o`
- `obj/reactions/orient_crds_to_template.o`
- `obj/reactions/omega_rotation.o`
- `obj/reactions/norm_function.o`
- `obj/reactions/measure_separations_to_identify_possible_reactions.o`
- `obj/reactions/measure_overlap_protein_interfaces.o`
- `obj/reactions/measure_overlap_free_protein_interfaces.o`
- `obj/reactions/measure_complex_displacement.o`
- `obj/reactions/isReactant.o`
- `obj/reactions/integrator.o`
- `obj/reactions/initialize_molecule_after_transmission_reaction.o`
- `obj/reactions/initialize_molecule_after_reaction.o`
- `obj/reactions/init_association_events.o`
- `obj/reactions/hasIntangibles.o`
- `obj/reactions/get_prevSurv.o`
- `obj/reactions/get_prevNorm.o`
- `obj/reactions/get_geodesic_distance.o`
- `obj/reactions/get_distance_to_surface.o`
- `obj/reactions/get_distance.o`
- `obj/reactions/functions_implicitlipid.o`
- `obj/reactions/functions_for_spherical_system.o`
- `obj/reactions/find_which_state_change_reaction.o`
- `obj/reactions/find_which_reaction.o`
- `obj/reactions/find_reaction_rate_state.o`
- `obj/reactions/evaluate_binding_within_complex.o`
- `obj/reactions/determine_rotation_angles.o`
- `obj/reactions/determine_parent_complex_IL.o`
- `obj/reactions/determine_parent_complex.o`
- `obj/reactions/determine_normal.o`
- `obj/reactions/determine_if_reaction_occurs.o`
- `obj/reactions/determine_exiting_compartment_probability.o`
- `obj/reactions/determine_entering_compartment_probability.o`
- `obj/reactions/determine_3D_implicitlipid_reaction_probability.o`
- `obj/reactions/determine_3D_bimolecular_reaction_probability.o`
- `obj/reactions/determine_2D_implicitlipid_reaction_probability.o`
- `obj/reactions/determine_2D_bimolecular_reaction_probability.o`
- `obj/reactions/determine_1D_bimolecular_reaction_probability.o`
- `obj/reactions/create_survMatrix.o`
- `obj/reactions/create_pirMatrix.o`
- `obj/reactions/create_normMatrix.o`
- `obj/reactions/create_molecule_and_complex_from_transmission_rxn.o`
- `obj/reactions/create_molecule_and_complex_from_rxn.o`
- `obj/reactions/create_arbitrary_vector.o`
- `obj/reactions/create_DDMatrices.o`
- `obj/reactions/correct_structutre.o`
- `obj/reactions/conservedRigid.o`
- `obj/reactions/conservedMags.o`
- `obj/reactions/com_of_two_tmp_complexes.o`
- `obj/reactions/check_perform_zeroth_first_order_reactions.o`
- `obj/reactions/check_overlap.o`
- `obj/reactions/check_implicit_reactions.o`
- `obj/reactions/check_for_zeroth_order_creation.o`
- `obj/reactions/check_for_unimolstatechange_reactions.o`
- `obj/reactions/check_for_unimolecular_reactions_population.o`
- `obj/reactions/check_for_unimolecular_reactions.o`
- `obj/reactions/check_for_structure_overlap_system.o`
- `obj/reactions/check_for_structure_overlap.o`
- `obj/reactions/check_dissociation_implicitlipid.o`
- `obj/reactions/check_dissociation.o`
- `obj/reactions/check_compartment_reaction.o`
- `obj/reactions/check_bimolecular_reactions.o`
- `obj/reactions/check_bases.o`
- `obj/reactions/calculate_phi.o`
- `obj/reactions/calculate_omega.o`
- `obj/reactions/calc_pirr.o`
- `obj/reactions/calc_one_angular_displacement.o`
- `obj/reactions/calc_angular_displacement.o`
- `obj/reactions/break_interaction_implicitlipid.o`
- `obj/reactions/break_interaction.o`
- `obj/reactions/associate_sphere.o`
- `obj/reactions/associate_box.o`
- `obj/reactions/associate_ImplicitLipid_sphere.o`
- `obj/reactions/associate_ImplicitLipid_box.o`
- `obj/reactions/associate_ImplicitLipid.o`
- `obj/reactions/associate.o`
- `obj/reactions/areSameAngle.o`
- `obj/reactions/areParallel.o`
- `obj/reactions/angleSignIsCorrect.o`
- `obj/reactions/DDpirr_pfree_ratio_ps.o`
- `obj/parser/write_mol_iface.o`
- `obj/parser/remove_comment.o`
- `obj/parser/read_internal_coordinates.o`
- `obj/parser/read_boolean.o`
- `obj/parser/read_bonds.o`
- `obj/parser/populate_reaction_lists.o`
- `obj/parser/parse_states.o`
- `obj/parser/parse_reaction.o`
- `obj/parser/parse_observable.o`
- `obj/parser/parse_molecule_bngl.o`
- `obj/parser/parse_molFile.o`
- `obj/parser/parse_input_for_add_file.o`
- `obj/parser/parse_input_for_a_restart_simulation.o`
- `obj/parser/parse_input_for_a_new_simulation.o`
- `obj/parser/parse_input_array.o`
- `obj/parser/parse_input.o`
- `obj/parser/parse_command.o`
- `obj/parser/display_all_reactions.o`
- `obj/parser/display_all_MolTemplates.o`
- `obj/parser/determine_iface_indices.o`
- `obj/parser/determine_bound_iface_index.o`
- `obj/parser/create_conjugate_reaction_itrs.o`
- `obj/parser/check_for_valid_states.o`
- `obj/parser/check_for_state_change.o`
- `obj/parser/areSameExceptState.o`
- `obj/math/rand_gsl.o`
- `obj/math/matrix_functions.o`
- `obj/math/math_functions.o`
- `obj/math/Faddeeva.o`
- `obj/io/write_xyz_assoc_cout.o`
- `obj/io/write_xyz_assoc.o`
- `obj/io/write_xyz.o`
- `obj/io/write_transition.o`
- `obj/io/write_traj.o`
- `obj/io/write_timestep_information.o`
- `obj/io/write_restart.o`
- `obj/io/write_psf.o`
- `obj/io/write_pdb.o`
- `obj/io/write_observables.o`
- `obj/io/write_crds.o`
- `obj/io/write_complex_crds.o`
- `obj/io/write_complex_components.o`
- `obj/io/write_bonded_complex_json.o`
- `obj/io/write_all_species.o`
- `obj/io/write_NboundPairs.o`
- `obj/io/read_restart.o`
- `obj/io/print_system_information.o`
- `obj/io/print_dimers.o`
- `obj/io/print_complex_hist.o`
- `obj/io/init_speciesFile.o`
- `obj/io/init_print_dimers.o`
- `obj/io/init_counterCopyNums.o`
- `obj/io/init_NboundPairs.o`
- `obj/io/eye_candy.o`
- `obj/error/error.o`
- `obj/classes/mpi_functions.o`
- `obj/classes/class_bngl_parser_functions.o`
- `obj/classes/class_Vector.o`
- `obj/classes/class_SimulVolume.o`
- `obj/classes/class_Rxns.o`
- `obj/classes/class_Quat.o`
- `obj/classes/class_Parameters.o`
- `obj/classes/class_Observable.o`
- `obj/classes/class_Molecule_Complex.o`
- `obj/classes/class_MolTemplate.o`
- `obj/classes/class_MDTimer.o`
- `obj/classes/class_Coord.o`
- `obj/classes/class_Cluster.o`
- `obj/boundary_conditions/reflect_traj_tmp_crds_sphere.o`
- `obj/boundary_conditions/reflect_traj_tmp_crds_compartment.o`
- `obj/boundary_conditions/reflect_traj_tmp_crds_box.o`
- `obj/boundary_conditions/reflect_traj_tmp_crds.o`
- `obj/boundary_conditions/reflect_traj_complex_rad_rot_sphere.o`
- `obj/boundary_conditions/reflect_traj_complex_rad_rot_nocheck_sphere.o`
- `obj/boundary_conditions/reflect_traj_complex_rad_rot_nocheck_box.o`
- `obj/boundary_conditions/reflect_traj_complex_rad_rot_nocheck.o`
- `obj/boundary_conditions/reflect_traj_complex_rad_rot_box.o`
- `obj/boundary_conditions/reflect_traj_complex_rad_rot.o`
- `obj/boundary_conditions/reflect_traj_complex_compartment.o`
- `obj/boundary_conditions/reflect_traj_check_span_sphere.o`
- `obj/boundary_conditions/reflect_traj_check_span_box.o`
- `obj/boundary_conditions/reflect_traj_check_span.o`
- `obj/boundary_conditions/reflect_complex_rad_rot_sphere.o`
- `obj/boundary_conditions/reflect_complex_rad_rot_box.o`
- `obj/boundary_conditions/reflect_complex_rad_rot.o`
- `obj/boundary_conditions/reflect_complex_compartment.o`
- `obj/boundary_conditions/check_if_spans_sphere.o`
- `obj/boundary_conditions/check_if_spans_box.o`
- `obj/boundary_conditions/check_if_spans.o`
- `images/.DS_Store`
- `bin/nerdss`
- `assets/.DS_Store`
- `PDB/99000.pdb`
- `PDB/98000.pdb`
- `PDB/97000.pdb`
- `PDB/96000.pdb`
- `PDB/95000.pdb`
- `PDB/94000.pdb`
- `PDB/93000.pdb`
- `PDB/92000.pdb`
- `PDB/91000.pdb`
- `PDB/90000.pdb`
- `PDB/9000.pdb`
- `PDB/89000.pdb`
- `PDB/88000.pdb`
- `PDB/87000.pdb`
- `PDB/86000.pdb`
- `PDB/85000.pdb`
- `PDB/84000.pdb`
- `PDB/83000.pdb`
- `PDB/82000.pdb`
- `PDB/81000.pdb`
- `PDB/80000.pdb`
- `PDB/8000.pdb`
- `PDB/79000.pdb`
- `PDB/78000.pdb`
- `PDB/77000.pdb`
- `PDB/76000.pdb`
- `PDB/75000.pdb`
- `PDB/74000.pdb`
- `PDB/73000.pdb`
- `PDB/72000.pdb`
- `PDB/71000.pdb`
- `PDB/70000.pdb`
- `PDB/7000.pdb`
- `PDB/69000.pdb`
- `PDB/68000.pdb`
- `PDB/67000.pdb`
- `PDB/66000.pdb`
- `PDB/65000.pdb`
- `PDB/64000.pdb`
- `PDB/63000.pdb`
- `PDB/62000.pdb`
- `PDB/61000.pdb`
- `PDB/60000.pdb`
- `PDB/6000.pdb`
- `PDB/59000.pdb`
- `PDB/58000.pdb`
- `PDB/57000.pdb`
- `PDB/56000.pdb`
- `PDB/55000.pdb`
- `PDB/54000.pdb`
- `PDB/53000.pdb`
- `PDB/52000.pdb`
- `PDB/51000.pdb`
- `PDB/50000.pdb`
- `PDB/5000.pdb`
- `PDB/49000.pdb`
- `PDB/48000.pdb`
- `PDB/47000.pdb`
- `PDB/46000.pdb`
- `PDB/45000.pdb`
- `PDB/44000.pdb`
- `PDB/43000.pdb`
- `PDB/42000.pdb`
- `PDB/41000.pdb`
- `PDB/40000.pdb`
- `PDB/4000.pdb`
- `PDB/39000.pdb`
- `PDB/38000.pdb`
- `PDB/37000.pdb`
- `PDB/36000.pdb`
- `PDB/35000.pdb`
- `PDB/34000.pdb`
- `PDB/33000.pdb`
- `PDB/32000.pdb`
- `PDB/31000.pdb`
- `PDB/30000.pdb`
- `PDB/3000.pdb`
- `PDB/29000.pdb`
- `PDB/289000.pdb`
- `PDB/288000.pdb`
- `PDB/287000.pdb`
- `PDB/286000.pdb`
- `PDB/285000.pdb`
- `PDB/284000.pdb`
- `PDB/283000.pdb`
- `PDB/282000.pdb`
- `PDB/281000.pdb`
- `PDB/280000.pdb`
- `PDB/28000.pdb`
- `PDB/279000.pdb`
- `PDB/278000.pdb`
- `PDB/277000.pdb`
- `PDB/276000.pdb`
- `PDB/275000.pdb`
- `PDB/274000.pdb`
- `PDB/273000.pdb`
- `PDB/272000.pdb`
- `PDB/271000.pdb`
- `PDB/270000.pdb`
- `PDB/27000.pdb`
- `PDB/269000.pdb`
- `PDB/268000.pdb`
- `PDB/267000.pdb`
- `PDB/266000.pdb`
- `PDB/265000.pdb`
- `PDB/264000.pdb`
- `PDB/263000.pdb`
- `PDB/262000.pdb`
- `PDB/261000.pdb`
- `PDB/260000.pdb`
- `PDB/26000.pdb`
- `PDB/259000.pdb`
- `PDB/258000.pdb`
- `PDB/257000.pdb`
- `PDB/256000.pdb`
- `PDB/255000.pdb`
- `PDB/254000.pdb`
- `PDB/253000.pdb`
- `PDB/252000.pdb`
- `PDB/251000.pdb`
- `PDB/250000.pdb`
- `PDB/25000.pdb`
- `PDB/249000.pdb`
- `PDB/248000.pdb`
- `PDB/247000.pdb`
- `PDB/246000.pdb`
- `PDB/245000.pdb`
- `PDB/244000.pdb`
- `PDB/243000.pdb`
- `PDB/242000.pdb`
- `PDB/241000.pdb`
- `PDB/240000.pdb`
- `PDB/24000.pdb`
- `PDB/239000.pdb`
- `PDB/238000.pdb`
- `PDB/237000.pdb`
- `PDB/236000.pdb`
- `PDB/235000.pdb`
- `PDB/234000.pdb`
- `PDB/233000.pdb`
- `PDB/232000.pdb`
- `PDB/231000.pdb`
- `PDB/230000.pdb`
- `PDB/23000.pdb`
- `PDB/229000.pdb`
- `PDB/228000.pdb`
- `PDB/227000.pdb`
- `PDB/226000.pdb`
- `PDB/225000.pdb`
- `PDB/224000.pdb`
- `PDB/223000.pdb`
- `PDB/222000.pdb`
- `PDB/221000.pdb`
- `PDB/220000.pdb`
- `PDB/22000.pdb`
- `PDB/219000.pdb`
- `PDB/218000.pdb`
- `PDB/217000.pdb`
- `PDB/216000.pdb`
- `PDB/215000.pdb`
- `PDB/214000.pdb`
- `PDB/213000.pdb`
- `PDB/212000.pdb`
- `PDB/211000.pdb`
- `PDB/210000.pdb`
- `PDB/21000.pdb`
- `PDB/209000.pdb`
- `PDB/208000.pdb`
- `PDB/207000.pdb`
- `PDB/206000.pdb`
- `PDB/205000.pdb`
- `PDB/204000.pdb`
- `PDB/203000.pdb`
- `PDB/202000.pdb`
- `PDB/201000.pdb`
- `PDB/200000.pdb`
- `PDB/20000.pdb`
- `PDB/2000.pdb`
- `PDB/199000.pdb`
- `PDB/198000.pdb`
- `PDB/197000.pdb`
- `PDB/196000.pdb`
- `PDB/195000.pdb`
- `PDB/194000.pdb`
- `PDB/193000.pdb`
- `PDB/192000.pdb`
- `PDB/191000.pdb`
- `PDB/190000.pdb`
- `PDB/19000.pdb`
- `PDB/189000.pdb`
- `PDB/188000.pdb`
- `PDB/187000.pdb`
- `PDB/186000.pdb`
- `PDB/185000.pdb`
- `PDB/184000.pdb`
- `PDB/183000.pdb`
- `PDB/182000.pdb`
- `PDB/181000.pdb`
- `PDB/180000.pdb`
- `PDB/18000.pdb`
- `PDB/179000.pdb`
- `PDB/178000.pdb`
- `PDB/177000.pdb`
- `PDB/176000.pdb`
- `PDB/175000.pdb`
- `PDB/174000.pdb`
- `PDB/173000.pdb`
- `PDB/172000.pdb`
- `PDB/171000.pdb`
- `PDB/170000.pdb`
- `PDB/17000.pdb`
- `PDB/169000.pdb`
- `PDB/168000.pdb`
- `PDB/167000.pdb`
- `PDB/166000.pdb`
- `PDB/165000.pdb`
- `PDB/164000.pdb`
- `PDB/163000.pdb`
- `PDB/162000.pdb`
- `PDB/161000.pdb`
- `PDB/160000.pdb`
- `PDB/16000.pdb`
- `PDB/159000.pdb`
- `PDB/158000.pdb`
- `PDB/157000.pdb`
- `PDB/156000.pdb`
- `PDB/155000.pdb`
- `PDB/154000.pdb`
- `PDB/153000.pdb`
- `PDB/152000.pdb`
- `PDB/151000.pdb`
- `PDB/150000.pdb`
- `PDB/15000.pdb`
- `PDB/149000.pdb`
- `PDB/148000.pdb`
- `PDB/147000.pdb`
- `PDB/146000.pdb`
- `PDB/145000.pdb`
- `PDB/144000.pdb`
- `PDB/143000.pdb`
- `PDB/142000.pdb`
- `PDB/141000.pdb`
- `PDB/140000.pdb`
- `PDB/14000.pdb`
- `PDB/139000.pdb`
- `PDB/138000.pdb`
- `PDB/137000.pdb`
- `PDB/136000.pdb`
- `PDB/135000.pdb`
- `PDB/134000.pdb`
- `PDB/133000.pdb`
- `PDB/132000.pdb`
- `PDB/131000.pdb`
- `PDB/130000.pdb`
- `PDB/13000.pdb`
- `PDB/129000.pdb`
- `PDB/128000.pdb`
- `PDB/127000.pdb`
- `PDB/126000.pdb`
- `PDB/125000.pdb`
- `PDB/124000.pdb`
- `PDB/123000.pdb`
- `PDB/122000.pdb`
- `PDB/121000.pdb`
- `PDB/120000.pdb`
- `PDB/12000.pdb`
- `PDB/119000.pdb`
- `PDB/118000.pdb`
- `PDB/117000.pdb`
- `PDB/116000.pdb`
- `PDB/115000.pdb`
- `PDB/114000.pdb`
- `PDB/113000.pdb`
- `PDB/112000.pdb`
- `PDB/111000.pdb`
- `PDB/110000.pdb`
- `PDB/11000.pdb`
- `PDB/109000.pdb`
- `PDB/108000.pdb`
- `PDB/107000.pdb`
- `PDB/106000.pdb`
- `PDB/105000.pdb`
- `PDB/104000.pdb`
- `PDB/103000.pdb`
- `PDB/102000.pdb`
- `PDB/101000.pdb`
- `PDB/100000.pdb`
- `PDB/10000.pdb`
- `PDB/1000.pdb`
- `PDB/0.pdb`
- `system.psf`
- `.DS_Store`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/RESTARTS/restart60000.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/RESTARTS/restart180000.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/RESTARTS/restart120000.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/PDB/0.pdb`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/transition_matrix_time.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/trajectory.xyz`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/system.psf`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/observables_time.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/mono_dimer_time.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/initial_crds.xyz`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/histogram_complexes_time.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/event_counters_time.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/copy_numbers_time.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/bound_pair_time.dat`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/assoc_dissoc_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart9000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart800000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart8000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart7000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart600000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart60000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart6000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart540000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart5000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart480000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart420000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart400000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart4000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart360000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart300000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart3000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart240000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart200000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart2000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1800000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart180000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1600000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1400000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1200000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart120000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1000000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/restart1000.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/9999.pdb`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/599999.pdb`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/1999999.pdb`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/1000000.pdb`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/0.pdb`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/transition_matrix_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/trajectory.xyz`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/system.psf`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/restart.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/observables_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/mono_dimer_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/initial_crds.xyz`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/histogram_complexes_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/final_coords.xyz`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/event_counters_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/copy_numbers_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/bound_pair_time.dat`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/assoc_dissoc_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/900000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/800000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/700000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/600000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/500000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4999999.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4900000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4800000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4700000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4600000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4500000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4400000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4300000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4200000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4100000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/4000000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/400000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3900000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3800000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3700000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3600000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3500000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3400000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3300000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3200000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3100000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/3000000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/300000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2900000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2800000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2700000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2600000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2500000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2400000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2300000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2200000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2100000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/2000000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/200000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1900000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1800000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1700000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1600000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1500000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1400000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1300000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1200000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1100000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/1000000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/100000.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/PDB/0.pdb`
- `sample_inputs/VALIDATE_SUITE/trimer/transition_matrix_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/trajectory.xyz`
- `sample_inputs/VALIDATE_SUITE/trimer/system.psf`
- `sample_inputs/VALIDATE_SUITE/trimer/restart.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/observables_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/mono_dimer_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/initial_crds.xyz`
- `sample_inputs/VALIDATE_SUITE/trimer/histogram_complexes_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/final_coords.xyz`
- `sample_inputs/VALIDATE_SUITE/trimer/event_counters_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/copy_numbers_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/bound_pair_time.dat`
- `sample_inputs/VALIDATE_SUITE/trimer/assoc_dissoc_time.dat`
- `src/.DS_Store`

## Content-Modified Files

- `.gitignore`
- `Makefile`
- `EXEs/nerdss.cpp`
- `include/boundary_conditions/reflect_functions.hpp`
- `include/classes/class_Molecule_Complex.hpp`
- `src/classes/class_Molecule_Complex.cpp`
- `src/classes/class_Parameters.cpp`
- `src/classes/class_Vector.cpp`
- `src/io/read_restart.cpp`
- `src/io/write_restart.cpp`
- `src/reactions/associate_box.cpp`
- `src/reactions/check_for_unimolecular_reactions_population.cpp`
- `src/reactions/determine_1D_bimolecular_reaction_probability.cpp`
- `src/reactions/pirr_pfree_ratio_psF_1D.cpp`
- `src/trajectory_functions/calculate_update_position_interface.cpp`
- `src/trajectory_functions/create_complex_propagation_vectors.cpp`

## Metadata-Only File Differences

- `CMakeLists.txt [.f..tp...]`
- `Dockerfile [.f..tp...]`
- `LICENSE [.f..tp...]`
- `Makefile_old [.f..tp...]`
- `NERDSS_USER_GUIDE_MASTER.docx [.f..tp...]`
- `NERDSS_USER_GUIDE_MASTER.pdf [.f..tp...]`
- `README.md [.f..tp...]`
- `gui.py [.f..tp...]`
- `text_search.bsh [.f..tp...]`
- `EXEs/nerdss_mpi.cpp [.f..tp...]`
- `MATLAB_functions/clath_crds.m [.f..tp...]`
- `MATLAB_functions/design_pucker_crds.m [.f..tp...]`
- `MATLAB_functions/max_timestep_nerdss.m [.f..tp...]`
- `MATLAB_functions/plot_clath_only.m [.f..tp...]`
- `MATLAB_functions/rotateEulerX.m [.f..tp...]`
- `MATLAB_functions/rotateEulerZ.m [.f..tp...]`
- `MATLAB_functions/rotation_axis_angle.m [.f..tp...]`
- `angle_calculation/PDB.py [.f..tp...]`
- `angle_calculation/README [.f..tp...]`
- `angle_calculation/angles.py [.f..tp...]`
- `angle_calculation/chain_int.py [.f..tp...]`
- `docs/Complexes.png [.f..tp...]`
- `docs/Fully_Concurrent_Model.png [.f..tp...]`
- `docs/NERDSS_Parallel_Developer_Guide.html [.f..tp...]`
- `docs/NERDSS_Parallel_Developer_Guide.md [.f..tp...]`
- `docs/Run_NERDSS_colab.ipynb [.f..tp...]`
- `docs/Tux2.png [.f..tp...]`
- `docs/angle_documentation.pdf [.f..tp...]`
- `docs/association.md [.f..tp...]`
- `docs/changelog.md [.f..tp...]`
- `docs/devel.md [.f..tp...]`
- `docs/doxygen_html_style.css [.f..tp...]`
- `docs/doxygen_info.md [.f..tp...]`
- `docs/doxygen_layout.xml [.f..tp...]`
- `docs/examples.md [.f..tp...]`
- `docs/index.html [.f..tp...]`
- `docs/input.md [.f..tp...]`
- `docs/mainpage.md [.f..tp...]`
- `docs/nerdss_banner_export.svg [.f..tp...]`
- `docs/associate_figs/algo_flowchart_horiz.pdf [.f..tp...]`
- `docs/associate_figs/algo_flowchart_horiz.tex [.f..tp...]`
- `docs/associate_figs/ap2_example.pdf [.f..tp...]`
- `docs/associate_figs/pucker_clat_example.pdf [.f..tp...]`
- `docs/associate_figs/angles/angles.pdf [.f..tp...]`
- `docs/associate_figs/angles/angles.tex [.f..tp...]`
- `docs/associate_figs/angles/assoc_before.pdf [.f..tp...]`
- `docs/associate_figs/angles/assoc_before.tex [.f..tp...]`
- `docs/associate_figs/angles/omega.pdf [.f..tp...]`
- `docs/associate_figs/angles/omega.tex [.f..tp...]`
- `docs/associate_figs/angles/phi.pdf [.f..tp...]`
- `docs/associate_figs/angles/phi.tex [.f..tp...]`
- `docs/associate_figs/angles/theta.eps [.f..tp...]`
- `docs/associate_figs/angles/theta.pdf [.f..tp...]`
- `docs/associate_figs/angles/theta.tex [.f..tp...]`
- `docs/examples/irrev_flat_clat/species.png [.f..tp...]`
- `docs/examples/irrev_flat_clat/species_kd100.png [.f..tp...]`
- `docs/examples/zero_create/zeroth_order_creation_avg.png [.f..tp...]`
- `docs/examples/zero_create/zeroth_order_creation_indiv.png [.f..tp...]`
- `docs/html/2_d__reaction__table__functions_8hpp.html [.f..tp...]`
- `docs/html/2_d__reaction__table__functions_8hpp_source.html [.f..tp...]`
- `docs/html/2d__reaction__table__functions_8cpp.html [.f..tp...]`
- `docs/html/_faddeeva_8hpp_source.html [.f..tp...]`
- `docs/html/algo_flowchart_horiz.pdf [.f..tp...]`
- `docs/html/angle_documentation.pdf [.f..tp...]`
- `docs/html/annotated.html [.f..tp...]`
- `docs/html/ap2_example.pdf [.f..tp...]`
- `docs/html/assoc_before.pdf [.f..tp...]`
- `docs/html/associate_8cpp.html [.f..tp...]`
- `docs/html/association.html [.f..tp...]`
- `docs/html/association_8hpp_source.html [.f..tp...]`
- `docs/html/association__functions_8cpp.html [.f..tp...]`
- `docs/html/association__functions_8cpp__incl.map [.f..tp...]`
- `docs/html/association__functions_8cpp__incl.md5 [.f..tp...]`
- `docs/html/association__functions_8cpp__incl.png [.f..tp...]`
- `docs/html/association__functions_8hpp.html [.f..tp...]`
- `docs/html/association__functions_8hpp__dep__incl.map [.f..tp...]`
- `docs/html/association__functions_8hpp__dep__incl.md5 [.f..tp...]`
- `docs/html/association__functions_8hpp__dep__incl.png [.f..tp...]`
- `docs/html/association__functions_8hpp__incl.map [.f..tp...]`
- `docs/html/association__functions_8hpp__incl.md5 [.f..tp...]`
- `docs/html/association__functions_8hpp__incl.png [.f..tp...]`
- `docs/html/association__functions_8hpp_source.html [.f..tp...]`
- `docs/html/association__quaternion__functions_8hpp_source.html [.f..tp...]`
- `docs/html/association__vector__functions_8hpp.html [.f..tp...]`
- `docs/html/association__vector__functions_8hpp_source.html [.f..tp...]`
- `docs/html/bc_s.png [.f..tp...]`
- `docs/html/bdwn.png [.f..tp...]`
- `docs/html/bimolecular__reactions_8cpp.html [.f..tp...]`
- `docs/html/bimolecular__reactions_8hpp.html [.f..tp...]`
- `docs/html/bimolecular__reactions_8hpp_source.html [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp.html [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp__dep__incl.map [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp__dep__incl.md5 [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp__dep__incl.png [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp__incl.map [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp__incl.md5 [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp__incl.png [.f..tp...]`
- `docs/html/bngl__parser__classes_8hpp_source.html [.f..tp...]`
- `docs/html/bngl__parser__member__functions_8cpp.html [.f..tp...]`
- `docs/html/bngl__parser__member__functions_8cpp__incl.map [.f..tp...]`
- `docs/html/bngl__parser__member__functions_8cpp__incl.md5 [.f..tp...]`
- `docs/html/bngl__parser__member__functions_8cpp__incl.png [.f..tp...]`
- `docs/html/changelog.html [.f..tp...]`
- `docs/html/class___coord_8cpp.html [.f..tp...]`
- `docs/html/class___coord_8hpp.html [.f..tp...]`
- `docs/html/class___coord_8hpp_source.html [.f..tp...]`
- `docs/html/class___m_d_timer_8hpp_source.html [.f..tp...]`
- `docs/html/class___mol_template_8hpp.html [.f..tp...]`
- `docs/html/class___mol_template_8hpp_source.html [.f..tp...]`
- `docs/html/class___molecule___complex_8hpp_source.html [.f..tp...]`
- `docs/html/class___observable_8cpp.html [.f..tp...]`
- `docs/html/class___observable_8hpp.html [.f..tp...]`
- `docs/html/class___observable_8hpp_source.html [.f..tp...]`
- `docs/html/class___parameters_8hpp_source.html [.f..tp...]`
- `docs/html/class___quat_8cpp.html [.f..tp...]`
- `docs/html/class___quat_8hpp.html [.f..tp...]`
- `docs/html/class___quat_8hpp_source.html [.f..tp...]`
- `docs/html/class___rxns_8hpp_source.html [.f..tp...]`
- `docs/html/class___simul_volume_8cpp.html [.f..tp...]`
- `docs/html/class___simul_volume_8hpp_source.html [.f..tp...]`
- `docs/html/class___vector_8cpp.html [.f..tp...]`
- `docs/html/class___vector_8hpp.html [.f..tp...]`
- `docs/html/class___vector_8hpp_source.html [.f..tp...]`
- `docs/html/class__bngl__parser_8hpp_source.html [.f..tp...]`
- `docs/html/class__coord_8cpp.html [.f..tp...]`
- `docs/html/class__coord_8cpp__incl.map [.f..tp...]`
- `docs/html/class__coord_8cpp__incl.md5 [.f..tp...]`
- `docs/html/class__coord_8cpp__incl.png [.f..tp...]`
- `docs/html/class__coord_8hpp.html [.f..tp...]`
- `docs/html/class__coord_8hpp__dep__incl.map [.f..tp...]`
- `docs/html/class__coord_8hpp__dep__incl.md5 [.f..tp...]`
- `docs/html/class__coord_8hpp__dep__incl.png [.f..tp...]`
- `docs/html/class__coord_8hpp__incl.map [.f..tp...]`
- `docs/html/class__coord_8hpp__incl.md5 [.f..tp...]`
- `docs/html/class__coord_8hpp__incl.png [.f..tp...]`
- `docs/html/class__coord_8hpp_source.html [.f..tp...]`
- `docs/html/class__mdtimer_8hpp_source.html [.f..tp...]`
- `docs/html/class__mol__containers_8hpp_source.html [.f..tp...]`
- `docs/html/class__moltemplate_8hpp.html [.f..tp...]`
- `docs/html/class__moltemplate_8hpp_source.html [.f..tp...]`
- `docs/html/class__parameter_8cpp.html [.f..tp...]`
- `docs/html/class__parameter_8cpp__incl.map [.f..tp...]`
- `docs/html/class__parameter_8cpp__incl.md5 [.f..tp...]`
- `docs/html/class__parameter_8cpp__incl.png [.f..tp...]`
- `docs/html/class__parameter_8hpp_source.html [.f..tp...]`
- `docs/html/class__quaternion_8hpp_source.html [.f..tp...]`
- `docs/html/class__rxns_8hpp_source.html [.f..tp...]`
- `docs/html/class__simulbox_8hpp.html [.f..tp...]`
- `docs/html/class__simulbox_8hpp_source.html [.f..tp...]`
- `docs/html/class__simulbox__functions_8cpp.html [.f..tp...]`
- `docs/html/class__subs_8cpp.html [.f..tp...]`
- `docs/html/class__subs_8cpp__incl.map [.f..tp...]`
- `docs/html/class__subs_8cpp__incl.md5 [.f..tp...]`
- `docs/html/class__subs_8cpp__incl.png [.f..tp...]`
- `docs/html/class__vector_8hpp.html [.f..tp...]`
- `docs/html/class__vector_8hpp__dep__incl.map [.f..tp...]`
- `docs/html/class__vector_8hpp__dep__incl.md5 [.f..tp...]`
- `docs/html/class__vector_8hpp__dep__incl.png [.f..tp...]`
- `docs/html/class__vector_8hpp__incl.map [.f..tp...]`
- `docs/html/class__vector_8hpp__incl.md5 [.f..tp...]`
- `docs/html/class__vector_8hpp__incl.png [.f..tp...]`
- `docs/html/class__vector_8hpp_source.html [.f..tp...]`
- `docs/html/class__vector__functions_8cpp.html [.f..tp...]`
- `docs/html/class__vector__functions_8cpp__incl.map [.f..tp...]`
- `docs/html/class__vector__functions_8cpp__incl.md5 [.f..tp...]`
- `docs/html/class__vector__functions_8cpp__incl.png [.f..tp...]`
- `docs/html/class_a.html [.f..tp...]`
- `docs/html/class_arg_exception.html [.f..tp...]`
- `docs/html/class_arg_exception.png [.f..tp...]`
- `docs/html/class_arg_exception__coll__graph.map [.f..tp...]`
- `docs/html/class_arg_exception__coll__graph.md5 [.f..tp...]`
- `docs/html/class_arg_exception__coll__graph.png [.f..tp...]`
- `docs/html/class_arg_exception__inherit__graph.map [.f..tp...]`
- `docs/html/class_arg_exception__inherit__graph.md5 [.f..tp...]`
- `docs/html/class_arg_exception__inherit__graph.png [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception.html [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception.png [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception__coll__graph.map [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception__coll__graph.md5 [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception__coll__graph.png [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception__inherit__graph.map [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception__inherit__graph.md5 [.f..tp...]`
- `docs/html/class_arg_o_o_f_exception__inherit__graph.png [.f..tp...]`
- `docs/html/class_coordinate.html [.f..tp...]`
- `docs/html/class_create_destruct_rxn.html [.f..tp...]`
- `docs/html/class_create_destruct_rxn.png [.f..tp...]`
- `docs/html/class_exception.html [.f..tp...]`
- `docs/html/class_exception.png [.f..tp...]`
- `docs/html/class_exception__coll__graph.map [.f..tp...]`
- `docs/html/class_exception__coll__graph.md5 [.f..tp...]`
- `docs/html/class_exception__coll__graph.png [.f..tp...]`
- `docs/html/class_exception__inherit__graph.map [.f..tp...]`
- `docs/html/class_exception__inherit__graph.md5 [.f..tp...]`
- `docs/html/class_exception__inherit__graph.png [.f..tp...]`
- `docs/html/class_parameters.html [.f..tp...]`
- `docs/html/class_parameters__coll__graph.map [.f..tp...]`
- `docs/html/class_parameters__coll__graph.md5 [.f..tp...]`
- `docs/html/class_parameters__coll__graph.png [.f..tp...]`
- `docs/html/class_rxn_base.html [.f..tp...]`
- `docs/html/class_rxn_base.png [.f..tp...]`
- `docs/html/class_rxn_base__coll__graph.map [.f..tp...]`
- `docs/html/class_rxn_base__coll__graph.md5 [.f..tp...]`
- `docs/html/class_rxn_base__coll__graph.png [.f..tp...]`
- `docs/html/class_rxn_base__inherit__graph.map [.f..tp...]`
- `docs/html/class_rxn_base__inherit__graph.md5 [.f..tp...]`
- `docs/html/class_rxn_base__inherit__graph.png [.f..tp...]`
- `docs/html/classdocument.html [.f..tp...]`
- `docs/html/classes.html [.f..tp...]`
- `docs/html/classes_8hpp.html [.f..tp...]`
- `docs/html/classes_8hpp__dep__incl.map [.f..tp...]`
- `docs/html/classes_8hpp__dep__incl.md5 [.f..tp...]`
- `docs/html/classes_8hpp__dep__incl.png [.f..tp...]`
- `docs/html/classes_8hpp__incl.map [.f..tp...]`
- `docs/html/classes_8hpp__incl.md5 [.f..tp...]`
- `docs/html/classes_8hpp__incl.png [.f..tp...]`
- `docs/html/classes_8hpp_source.html [.f..tp...]`
- `docs/html/closed.png [.f..tp...]`
- `docs/html/constants_8hpp_source.html [.f..tp...]`
- `docs/html/creation_8cpp.html [.f..tp...]`
- `docs/html/creation_8hpp.html [.f..tp...]`
- `docs/html/creation_8hpp_source.html [.f..tp...]`
- `docs/html/development.html [.f..tp...]`
- `docs/html/dir_000004_000003.html [.f..tp...]`
- `docs/html/dir_000005_000003.html [.f..tp...]`
- `docs/html/dir_000006_000003.html [.f..tp...]`
- `docs/html/dir_000007_000003.html [.f..tp...]`
- `docs/html/dir_0435e8edf355b71c77084e5baa379b52.html [.f..tp...]`
- `docs/html/dir_157e42751e6915f3d0604c2ba7028392.html [.f..tp...]`
- `docs/html/dir_1aa96dc6d53c936ca0f7b620e65d6df8.html [.f..tp...]`
- `docs/html/dir_238c0345f0059a07202ddaa36b6fbeaf.html [.f..tp...]`
- `docs/html/dir_24a8b6e96e8d7e949def50b6991039f4.html [.f..tp...]`
- `docs/html/dir_31a61863f657b4eebe6798c8f01c8111.html [.f..tp...]`
- `docs/html/dir_3f2e7598d726066a8e307b1ca0a6cf8a.html [.f..tp...]`
- `docs/html/dir_429fdd4f192b9ccb33740819c5d3a63b.html [.f..tp...]`
- `docs/html/dir_49e56c817e5e54854c35e136979f97ca.html [.f..tp...]`
- `docs/html/dir_5148c598f6de7e391fbd67a54f93a66b.html [.f..tp...]`
- `docs/html/dir_5148c598f6de7e391fbd67a54f93a66b_dep.map [.f..tp...]`
- `docs/html/dir_5148c598f6de7e391fbd67a54f93a66b_dep.md5 [.f..tp...]`
- `docs/html/dir_5148c598f6de7e391fbd67a54f93a66b_dep.png [.f..tp...]`
- `docs/html/dir_609ea522d7846d54795379d3181f3c62.html [.f..tp...]`
- `docs/html/dir_64beed06bc3232e38dd1bdac6ce506fe.html [.f..tp...]`
- `docs/html/dir_64beed06bc3232e38dd1bdac6ce506fe_dep.map [.f..tp...]`
- `docs/html/dir_64beed06bc3232e38dd1bdac6ce506fe_dep.md5 [.f..tp...]`
- `docs/html/dir_64beed06bc3232e38dd1bdac6ce506fe_dep.png [.f..tp...]`
- `docs/html/dir_676968bf5a8bfb646ff33491a62f6ade.html [.f..tp...]`
- `docs/html/dir_6ba07254a5c7884c974b0338a69e0bc8.html [.f..tp...]`
- `docs/html/dir_7436badfbf8a7ed1cf5719a5f7b7b06b.html [.f..tp...]`
- `docs/html/dir_785c93a950cc12bc56dfc18c98c3048e.html [.f..tp...]`
- `docs/html/dir_7bdd6bd8b616be3d42595e1ea8da8ea4.html [.f..tp...]`
- `docs/html/dir_84f6f846b7c4db3f9285e5437e173d5d.html [.f..tp...]`
- `docs/html/dir_8835e511f42c1eeb82adeb31d8477a4e.html [.f..tp...]`
- `docs/html/dir_8835e511f42c1eeb82adeb31d8477a4e_dep.map [.f..tp...]`
- `docs/html/dir_8835e511f42c1eeb82adeb31d8477a4e_dep.md5 [.f..tp...]`
- `docs/html/dir_8835e511f42c1eeb82adeb31d8477a4e_dep.png [.f..tp...]`
- `docs/html/dir_9199c10b3bca305f0942bd0c4bc154a9.html [.f..tp...]`
- `docs/html/dir_973050d4f82d7fdb7cb9757be51fb98d.html [.f..tp...]`
- `docs/html/dir_a0fed543900c5d09f0008580da65bf73.html [.f..tp...]`
- `docs/html/dir_a0fed543900c5d09f0008580da65bf73_dep.map [.f..tp...]`
- `docs/html/dir_a0fed543900c5d09f0008580da65bf73_dep.md5 [.f..tp...]`
- `docs/html/dir_a0fed543900c5d09f0008580da65bf73_dep.png [.f..tp...]`
- `docs/html/dir_a545885098ec48be19b11df73fe56da4.html [.f..tp...]`
- `docs/html/dir_a545885098ec48be19b11df73fe56da4_dep.map [.f..tp...]`
- `docs/html/dir_a545885098ec48be19b11df73fe56da4_dep.md5 [.f..tp...]`
- `docs/html/dir_a545885098ec48be19b11df73fe56da4_dep.png [.f..tp...]`
- `docs/html/dir_ac33cf0d87067c2f5b205ac434b39f87.html [.f..tp...]`
- `docs/html/dir_b4019b81bf41b317a7c625f352952ae0.html [.f..tp...]`
- `docs/html/dir_c54eb78537d5dddf4290e9c0463798ed.html [.f..tp...]`
- `docs/html/dir_c54eb78537d5dddf4290e9c0463798ed_dep.map [.f..tp...]`
- `docs/html/dir_c54eb78537d5dddf4290e9c0463798ed_dep.md5 [.f..tp...]`
- `docs/html/dir_c54eb78537d5dddf4290e9c0463798ed_dep.png [.f..tp...]`
- `docs/html/dir_dfc3d69610b75d059cb6f1cd4a8d10ed.html [.f..tp...]`
- `docs/html/dir_ee15a1d2330c8f88c1d64fa183f288c3.html [.f..tp...]`
- `docs/html/dir_ee15a1d2330c8f88c1d64fa183f288c3_dep.map [.f..tp...]`
- `docs/html/dir_ee15a1d2330c8f88c1d64fa183f288c3_dep.md5 [.f..tp...]`
- `docs/html/dir_ee15a1d2330c8f88c1d64fa183f288c3_dep.png [.f..tp...]`
- `docs/html/dissociation_8cpp.html [.f..tp...]`
- `docs/html/dissociation_8hpp.html [.f..tp...]`
- `docs/html/dissociation_8hpp_source.html [.f..tp...]`
- `docs/html/doc.png [.f..tp...]`
- `docs/html/doxygen.css [.f..tp...]`
- `docs/html/doxygen.html [.f..tp...]`
- `docs/html/doxygen.png [.f..tp...]`
- `docs/html/doxygen_html_style.css [.f..tp...]`
- `docs/html/dynsections.js [.f..tp...]`
- `docs/html/error__handling_8hpp.html [.f..tp...]`
- `docs/html/error__handling_8hpp__dep__incl.map [.f..tp...]`
- `docs/html/error__handling_8hpp__dep__incl.md5 [.f..tp...]`
- `docs/html/error__handling_8hpp__dep__incl.png [.f..tp...]`
- `docs/html/error__handling_8hpp__incl.map [.f..tp...]`
- `docs/html/error__handling_8hpp__incl.md5 [.f..tp...]`
- `docs/html/error__handling_8hpp__incl.png [.f..tp...]`
- `docs/html/error__handling_8hpp_source.html [.f..tp...]`
- `docs/html/evaluate__binding__pair_8cpp.html [.f..tp...]`
- `docs/html/examples_page.html [.f..tp...]`
- `docs/html/files.html [.f..tp...]`
- `docs/html/folderclosed.png [.f..tp...]`
- `docs/html/folderopen.png [.f..tp...]`
- `docs/html/form_0.png [.f..tp...]`
- `docs/html/form_1.png [.f..tp...]`
- `docs/html/form_10.png [.f..tp...]`
- `docs/html/form_11.png [.f..tp...]`
- `docs/html/form_12.png [.f..tp...]`
- `docs/html/form_13.png [.f..tp...]`
- `docs/html/form_14.png [.f..tp...]`
- `docs/html/form_15.png [.f..tp...]`
- `docs/html/form_16.png [.f..tp...]`
- `docs/html/form_17.png [.f..tp...]`
- `docs/html/form_18.png [.f..tp...]`
- `docs/html/form_19.png [.f..tp...]`
- `docs/html/form_2.png [.f..tp...]`
- `docs/html/form_20.png [.f..tp...]`
- `docs/html/form_21.png [.f..tp...]`
- `docs/html/form_22.png [.f..tp...]`
- `docs/html/form_23.png [.f..tp...]`
- `docs/html/form_24.png [.f..tp...]`
- `docs/html/form_25.png [.f..tp...]`
- `docs/html/form_26.png [.f..tp...]`
- `docs/html/form_27.png [.f..tp...]`
- `docs/html/form_28.png [.f..tp...]`
- `docs/html/form_29.png [.f..tp...]`
- `docs/html/form_3.png [.f..tp...]`
- `docs/html/form_30.png [.f..tp...]`
- `docs/html/form_31.png [.f..tp...]`
- `docs/html/form_32.png [.f..tp...]`
- `docs/html/form_4.png [.f..tp...]`
- `docs/html/form_5.png [.f..tp...]`
- `docs/html/form_6.png [.f..tp...]`
- `docs/html/form_7.png [.f..tp...]`
- `docs/html/form_8.png [.f..tp...]`
- `docs/html/form_9.png [.f..tp...]`
- `docs/html/formula.repository [.f..tp...]`
- `docs/html/functions.html [.f..tp...]`
- `docs/html/functions_b.html [.f..tp...]`
- `docs/html/functions_c.html [.f..tp...]`
- `docs/html/functions_d.html [.f..tp...]`
- `docs/html/functions_e.html [.f..tp...]`
- `docs/html/functions_f.html [.f..tp...]`
- `docs/html/functions_func.html [.f..tp...]`
- `docs/html/functions_h.html [.f..tp...]`
- `docs/html/functions_i.html [.f..tp...]`
- `docs/html/functions_k.html [.f..tp...]`
- `docs/html/functions_m.html [.f..tp...]`
- `docs/html/functions_n.html [.f..tp...]`
- `docs/html/functions_o.html [.f..tp...]`
- `docs/html/functions_p.html [.f..tp...]`
- `docs/html/functions_r.html [.f..tp...]`
- `docs/html/functions_s.html [.f..tp...]`
- `docs/html/functions_t.html [.f..tp...]`
- `docs/html/functions_u.html [.f..tp...]`
- `docs/html/functions_vars.html [.f..tp...]`
- `docs/html/functions_w.html [.f..tp...]`
- `docs/html/functions_x.html [.f..tp...]`
- `docs/html/functions_y.html [.f..tp...]`
- `docs/html/functions_z.html [.f..tp...]`
- `docs/html/get__distance_8cpp.html [.f..tp...]`
- `docs/html/globals.html [.f..tp...]`
- `docs/html/globals_8hpp_source.html [.f..tp...]`
- `docs/html/globals_enum.html [.f..tp...]`
- `docs/html/globals_func.html [.f..tp...]`
- `docs/html/graph_legend.html [.f..tp...]`
- `docs/html/graph_legend.md5 [.f..tp...]`
- `docs/html/graph_legend.png [.f..tp...]`
- `docs/html/group___associate.html [.f..tp...]`
- `docs/html/group___associate.map [.f..tp...]`
- `docs/html/group___associate.md5 [.f..tp...]`
- `docs/html/group___associate.png [.f..tp...]`
- `docs/html/group___d_reactions.html [.f..tp...]`
- `docs/html/group___i_o.html [.f..tp...]`
- `docs/html/group___parser.html [.f..tp...]`
- `docs/html/group___reactions.html [.f..tp...]`
- `docs/html/group___reactions.map [.f..tp...]`
- `docs/html/group___reactions.md5 [.f..tp...]`
- `docs/html/group___reactions.png [.f..tp...]`
- `docs/html/group___simul_classes.html [.f..tp...]`
- `docs/html/group___simul_classes.map [.f..tp...]`
- `docs/html/group___simul_classes.md5 [.f..tp...]`
- `docs/html/group___simul_classes.png [.f..tp...]`
- `docs/html/group___species_tracker.html [.f..tp...]`
- `docs/html/group___system_setup.html [.f..tp...]`
- `docs/html/group___templates.html [.f..tp...]`
- `docs/html/group___text.html [.f..tp...]`
- `docs/html/hierarchy.html [.f..tp...]`
- `docs/html/index.html [.f..tp...]`
- `docs/html/inherit_graph_0.map [.f..tp...]`
- `docs/html/inherit_graph_0.md5 [.f..tp...]`
- `docs/html/inherit_graph_0.png [.f..tp...]`
- `docs/html/inherit_graph_1.map [.f..tp...]`
- `docs/html/inherit_graph_1.md5 [.f..tp...]`
- `docs/html/inherit_graph_1.png [.f..tp...]`
- `docs/html/inherit_graph_10.map [.f..tp...]`
- `docs/html/inherit_graph_10.md5 [.f..tp...]`
- `docs/html/inherit_graph_10.png [.f..tp...]`
- `docs/html/inherit_graph_11.map [.f..tp...]`
- `docs/html/inherit_graph_11.md5 [.f..tp...]`
- `docs/html/inherit_graph_11.png [.f..tp...]`
- `docs/html/inherit_graph_12.map [.f..tp...]`
- `docs/html/inherit_graph_12.md5 [.f..tp...]`
- `docs/html/inherit_graph_12.png [.f..tp...]`
- `docs/html/inherit_graph_13.map [.f..tp...]`
- `docs/html/inherit_graph_13.md5 [.f..tp...]`
- `docs/html/inherit_graph_13.png [.f..tp...]`
- `docs/html/inherit_graph_14.map [.f..tp...]`
- `docs/html/inherit_graph_14.md5 [.f..tp...]`
- `docs/html/inherit_graph_14.png [.f..tp...]`
- `docs/html/inherit_graph_15.map [.f..tp...]`
- `docs/html/inherit_graph_15.md5 [.f..tp...]`
- `docs/html/inherit_graph_15.png [.f..tp...]`
- `docs/html/inherit_graph_16.map [.f..tp...]`
- `docs/html/inherit_graph_16.md5 [.f..tp...]`
- `docs/html/inherit_graph_16.png [.f..tp...]`
- `docs/html/inherit_graph_17.map [.f..tp...]`
- `docs/html/inherit_graph_17.md5 [.f..tp...]`
- `docs/html/inherit_graph_17.png [.f..tp...]`
- `docs/html/inherit_graph_18.map [.f..tp...]`
- `docs/html/inherit_graph_18.md5 [.f..tp...]`
- `docs/html/inherit_graph_18.png [.f..tp...]`
- `docs/html/inherit_graph_19.map [.f..tp...]`
- `docs/html/inherit_graph_19.md5 [.f..tp...]`
- `docs/html/inherit_graph_19.png [.f..tp...]`
- `docs/html/inherit_graph_2.map [.f..tp...]`
- `docs/html/inherit_graph_2.md5 [.f..tp...]`
- `docs/html/inherit_graph_2.png [.f..tp...]`
- `docs/html/inherit_graph_20.map [.f..tp...]`
- `docs/html/inherit_graph_20.md5 [.f..tp...]`
- `docs/html/inherit_graph_20.png [.f..tp...]`
- `docs/html/inherit_graph_21.map [.f..tp...]`
- `docs/html/inherit_graph_21.md5 [.f..tp...]`
- `docs/html/inherit_graph_21.png [.f..tp...]`
- `docs/html/inherit_graph_22.map [.f..tp...]`
- `docs/html/inherit_graph_22.md5 [.f..tp...]`
- `docs/html/inherit_graph_22.png [.f..tp...]`
- `docs/html/inherit_graph_23.map [.f..tp...]`
- `docs/html/inherit_graph_23.md5 [.f..tp...]`
- `docs/html/inherit_graph_23.png [.f..tp...]`
- `docs/html/inherit_graph_24.map [.f..tp...]`
- `docs/html/inherit_graph_24.md5 [.f..tp...]`
- `docs/html/inherit_graph_24.png [.f..tp...]`
- `docs/html/inherit_graph_3.map [.f..tp...]`
- `docs/html/inherit_graph_3.md5 [.f..tp...]`
- `docs/html/inherit_graph_3.png [.f..tp...]`
- `docs/html/inherit_graph_4.map [.f..tp...]`
- `docs/html/inherit_graph_4.md5 [.f..tp...]`
- `docs/html/inherit_graph_4.png [.f..tp...]`
- `docs/html/inherit_graph_5.map [.f..tp...]`
- `docs/html/inherit_graph_5.md5 [.f..tp...]`
- `docs/html/inherit_graph_5.png [.f..tp...]`
- `docs/html/inherit_graph_6.map [.f..tp...]`
- `docs/html/inherit_graph_6.md5 [.f..tp...]`
- `docs/html/inherit_graph_6.png [.f..tp...]`
- `docs/html/inherit_graph_7.map [.f..tp...]`
- `docs/html/inherit_graph_7.md5 [.f..tp...]`
- `docs/html/inherit_graph_7.png [.f..tp...]`
- `docs/html/inherit_graph_8.map [.f..tp...]`
- `docs/html/inherit_graph_8.md5 [.f..tp...]`
- `docs/html/inherit_graph_8.png [.f..tp...]`
- `docs/html/inherit_graph_9.map [.f..tp...]`
- `docs/html/inherit_graph_9.md5 [.f..tp...]`
- `docs/html/inherit_graph_9.png [.f..tp...]`
- `docs/html/inherits.html [.f..tp...]`
- `docs/html/input.html [.f..tp...]`
- `docs/html/io_8cpp.html [.f..tp...]`
- `docs/html/io_8hpp.html [.f..tp...]`
- `docs/html/io_8hpp_source.html [.f..tp...]`
- `docs/html/jquery.js [.f..tp...]`
- `docs/html/math__functions_8hpp.html [.f..tp...]`
- `docs/html/math__functions_8hpp_source.html [.f..tp...]`
- `docs/html/matrix_8hpp.html [.f..tp...]`
- `docs/html/matrix_8hpp_source.html [.f..tp...]`
- `docs/html/matrix__functions_8cpp.html [.f..tp...]`
- `docs/html/menu.js [.f..tp...]`
- `docs/html/menudata.js [.f..tp...]`
- `docs/html/modules.html [.f..tp...]`
- `docs/html/namespacedocument.html [.f..tp...]`
- `docs/html/namespaces.html [.f..tp...]`
- `docs/html/nav_f.png [.f..tp...]`
- `docs/html/nav_g.png [.f..tp...]`
- `docs/html/nav_h.png [.f..tp...]`
- `docs/html/omega.pdf [.f..tp...]`
- `docs/html/open.png [.f..tp...]`
- `docs/html/pages.html [.f..tp...]`
- `docs/html/parse__input_8cpp.html [.f..tp...]`
- `docs/html/parser__functions_8hpp_source.html [.f..tp...]`
- `docs/html/parser__math__functions_8hpp.html [.f..tp...]`
- `docs/html/parser__math__functions_8hpp_source.html [.f..tp...]`
- `docs/html/phi.pdf [.f..tp...]`
- `docs/html/pucker_clat_example.pdf [.f..tp...]`
- `docs/html/rand__gsl_8h_source.html [.f..tp...]`
- `docs/html/rand__gsl_8hpp_source.html [.f..tp...]`
- `docs/html/rotation.pdf [.f..tp...]`
- `docs/html/shared__reaction__functions_8cpp.html [.f..tp...]`
- `docs/html/shared__reaction__functions_8hpp.html [.f..tp...]`
- `docs/html/shared__reaction__functions_8hpp_source.html [.f..tp...]`
- `docs/html/species.png [.f..tp...]`
- `docs/html/splitbar.png [.f..tp...]`
- `docs/html/struct_angles.html [.f..tp...]`
- `docs/html/struct_back_rxn.html [.f..tp...]`
- `docs/html/struct_back_rxn.png [.f..tp...]`
- `docs/html/struct_back_rxn__coll__graph.map [.f..tp...]`
- `docs/html/struct_back_rxn__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_back_rxn__coll__graph.png [.f..tp...]`
- `docs/html/struct_back_rxn__inherit__graph.map [.f..tp...]`
- `docs/html/struct_back_rxn__inherit__graph.md5 [.f..tp...]`
- `docs/html/struct_back_rxn__inherit__graph.png [.f..tp...]`
- `docs/html/struct_bi_mol_data.html [.f..tp...]`
- `docs/html/struct_complex.html [.f..tp...]`
- `docs/html/struct_complex__coll__graph.map [.f..tp...]`
- `docs/html/struct_complex__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_complex__coll__graph.png [.f..tp...]`
- `docs/html/struct_coord.html [.f..tp...]`
- `docs/html/struct_coord.png [.f..tp...]`
- `docs/html/struct_coord__inherit__graph.map [.f..tp...]`
- `docs/html/struct_coord__inherit__graph.md5 [.f..tp...]`
- `docs/html/struct_coord__inherit__graph.png [.f..tp...]`
- `docs/html/struct_create_destruct_rxn.html [.f..tp...]`
- `docs/html/struct_create_destruct_rxn.png [.f..tp...]`
- `docs/html/struct_create_destruct_rxn_1_1_create_destruct_mol.html [.f..tp...]`
- `docs/html/struct_create_destruct_rxn__coll__graph.map [.f..tp...]`
- `docs/html/struct_create_destruct_rxn__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_create_destruct_rxn__coll__graph.png [.f..tp...]`
- `docs/html/struct_create_destruct_rxn__inherit__graph.map [.f..tp...]`
- `docs/html/struct_create_destruct_rxn__inherit__graph.md5 [.f..tp...]`
- `docs/html/struct_create_destruct_rxn__inherit__graph.png [.f..tp...]`
- `docs/html/struct_create_destruct_rxn_new.html [.f..tp...]`
- `docs/html/struct_forward_rxn.html [.f..tp...]`
- `docs/html/struct_forward_rxn.png [.f..tp...]`
- `docs/html/struct_forward_rxn_1_1_angles.html [.f..tp...]`
- `docs/html/struct_forward_rxn__coll__graph.map [.f..tp...]`
- `docs/html/struct_forward_rxn__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_forward_rxn__coll__graph.png [.f..tp...]`
- `docs/html/struct_forward_rxn__inherit__graph.map [.f..tp...]`
- `docs/html/struct_forward_rxn__inherit__graph.md5 [.f..tp...]`
- `docs/html/struct_forward_rxn__inherit__graph.png [.f..tp...]`
- `docs/html/struct_iface_info.html [.f..tp...]`
- `docs/html/struct_int_coord_cont.html [.f..tp...]`
- `docs/html/struct_int_coord_cont_1_1_mol.html [.f..tp...]`
- `docs/html/struct_int_coord_cont_1_1_mol__coll__graph.map [.f..tp...]`
- `docs/html/struct_int_coord_cont_1_1_mol__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_int_coord_cont_1_1_mol__coll__graph.png [.f..tp...]`
- `docs/html/struct_int_coord_cont__coll__graph.map [.f..tp...]`
- `docs/html/struct_int_coord_cont__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_int_coord_cont__coll__graph.png [.f..tp...]`
- `docs/html/struct_integrand_params.html [.f..tp...]`
- `docs/html/struct_interaction.html [.f..tp...]`
- `docs/html/struct_interface.html [.f..tp...]`
- `docs/html/struct_interface_1_1_index.html [.f..tp...]`
- `docs/html/struct_interface_1_1_rxn_state.html [.f..tp...]`
- `docs/html/struct_interface_1_1_state.html [.f..tp...]`
- `docs/html/struct_interface_1_1_state__coll__graph.map [.f..tp...]`
- `docs/html/struct_interface_1_1_state__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_interface_1_1_state__coll__graph.png [.f..tp...]`
- `docs/html/struct_interface__coll__graph.map [.f..tp...]`
- `docs/html/struct_interface__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_interface__coll__graph.png [.f..tp...]`
- `docs/html/struct_m_d_timer.html [.f..tp...]`
- `docs/html/struct_mol_template.html [.f..tp...]`
- `docs/html/struct_mol_template__coll__graph.map [.f..tp...]`
- `docs/html/struct_mol_template__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_mol_template__coll__graph.png [.f..tp...]`
- `docs/html/struct_molecule.html [.f..tp...]`
- `docs/html/struct_molecule_1_1_encounter.html [.f..tp...]`
- `docs/html/struct_molecule_1_1_iface.html [.f..tp...]`
- `docs/html/struct_molecule_1_1_iface__coll__graph.map [.f..tp...]`
- `docs/html/struct_molecule_1_1_iface__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_molecule_1_1_iface__coll__graph.png [.f..tp...]`
- `docs/html/struct_molecule_1_1_interaction.html [.f..tp...]`
- `docs/html/struct_molecule__coll__graph.map [.f..tp...]`
- `docs/html/struct_molecule__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_molecule__coll__graph.png [.f..tp...]`
- `docs/html/struct_observable.html [.f..tp...]`
- `docs/html/struct_observable_1_1_constituents.html [.f..tp...]`
- `docs/html/struct_observable_1_1_constituents_1_1_iface.html [.f..tp...]`
- `docs/html/struct_parameters.html [.f..tp...]`
- `docs/html/struct_parameters_1_1_debug.html [.f..tp...]`
- `docs/html/struct_parameters_1_1_water_box.html [.f..tp...]`
- `docs/html/struct_parsed_mol.html [.f..tp...]`
- `docs/html/struct_parsed_mol_1_1_iface_info.html [.f..tp...]`
- `docs/html/struct_parsed_mol_1_1_iface_info__coll__graph.map [.f..tp...]`
- `docs/html/struct_parsed_mol_1_1_iface_info__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_parsed_mol_1_1_iface_info__coll__graph.png [.f..tp...]`
- `docs/html/struct_parsed_mol__coll__graph.map [.f..tp...]`
- `docs/html/struct_parsed_mol__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_parsed_mol__coll__graph.png [.f..tp...]`
- `docs/html/struct_parsed_rxn.html [.f..tp...]`
- `docs/html/struct_parsed_rxn.png [.f..tp...]`
- `docs/html/struct_parsed_rxn__coll__graph.map [.f..tp...]`
- `docs/html/struct_parsed_rxn__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_parsed_rxn__coll__graph.png [.f..tp...]`
- `docs/html/struct_parsed_rxn__inherit__graph.map [.f..tp...]`
- `docs/html/struct_parsed_rxn__inherit__graph.md5 [.f..tp...]`
- `docs/html/struct_parsed_rxn__inherit__graph.png [.f..tp...]`
- `docs/html/struct_quat.html [.f..tp...]`
- `docs/html/struct_rxn_base.html [.f..tp...]`
- `docs/html/struct_rxn_base.png [.f..tp...]`
- `docs/html/struct_rxn_base_1_1_coupled_rxn.html [.f..tp...]`
- `docs/html/struct_rxn_base_1_1_rate_state.html [.f..tp...]`
- `docs/html/struct_rxn_base_1_1_rate_state__coll__graph.map [.f..tp...]`
- `docs/html/struct_rxn_base_1_1_rate_state__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_rxn_base_1_1_rate_state__coll__graph.png [.f..tp...]`
- `docs/html/struct_rxn_iface.html [.f..tp...]`
- `docs/html/struct_rxn_iface__coll__graph.map [.f..tp...]`
- `docs/html/struct_rxn_iface__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_rxn_iface__coll__graph.png [.f..tp...]`
- `docs/html/struct_rxn_state.html [.f..tp...]`
- `docs/html/struct_simul_box.html [.f..tp...]`
- `docs/html/struct_simul_box_1_1_dimensions.html [.f..tp...]`
- `docs/html/struct_simul_box_1_1_sub_box.html [.f..tp...]`
- `docs/html/struct_simul_volume.html [.f..tp...]`
- `docs/html/struct_simul_volume_1_1_dimensions.html [.f..tp...]`
- `docs/html/struct_simul_volume_1_1_sub_volume.html [.f..tp...]`
- `docs/html/struct_species_tracker_1_1_observable.html [.f..tp...]`
- `docs/html/struct_species_tracker_1_1_observable_1_1_constituent.html [.f..tp...]`
- `docs/html/struct_species_tracker_1_1_observable_1_1_constituent_1_1_iface.html [.f..tp...]`
- `docs/html/struct_vector.html [.f..tp...]`
- `docs/html/struct_vector.png [.f..tp...]`
- `docs/html/struct_vector__coll__graph.map [.f..tp...]`
- `docs/html/struct_vector__coll__graph.md5 [.f..tp...]`
- `docs/html/struct_vector__coll__graph.png [.f..tp...]`
- `docs/html/struct_vector__inherit__graph.map [.f..tp...]`
- `docs/html/struct_vector__inherit__graph.md5 [.f..tp...]`
- `docs/html/struct_vector__inherit__graph.png [.f..tp...]`
- `docs/html/structor.html [.f..tp...]`
- `docs/html/sync_off.png [.f..tp...]`
- `docs/html/sync_on.png [.f..tp...]`
- `docs/html/system__setup_8hpp.html [.f..tp...]`
- `docs/html/system__setup_8hpp__incl.map [.f..tp...]`
- `docs/html/system__setup_8hpp__incl.md5 [.f..tp...]`
- `docs/html/system__setup_8hpp__incl.png [.f..tp...]`
- `docs/html/system__setup_8hpp_source.html [.f..tp...]`
- `docs/html/tab_a.png [.f..tp...]`
- `docs/html/tab_b.png [.f..tp...]`
- `docs/html/tab_h.png [.f..tp...]`
- `docs/html/tab_s.png [.f..tp...]`
- `docs/html/tabs.css [.f..tp...]`
- `docs/html/theta.pdf [.f..tp...]`
- `docs/html/trajectory__functions_8cpp.html [.f..tp...]`
- `docs/html/trajectory__functions_8hpp.html [.f..tp...]`
- `docs/html/trajectory__functions_8hpp_source.html [.f..tp...]`
- `docs/html/unimolecular__reactions_8hpp_source.html [.f..tp...]`
- `docs/html/zeroth_order_creation_avg.png [.f..tp...]`
- `docs/html/zeroth_order_creation_indiv.png [.f..tp...]`
- `docs/html/search/all_0.html [.f..tp...]`
- `docs/html/search/all_0.js [.f..tp...]`
- `docs/html/search/all_1.html [.f..tp...]`
- `docs/html/search/all_1.js [.f..tp...]`
- `docs/html/search/all_10.html [.f..tp...]`
- `docs/html/search/all_10.js [.f..tp...]`
- `docs/html/search/all_11.html [.f..tp...]`
- `docs/html/search/all_11.js [.f..tp...]`
- `docs/html/search/all_12.html [.f..tp...]`
- `docs/html/search/all_12.js [.f..tp...]`
- `docs/html/search/all_13.html [.f..tp...]`
- `docs/html/search/all_13.js [.f..tp...]`
- `docs/html/search/all_14.html [.f..tp...]`
- `docs/html/search/all_14.js [.f..tp...]`
- `docs/html/search/all_15.html [.f..tp...]`
- `docs/html/search/all_15.js [.f..tp...]`
- `docs/html/search/all_16.html [.f..tp...]`
- `docs/html/search/all_16.js [.f..tp...]`
- `docs/html/search/all_17.html [.f..tp...]`
- `docs/html/search/all_17.js [.f..tp...]`
- `docs/html/search/all_18.html [.f..tp...]`
- `docs/html/search/all_18.js [.f..tp...]`
- `docs/html/search/all_19.html [.f..tp...]`
- `docs/html/search/all_19.js [.f..tp...]`
- `docs/html/search/all_2.html [.f..tp...]`
- `docs/html/search/all_2.js [.f..tp...]`
- `docs/html/search/all_3.html [.f..tp...]`
- `docs/html/search/all_3.js [.f..tp...]`
- `docs/html/search/all_4.html [.f..tp...]`
- `docs/html/search/all_4.js [.f..tp...]`
- `docs/html/search/all_5.html [.f..tp...]`
- `docs/html/search/all_5.js [.f..tp...]`
- `docs/html/search/all_6.html [.f..tp...]`
- `docs/html/search/all_6.js [.f..tp...]`
- `docs/html/search/all_7.html [.f..tp...]`
- `docs/html/search/all_7.js [.f..tp...]`
- `docs/html/search/all_8.html [.f..tp...]`
- `docs/html/search/all_8.js [.f..tp...]`
- `docs/html/search/all_9.html [.f..tp...]`
- `docs/html/search/all_9.js [.f..tp...]`
- `docs/html/search/all_a.html [.f..tp...]`
- `docs/html/search/all_a.js [.f..tp...]`
- `docs/html/search/all_b.html [.f..tp...]`
- `docs/html/search/all_b.js [.f..tp...]`
- `docs/html/search/all_c.html [.f..tp...]`
- `docs/html/search/all_c.js [.f..tp...]`
- `docs/html/search/all_d.html [.f..tp...]`
- `docs/html/search/all_d.js [.f..tp...]`
- `docs/html/search/all_e.html [.f..tp...]`
- `docs/html/search/all_e.js [.f..tp...]`
- `docs/html/search/all_f.html [.f..tp...]`
- `docs/html/search/all_f.js [.f..tp...]`
- `docs/html/search/classes_0.html [.f..tp...]`
- `docs/html/search/classes_0.js [.f..tp...]`
- `docs/html/search/classes_1.html [.f..tp...]`
- `docs/html/search/classes_1.js [.f..tp...]`
- `docs/html/search/classes_2.html [.f..tp...]`
- `docs/html/search/classes_2.js [.f..tp...]`
- `docs/html/search/classes_3.html [.f..tp...]`
- `docs/html/search/classes_3.js [.f..tp...]`
- `docs/html/search/classes_4.html [.f..tp...]`
- `docs/html/search/classes_4.js [.f..tp...]`
- `docs/html/search/classes_5.html [.f..tp...]`
- `docs/html/search/classes_5.js [.f..tp...]`
- `docs/html/search/classes_6.html [.f..tp...]`
- `docs/html/search/classes_6.js [.f..tp...]`
- `docs/html/search/classes_7.html [.f..tp...]`
- `docs/html/search/classes_7.js [.f..tp...]`
- `docs/html/search/classes_8.html [.f..tp...]`
- `docs/html/search/classes_8.js [.f..tp...]`
- `docs/html/search/classes_9.html [.f..tp...]`
- `docs/html/search/classes_9.js [.f..tp...]`
- `docs/html/search/classes_a.html [.f..tp...]`
- `docs/html/search/classes_a.js [.f..tp...]`
- `docs/html/search/classes_b.html [.f..tp...]`
- `docs/html/search/classes_b.js [.f..tp...]`
- `docs/html/search/classes_c.html [.f..tp...]`
- `docs/html/search/classes_c.js [.f..tp...]`
- `docs/html/search/classes_d.html [.f..tp...]`
- `docs/html/search/classes_d.js [.f..tp...]`
- `docs/html/search/classes_e.html [.f..tp...]`
- `docs/html/search/classes_e.js [.f..tp...]`
- `docs/html/search/close.png [.f..tp...]`
- `docs/html/search/enums_0.html [.f..tp...]`
- `docs/html/search/enums_0.js [.f..tp...]`
- `docs/html/search/enums_1.html [.f..tp...]`
- `docs/html/search/enums_1.js [.f..tp...]`
- `docs/html/search/enums_2.html [.f..tp...]`
- `docs/html/search/enums_2.js [.f..tp...]`
- `docs/html/search/enums_3.html [.f..tp...]`
- `docs/html/search/enums_3.js [.f..tp...]`
- `docs/html/search/enums_4.html [.f..tp...]`
- `docs/html/search/enums_4.js [.f..tp...]`
- `docs/html/search/enumvalues_0.html [.f..tp...]`
- `docs/html/search/enumvalues_0.js [.f..tp...]`
- `docs/html/search/enumvalues_1.html [.f..tp...]`
- `docs/html/search/enumvalues_1.js [.f..tp...]`
- `docs/html/search/enumvalues_2.html [.f..tp...]`
- `docs/html/search/enumvalues_2.js [.f..tp...]`
- `docs/html/search/enumvalues_3.html [.f..tp...]`
- `docs/html/search/enumvalues_3.js [.f..tp...]`
- `docs/html/search/enumvalues_4.html [.f..tp...]`
- `docs/html/search/enumvalues_4.js [.f..tp...]`
- `docs/html/search/enumvalues_5.html [.f..tp...]`
- `docs/html/search/enumvalues_5.js [.f..tp...]`
- `docs/html/search/enumvalues_6.html [.f..tp...]`
- `docs/html/search/enumvalues_6.js [.f..tp...]`
- `docs/html/search/enumvalues_7.html [.f..tp...]`
- `docs/html/search/enumvalues_7.js [.f..tp...]`
- `docs/html/search/enumvalues_8.html [.f..tp...]`
- `docs/html/search/enumvalues_8.js [.f..tp...]`
- `docs/html/search/files_0.html [.f..tp...]`
- `docs/html/search/files_0.js [.f..tp...]`
- `docs/html/search/files_1.html [.f..tp...]`
- `docs/html/search/files_1.js [.f..tp...]`
- `docs/html/search/files_2.html [.f..tp...]`
- `docs/html/search/files_2.js [.f..tp...]`
- `docs/html/search/files_3.html [.f..tp...]`
- `docs/html/search/files_3.js [.f..tp...]`
- `docs/html/search/files_4.html [.f..tp...]`
- `docs/html/search/files_4.js [.f..tp...]`
- `docs/html/search/files_5.html [.f..tp...]`
- `docs/html/search/files_5.js [.f..tp...]`
- `docs/html/search/files_6.html [.f..tp...]`
- `docs/html/search/files_6.js [.f..tp...]`
- `docs/html/search/files_7.html [.f..tp...]`
- `docs/html/search/files_7.js [.f..tp...]`
- `docs/html/search/files_8.html [.f..tp...]`
- `docs/html/search/files_8.js [.f..tp...]`
- `docs/html/search/files_9.html [.f..tp...]`
- `docs/html/search/files_9.js [.f..tp...]`
- `docs/html/search/functions_0.html [.f..tp...]`
- `docs/html/search/functions_0.js [.f..tp...]`
- `docs/html/search/functions_1.html [.f..tp...]`
- `docs/html/search/functions_1.js [.f..tp...]`
- `docs/html/search/functions_10.html [.f..tp...]`
- `docs/html/search/functions_10.js [.f..tp...]`
- `docs/html/search/functions_11.html [.f..tp...]`
- `docs/html/search/functions_11.js [.f..tp...]`
- `docs/html/search/functions_12.html [.f..tp...]`
- `docs/html/search/functions_12.js [.f..tp...]`
- `docs/html/search/functions_13.html [.f..tp...]`
- `docs/html/search/functions_13.js [.f..tp...]`
- `docs/html/search/functions_14.html [.f..tp...]`
- `docs/html/search/functions_14.js [.f..tp...]`
- `docs/html/search/functions_2.html [.f..tp...]`
- `docs/html/search/functions_2.js [.f..tp...]`
- `docs/html/search/functions_3.html [.f..tp...]`
- `docs/html/search/functions_3.js [.f..tp...]`
- `docs/html/search/functions_4.html [.f..tp...]`
- `docs/html/search/functions_4.js [.f..tp...]`
- `docs/html/search/functions_5.html [.f..tp...]`
- `docs/html/search/functions_5.js [.f..tp...]`
- `docs/html/search/functions_6.html [.f..tp...]`
- `docs/html/search/functions_6.js [.f..tp...]`
- `docs/html/search/functions_7.html [.f..tp...]`
- `docs/html/search/functions_7.js [.f..tp...]`
- `docs/html/search/functions_8.html [.f..tp...]`
- `docs/html/search/functions_8.js [.f..tp...]`
- `docs/html/search/functions_9.html [.f..tp...]`
- `docs/html/search/functions_9.js [.f..tp...]`
- `docs/html/search/functions_a.html [.f..tp...]`
- `docs/html/search/functions_a.js [.f..tp...]`
- `docs/html/search/functions_b.html [.f..tp...]`
- `docs/html/search/functions_b.js [.f..tp...]`
- `docs/html/search/functions_c.html [.f..tp...]`
- `docs/html/search/functions_c.js [.f..tp...]`
- `docs/html/search/functions_d.html [.f..tp...]`
- `docs/html/search/functions_d.js [.f..tp...]`
- `docs/html/search/functions_e.html [.f..tp...]`
- `docs/html/search/functions_e.js [.f..tp...]`
- `docs/html/search/functions_f.html [.f..tp...]`
- `docs/html/search/functions_f.js [.f..tp...]`
- `docs/html/search/groups_0.html [.f..tp...]`
- `docs/html/search/groups_0.js [.f..tp...]`
- `docs/html/search/groups_1.html [.f..tp...]`
- `docs/html/search/groups_1.js [.f..tp...]`
- `docs/html/search/groups_2.html [.f..tp...]`
- `docs/html/search/groups_2.js [.f..tp...]`
- `docs/html/search/groups_3.html [.f..tp...]`
- `docs/html/search/groups_3.js [.f..tp...]`
- `docs/html/search/groups_4.html [.f..tp...]`
- `docs/html/search/groups_4.js [.f..tp...]`
- `docs/html/search/groups_5.html [.f..tp...]`
- `docs/html/search/groups_5.js [.f..tp...]`
- `docs/html/search/groups_6.html [.f..tp...]`
- `docs/html/search/groups_6.js [.f..tp...]`
- `docs/html/search/mag_sel.png [.f..tp...]`
- `docs/html/search/namespaces_0.html [.f..tp...]`
- `docs/html/search/namespaces_0.js [.f..tp...]`
- `docs/html/search/nomatches.html [.f..tp...]`
- `docs/html/search/pages_0.html [.f..tp...]`
- `docs/html/search/pages_0.js [.f..tp...]`
- `docs/html/search/pages_1.html [.f..tp...]`
- `docs/html/search/pages_1.js [.f..tp...]`
- `docs/html/search/pages_2.html [.f..tp...]`
- `docs/html/search/pages_2.js [.f..tp...]`
- `docs/html/search/pages_3.html [.f..tp...]`
- `docs/html/search/pages_3.js [.f..tp...]`
- `docs/html/search/pages_4.html [.f..tp...]`
- `docs/html/search/pages_4.js [.f..tp...]`
- `docs/html/search/pages_5.html [.f..tp...]`
- `docs/html/search/pages_5.js [.f..tp...]`
- `docs/html/search/search.css [.f..tp...]`
- `docs/html/search/search.js [.f..tp...]`
- `docs/html/search/search_l.png [.f..tp...]`
- `docs/html/search/search_m.png [.f..tp...]`
- `docs/html/search/search_r.png [.f..tp...]`
- `docs/html/search/searchdata.js [.f..tp...]`
- `docs/html/search/variables_0.html [.f..tp...]`
- `docs/html/search/variables_0.js [.f..tp...]`
- `docs/html/search/variables_1.html [.f..tp...]`
- `docs/html/search/variables_1.js [.f..tp...]`
- `docs/html/search/variables_10.html [.f..tp...]`
- `docs/html/search/variables_10.js [.f..tp...]`
- `docs/html/search/variables_11.html [.f..tp...]`
- `docs/html/search/variables_11.js [.f..tp...]`
- `docs/html/search/variables_12.html [.f..tp...]`
- `docs/html/search/variables_12.js [.f..tp...]`
- `docs/html/search/variables_13.html [.f..tp...]`
- `docs/html/search/variables_13.js [.f..tp...]`
- `docs/html/search/variables_2.html [.f..tp...]`
- `docs/html/search/variables_2.js [.f..tp...]`
- `docs/html/search/variables_3.html [.f..tp...]`
- `docs/html/search/variables_3.js [.f..tp...]`
- `docs/html/search/variables_4.html [.f..tp...]`
- `docs/html/search/variables_4.js [.f..tp...]`
- `docs/html/search/variables_5.html [.f..tp...]`
- `docs/html/search/variables_5.js [.f..tp...]`
- `docs/html/search/variables_6.html [.f..tp...]`
- `docs/html/search/variables_6.js [.f..tp...]`
- `docs/html/search/variables_7.html [.f..tp...]`
- `docs/html/search/variables_7.js [.f..tp...]`
- `docs/html/search/variables_8.html [.f..tp...]`
- `docs/html/search/variables_8.js [.f..tp...]`
- `docs/html/search/variables_9.html [.f..tp...]`
- `docs/html/search/variables_9.js [.f..tp...]`
- `docs/html/search/variables_a.html [.f..tp...]`
- `docs/html/search/variables_a.js [.f..tp...]`
- `docs/html/search/variables_b.html [.f..tp...]`
- `docs/html/search/variables_b.js [.f..tp...]`
- `docs/html/search/variables_c.html [.f..tp...]`
- `docs/html/search/variables_c.js [.f..tp...]`
- `docs/html/search/variables_d.html [.f..tp...]`
- `docs/html/search/variables_d.js [.f..tp...]`
- `docs/html/search/variables_e.html [.f..tp...]`
- `docs/html/search/variables_e.js [.f..tp...]`
- `docs/html/search/variables_f.html [.f..tp...]`
- `docs/html/search/variables_f.js [.f..tp...]`
- `helper_scripts/clat_complex.py [.f..tp...]`
- `helper_scripts/clat_on_il.py [.f..tp...]`
- `helper_scripts/clat_on_mem.py [.f..tp...]`
- `helper_scripts/clean_push.bsh [.f..t....]`
- `helper_scripts/combine_MABMruns.bsh [.f..tp...]`
- `helper_scripts/concatenate_files.py [.f..tp...]`
- `helper_scripts/gen_delete_inp.bsh [.f..t....]`
- `helper_scripts/gen_delete_puc.bsh [.f..t....]`
- `helper_scripts/get_crds_complex.py [.f..tp...]`
- `helper_scripts/pull_1col.py [.f..tp...]`
- `helper_scripts/pull_2col.py [.f..tp...]`
- `helper_scripts/pull_complex.py [.f..tp...]`
- `helper_scripts/pull_complex_all.py [.f..tp...]`
- `helper_scripts/pull_complex_allMABM.py [.f..tp...]`
- `helper_scripts/pull_complex_allMABM_time.py [.f..tp...]`
- `helper_scripts/remove_empty1.py [.f..tp...]`
- `helper_scripts/remove_empty2.py [.f..tp...]`
- `helper_scripts/remove_molecules_4inputs.py [.f..tp...]`
- `helper_scripts/remove_molecules_step1.py [.f..tp...]`
- `helper_scripts/text_search.bsh [.f..t....]`
- `helper_scripts/write_pdb_from_xyz.py [.f..tp...]`
- `include/error_handling.hpp [.f..tp...]`
- `include/json.hpp [.f..tp...]`
- `include/macro.hpp [.f..tp...]`
- `include/split.cpp [.f..tp...]`
- `include/tracing.hpp [.f..tp...]`
- `include/classes/class_Cluster.hpp [.f..tp...]`
- `include/classes/class_Coord.hpp [.f..tp...]`
- `include/classes/class_MDTimer.hpp [.f..t....]`
- `include/classes/class_Membrane.hpp [.f..tp...]`
- `include/classes/class_MolTemplate.hpp [.f..tp...]`
- `include/classes/class_Observable.hpp [.f..tp...]`
- `include/classes/class_Parameters.hpp [.f..tp...]`
- `include/classes/class_Quat.hpp [.f..tp...]`
- `include/classes/class_Rxns.hpp [.f..tp...]`
- `include/classes/class_SimulVolume.hpp [.f..tp...]`
- `include/classes/class_Vector.hpp [.f..tp...]`
- `include/classes/class_bngl_parser.hpp [.f..tp...]`
- `include/classes/class_copyCounters.hpp [.f..tp...]`
- `include/classes/mpi_functions.hpp [.f..tp...]`
- `include/debug/debug.hpp [.f..tp...]`
- `include/error/error.hpp [.f..tp...]`
- `include/io/io.hpp [.f..tp...]`
- `include/math/Faddeeva.hpp [.f..t....]`
- `include/math/constants.hpp [.f..tp...]`
- `include/math/math_functions.hpp [.f..tp...]`
- `include/math/matrix.hpp [.f..tp...]`
- `include/math/rand_gsl.hpp [.f..t....]`
- `include/mpi/mpi_function.hpp [.f..tp...]`
- `include/parser/parser_functions.hpp [.f..tp...]`
- `include/reactions/shared_reaction_functions.hpp [.f..tp...]`
- `include/reactions/association/association.hpp [.f..tp...]`
- `include/reactions/association/functions_for_spherical_system.hpp [.f..tp...]`
- `include/reactions/bimolecular/2D_reaction_table_functions.hpp [.f..tp...]`
- `include/reactions/bimolecular/bimolecular_reactions.hpp [.f..tp...]`
- `include/reactions/implicitlipid/implicitlipid_reactions.hpp [.f..tp...]`
- `include/reactions/unimolecular/unimolecular_reactions.hpp [.f..tp...]`
- `include/system_setup/system_setup.hpp [.f..tp...]`
- `include/trajectory_functions/trajectory_functions.hpp [.f..tp...]`
- `mpi_proto/complex.cpp [.f..tp...]`
- `mpi_proto/full_par_comm.cpp [.f..tp...]`
- `python_scripts/mergePDB.py [.f..tp...]`
- `python_scripts/plotPDB.py [.f..tp...]`
- `python_scripts/plot_copy_numbers_time.py [.f..tp...]`
- `run_code_tests/BimolecularTest3D/A.mol [.f..tp...]`
- `run_code_tests/BimolecularTest3D/BimolecularTest3D.ipynb [.f..tp...]`
- `run_code_tests/BimolecularTest3D/R.mol [.f..tp...]`
- `run_code_tests/BimolecularTest3D/parms3d.inp [.f..tp...]`
- `run_code_tests/ClathrinAssemblyFlat/flatClathrinAssemblyTests.ipynb [.f..tp...]`
- `run_code_tests/HomoTimerTest/A.mol [.f..tp...]`
- `run_code_tests/HomoTimerTest/parallel_debug_trimer-6traces.ipynb [.f..tp...]`
- `run_code_tests/HomoTimerTest/parms.inp [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/ParticleDiffusionOnFlatMembrane/A.mol [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/ParticleDiffusionOnFlatMembrane/IL.mol [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/ParticleDiffusionOnFlatMembrane/parms.inp [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/ParticleDiffusionOnFlatMembrane/singleParticlesPropagationOnFlatTest_previous.ipynb [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/SingleParticleBindingToImplicitLipid/A.mol [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/SingleParticleBindingToImplicitLipid/IL.mol [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/SingleParticleBindingToImplicitLipid/parms.inp [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/SingleParticleBindingToImplicitLipid/singleParticlesPropagationOnSphereTest.ipynb [.f..tp...]`
- `run_code_tests/PropagationOnSphereTest/SingleParticleBindingToImplicitLipid/singleParticlesPropagationOnSphereTest_previous.ipynb [.f..tp...]`
- `sample_inputs/#pucadyil#/IL.mol [.f..tp...]`
- `sample_inputs/#pucadyil#/Pucadyil_time_mean_std_sem [.f..tp...]`
- `sample_inputs/#pucadyil#/README [.f..tp...]`
- `sample_inputs/#pucadyil#/clat.mol [.f..tp...]`
- `sample_inputs/#pucadyil#/fitExpLagSimData_Time.py [.f..tp...]`
- `sample_inputs/#pucadyil#/fitPlot.m [.f..tp...]`
- `sample_inputs/#pucadyil#/get_clath_on_mem_time_from_histfile.py [.f..tp...]`
- `sample_inputs/#pucadyil#/myexp_lag.m [.f..tp...]`
- `sample_inputs/#pucadyil#/parms.inp [.f..tp...]`
- `sample_inputs/Gag/A.mol [.f..tp...]`
- `sample_inputs/Gag/nerdss [.f..t....]`
- `sample_inputs/Gag/parm.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/README [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/ReversibleBindingTests.pdf [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/collectData.sh [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/nerdss_reversibleBimolecular_3domains_WRedNERDSS.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/riccati.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_2D/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_2D/R.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_2D/README [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_2D/kon2D_value.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_2D/parms2D.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_2D/parms2D.inp~ [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D/R.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D/parms3d.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/M.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/parms3dto2d.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/.ipynb_checkpoints/parms3dto2d-checkpoint.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/ClathrinFlatLattice_Description.docx [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/ClathrinFlatLattice_Description.pdf [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/clat.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/clat8.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/clathModelNFsim_Nointra_Kd100uM_utl.bngl [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/clathModelNFsim_intra_Kd1_utl.bngl [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/design_pucker_crds.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/nerdss_clathrin [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/clathrin/parms_clath_kon0.2uM.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/parms_clath_kon100uM.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/parms_clath_kon1uM.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/plot_nerdss_traj.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/plot_nfsim_clath_only.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/time_vs_boundpairs_DeltatTwice_kd1_f5.9e6_VTHEORY.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/.ipynb_checkpoints/clat-checkpoint.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/.ipynb_checkpoints/parms_clath_kon0.2uM-checkpoint.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clathrin/.ipynb_checkpoints/parms_clath_kon100uM-checkpoint.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/ClockModel_DescriptionResults.docx [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/ClockModel_DescriptionResults.pdf [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/PrA.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/PrR.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/R.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/RNA.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/RNR.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/calc_peak_sepsAR.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/clock_model.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/clock_model.inp~ [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/findOscillationPeriodFFTZeroPad.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/clock_model/nerdss_clock_model [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/closed_homoTrimer/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/closed_homoTrimer/compare_ODE.ipynb [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/closed_homoTrimer/parmTri6.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/Phos.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/README [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/create.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/createDestroy_test.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/createDestroy_test.pdf [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/createDestroy_test_timeStep.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/.ipynb_checkpoints/A-checkpoint.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/.ipynb_checkpoints/create-checkpoint.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/create_destroy/.ipynb_checkpoints/createDestroy_test-checkpoint.pdf [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/C.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/hetTrimer/parm_autodiff_hetTri.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/hexamer/hex.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/hexamer/parms_phex.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/homoTrimer/parmTri6.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/implicit_lipid/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/implicit_lipid/IL.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/implicit_lipid/Implicit_Lipid model.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/implicit_lipid/Implicit_Lipid model_timeStep_48traces.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/implicit_lipid/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/implicit_lipid/parms.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/README_ModelDescriptions.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/ode_9species_memlocalize.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/solveKaeff_converge.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/ODEMatlab_LargeBox.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/IL_fixD_largbox_Avg_48traj.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/IL_fixD_largbox_SEM_48traj.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/M.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/nerdss_IL [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/parms.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/ODEMatlab_SmallBox.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/CompELILEQABM.tif [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/M.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/NERDSS_box752_fastD_64traj_AvgsEL.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/NERDSS_box752_fastD_64traj_SEM.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/parms.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/M.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/nerdss_EL [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/parms.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/M.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/nerdss_IL [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/parms.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/FPR_slowD_Avgs_48traj.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/FPR_slowD_SEM_48traj.txt [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/M.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/nerdss_SlowDsol [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/parms.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/E.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/MMmodel_Rxn2.eps [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/MMmodel_Rxn2.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/README [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/S.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/michaelis.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/michaelis.inp~ [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/nerdss_michaelis_menten [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/michaelis_menten/ode_michaelisMenten.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/sphere/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/sphere/IL.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/sphere/parms_sphere.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/sphere/sphere.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/trimer/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/trimer/Analyze_Trimer.ipynb [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/trimer/B.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/trimer/C.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/trimer/parm_autodiff_hetTri.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/trimer/parm_hetTri.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/A.mol [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/Avstime_vsTheory.fig [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/OUTPUT [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/UnimolecularReactions.pdf [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/nerdss_unimolecular_reverse [.f..t....]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/reversible_unimolecular.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/reversible_unimolecular.m~ [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/uni_reverse.inp [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/.ipynb_checkpoints/reversible_unimolecular-checkpoint.m [.f..tp...]`
- `sample_inputs/VALIDATE_SUITE/unimolecular_reverse/.ipynb_checkpoints/uni_reverse-checkpoint.inp [.f..tp...]`
- `sample_inputs/auto_phos/A.mol [.f..tp...]`
- `sample_inputs/auto_phos/Phos.mol [.f..tp...]`
- `sample_inputs/auto_phos/autophos_D10.inp [.f..tp...]`
- `sample_inputs/clathrin_coat/cooperative_flat_clat.dir/00ReadMe.md [.f..tp...]`
- `sample_inputs/clathrin_coat/cooperative_flat_clat.dir/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_coat/cooperative_flat_clat.dir/parms_kd100_10x.inp [.f..t....]`
- `sample_inputs/clathrin_coat/cooperative_flat_clat.dir/parms_kd100_coop_none.inp [.f..t....]`
- `sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir/nerdss -> /home/local/WIN/msang2/mankun/nerdss_development/bin/nerdss [.L...p...]`
- `sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir/parms_clath_ap_pip.inp [.f..tp...]`
- `sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir/pip2.mol [.f..tp...]`
- `sample_inputs/clathrin_coat/flat_clat.dir/.clath_kon1uM.inp.swp [.f..tp...]`
- `sample_inputs/clathrin_coat/flat_clat.dir/00ReadMe.md [.f..tp...]`
- `sample_inputs/clathrin_coat/flat_clat.dir/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_coat/flat_clat.dir/clath_kon1uM.inp [.f..tp...]`
- `sample_inputs/clathrin_coat/puckered_clat.dir/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_coat/puckered_clat.dir/nerdss [.f..t....]`
- `sample_inputs/clathrin_coat/puckered_clat.dir/parms_pucker.inp [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/curved_clathrin_solution/README [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/curved_clathrin_solution/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/curved_clathrin_solution/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/curved_clathrin_solution/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invitro/IL.mol [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invitro/README [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invitro/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invitro/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invivo/README [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invivo/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invivo/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invivo/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invivo/parmsMacroRate.inp [.f..tp...]`
- `sample_inputs/clathrin_invitro_invivo/invivo/pip2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/README.md [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2.mol.old [.f..tp...]`
- `sample_inputs/clathrin_pioneer/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/e.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/f.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/p.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/parmsMacroRate.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/parmscross.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/parmsfel.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/parmsffl.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2cla/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2cla/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2cla/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2cla/pip2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafcho/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafcho/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafcho/fcho.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafcho/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafcho/pip2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps/eps.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps/fcho.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps/pip2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps4site/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps4site/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps4site/eps.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps4site/fcho.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps4site/parms.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ap2clafchoeps4site/pip2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fel/e.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fel/f.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fel/p.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fel/parmsfel.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fela/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fela/e.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fela/f.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fela/p.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/fela/parmsfela.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ffl/f.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ffl/p.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/ffl/parmsffl.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/figs/Model.png [.f..tp...]`
- `sample_inputs/clathrin_pioneer/macro/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/macro/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/macro/e.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/macro/f.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/macro/parmsMacroRate.inp [.f..tp...]`
- `sample_inputs/clathrin_pioneer/macro/pip2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/test/ap2.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/test/clat.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/test/e.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/test/f.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/test/p.mol [.f..tp...]`
- `sample_inputs/clathrin_pioneer/test/parmscross.inp [.f..tp...]`
- `sample_inputs/compartment/.histogram_complexes_time.dat.swp [.f..tp...]`
- `sample_inputs/compartment/clat.mol [.f..tp...]`
- `sample_inputs/compartment/clatComp.mol [.f..tp...]`
- `sample_inputs/compartment/clath_compartment.inp [.f..tp...]`
- `sample_inputs/compartment/parms.inp~ [.f..tp...]`
- `sample_inputs/enzyme/README [.f..tp...]`
- `sample_inputs/enzyme/ap2.mol [.f..tp...]`
- `sample_inputs/enzyme/clat.mol [.f..tp...]`
- `sample_inputs/enzyme/parms_clat_enzyme.inp [.f..tp...]`
- `sample_inputs/enzyme/pip2.mol [.f..tp...]`
- `sample_inputs/enzyme/syn.mol [.f..tp...]`
- `sample_inputs/gag_coat/gagpol/gag.mol [.f..tp...]`
- `sample_inputs/gag_coat/gagpol/parms_gag_pol.inp [.f..tp...]`
- `sample_inputs/gag_coat/gagpol/pol.mol [.f..tp...]`
- `sample_inputs/gag_coat/solution/gag.mol [.f..tp...]`
- `sample_inputs/gag_coat/solution/parms.inp [.f..tp...]`
- `sample_inputs/gagsphere/gag.mol [.f..tp...]`
- `sample_inputs/gagsphere/gagOriginalModelSolution.inp [.f..tp...]`
- `sample_inputs/gagsphere/gagOriginalModelSolution.inp~ [.f..tp...]`
- `sample_inputs/genetic_oscillator/A.mol [.f..tp...]`
- `sample_inputs/genetic_oscillator/PrA.mol [.f..tp...]`
- `sample_inputs/genetic_oscillator/PrR.mol [.f..tp...]`
- `sample_inputs/genetic_oscillator/R.mol [.f..tp...]`
- `sample_inputs/genetic_oscillator/README [.f..tp...]`
- `sample_inputs/genetic_oscillator/RNA.mol [.f..tp...]`
- `sample_inputs/genetic_oscillator/RNR.mol [.f..tp...]`
- `sample_inputs/genetic_oscillator/clock_model.inp [.f..tp...]`
- `sample_inputs/genetic_oscillator/clock_model.inp~ [.f..tp...]`
- `sample_inputs/implicit_lipid/B.mol [.f..tp...]`
- `sample_inputs/implicit_lipid/IL.mol [.f..tp...]`
- `sample_inputs/implicit_lipid/parms.inp [.f..tp...]`
- `sample_inputs/membrane_rev_localization/README_ModelDescriptions.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/ode_9species_memlocalize.m [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/ODEMatlab_LargeBox.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/explicit/A.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/explicit/AvgsMemloc_FPRlargbox_10traj.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/explicit/B.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/explicit/M.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/explicit/SEMMemloc_FPRlargbox_10traj.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/explicit/parms.inp [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/implicitLipid/A.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/implicitLipid/B.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/implicitLipid/IL_fixD_largbox_Avg_48traj.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/implicitLipid/IL_fixD_largbox_SEM_48traj.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/implicitLipid/M.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/LargeBox/implicitLipid/parms.inp [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/ODEMatlab_SmallBox.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/A.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/B.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/CompELILEQABM.tif [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/M.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/NERDSS_box752_fastD_64traj_Avgs.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/NERDSS_box752_fastD_64traj_SEM.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/parms.inp [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/EL/A.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/EL/B.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/EL/M.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/EL/parms.inp [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/IL/A.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/IL/B.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/IL/M.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/FastDsol/IL/parms.inp [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/A.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/B.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/FPR_slowD_Avgs_48traj.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/FPR_slowD_SEM_48traj.txt [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/M.mol [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/nerdss [.f..t....]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/output [.f..tp...]`
- `sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/parms.inp [.f..tp...]`
- `sample_inputs/sphere/A.mol [.f..tp...]`
- `sample_inputs/sphere/IL.mol [.f..tp...]`
- `sample_inputs/sphere/parms_sphere.inp [.f..tp...]`
- `sample_inputs/testAdd/A.mol [.f..tp...]`
- `sample_inputs/testAdd/B.mol [.f..tp...]`
- `sample_inputs/testAdd/M.mol [.f..tp...]`
- `sample_inputs/testAdd/add.inp [.f..tp...]`
- `sample_inputs/testAdd/parms.inp [.f..tp...]`
- `sample_inputs/testAdd/box/A.mol [.f..tp...]`
- `sample_inputs/testAdd/box/B.mol [.f..tp...]`
- `sample_inputs/testAdd/box/M.mol [.f..tp...]`
- `sample_inputs/testAdd/box/add.inp [.f..tp...]`
- `sample_inputs/testAdd/box/parms.inp [.f..tp...]`
- `sample_inputs/testAdd/box/1/A.mol [.f..tp...]`
- `sample_inputs/testAdd/box/1/B.mol [.f..tp...]`
- `sample_inputs/testAdd/box/1/M.mol [.f..tp...]`
- `sample_inputs/testAdd/box/1/add.inp [.f..tp...]`
- `sample_inputs/testAdd/box/1/parms.inp [.f..tp...]`
- `sample_inputs/testAdd/sphere/A.mol [.f..tp...]`
- `sample_inputs/testAdd/sphere/B.mol [.f..tp...]`
- `sample_inputs/testAdd/sphere/M.mol [.f..tp...]`
- `sample_inputs/testAdd/sphere/add.inp [.f..tp...]`
- `sample_inputs/testAdd/sphere/parms.inp [.f..tp...]`
- `sample_inputs/testAdd/sphere/1/A.mol [.f..tp...]`
- `sample_inputs/testAdd/sphere/1/B.mol [.f..tp...]`
- `sample_inputs/testAdd/sphere/1/M.mol [.f..tp...]`
- `sample_inputs/testAdd/sphere/1/add.inp [.f..tp...]`
- `sample_inputs/testAdd/sphere/1/parms.inp [.f..tp...]`
- `src/boundary_conditions/check_if_spans.cpp [.f..tp...]`
- `src/boundary_conditions/check_if_spans_box.cpp [.f..tp...]`
- `src/boundary_conditions/check_if_spans_sphere.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_complex_compartment.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_complex_rad_rot.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_complex_rad_rot_box.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_complex_rad_rot_sphere.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_check_span.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_check_span_box.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_check_span_sphere.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_complex_compartment.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_complex_rad_rot.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_complex_rad_rot_box.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_complex_rad_rot_nocheck.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_complex_rad_rot_nocheck_box.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_complex_rad_rot_nocheck_sphere.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_complex_rad_rot_sphere.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_tmp_crds.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_tmp_crds_box.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_tmp_crds_compartment.cpp [.f..tp...]`
- `src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp [.f..tp...]`
- `src/classes/class_Cluster.cpp [.f..tp...]`
- `src/classes/class_Coord.cpp [.f..tp...]`
- `src/classes/class_MDTimer.cpp [.f..t....]`
- `src/classes/class_MolTemplate.cpp [.f..tp...]`
- `src/classes/class_Observable.cpp [.f..tp...]`
- `src/classes/class_Quat.cpp [.f..tp...]`
- `src/classes/class_Rxns.cpp [.f..tp...]`
- `src/classes/class_SimulVolume.cpp [.f..tp...]`
- `src/classes/class_bngl_parser_functions.cpp [.f..tp...]`
- `src/classes/mpi_functions.cpp [.f..tp...]`
- `src/debug/debug.cpp [.f..tp...]`
- `src/debug/debug_molecule_complex_missmatch.cpp [.f..tp...]`
- `src/debug/write_debug_information.cpp [.f..tp...]`
- `src/error/error.cpp [.f..tp...]`
- `src/io/eye_candy.cpp [.f..tp...]`
- `src/io/init_NboundPairs.cpp [.f..tp...]`
- `src/io/init_counterCopyNums.cpp [.f..tp...]`
- `src/io/init_print_dimers.cpp [.f..tp...]`
- `src/io/init_speciesFile.cpp [.f..tp...]`
- `src/io/print_complex_hist.cpp [.f..tp...]`
- `src/io/print_dimers.cpp [.f..tp...]`
- `src/io/print_system_information.cpp [.f..tp...]`
- `src/io/write_NboundPairs.cpp [.f..tp...]`
- `src/io/write_all_species.cpp [.f..tp...]`
- `src/io/write_bonded_complex_json.cpp [.f..tp...]`
- `src/io/write_complex_components.cpp [.f..tp...]`
- `src/io/write_complex_crds.cpp [.f..tp...]`
- `src/io/write_crds.cpp [.f..tp...]`
- `src/io/write_observables.cpp [.f..tp...]`
- `src/io/write_pdb.cpp [.f..tp...]`
- `src/io/write_psf.cpp [.f..tp...]`
- `src/io/write_timestep_information.cpp [.f..tp...]`
- `src/io/write_traj.cpp [.f..tp...]`
- `src/io/write_transition.cpp [.f..tp...]`
- `src/io/write_xyz.cpp [.f..tp...]`
- `src/io/write_xyz_assoc.cpp [.f..tp...]`
- `src/io/write_xyz_assoc_cout.cpp [.f..tp...]`
- `src/io_mpi/merge_outputs.cpp [.f..tp...]`
- `src/io_mpi/print_complex_hist.cpp [.f..tp...]`
- `src/io_mpi/write_all_species.cpp [.f..tp...]`
- `src/io_mpi/write_output.cpp [.f..tp...]`
- `src/io_mpi/write_pdb.cpp [.f..tp...]`
- `src/math/Faddeeva.cpp [.f..t....]`
- `src/math/math_functions.cpp [.f..tp...]`
- `src/math/matrix_functions.cpp [.f..tp...]`
- `src/math/rand_gsl.cpp [.f..tp...]`
- `src/mpi/check_molecule_coordinates.cpp [.f..tp...]`
- `src/mpi/delete_disappeared.cpp [.f..tp...]`
- `src/mpi/deserialize.cpp [.f..tp...]`
- `src/mpi/id_index.cpp [.f..tp...]`
- `src/mpi/prepare.cpp [.f..tp...]`
- `src/mpi/send_receive.cpp [.f..tp...]`
- `src/mpi/serialize.cpp [.f..tp...]`
- `src/parser/areSameExceptState.cpp [.f..tp...]`
- `src/parser/check_for_state_change.cpp [.f..tp...]`
- `src/parser/check_for_valid_states.cpp [.f..tp...]`
- `src/parser/create_conjugate_reaction_itrs.cpp [.f..tp...]`
- `src/parser/determine_bound_iface_index.cpp [.f..tp...]`
- `src/parser/determine_iface_indices.cpp [.f..tp...]`
- `src/parser/display_all_MolTemplates.cpp [.f..tp...]`
- `src/parser/display_all_reactions.cpp [.f..tp...]`
- `src/parser/parse_command.cpp [.f..tp...]`
- `src/parser/parse_input.cpp [.f..tp...]`
- `src/parser/parse_input_array.cpp [.f..tp...]`
- `src/parser/parse_input_for_a_new_simulation.cpp [.f..tp...]`
- `src/parser/parse_input_for_a_restart_simulation.cpp [.f..tp...]`
- `src/parser/parse_input_for_add_file.cpp [.f..tp...]`
- `src/parser/parse_molFile.cpp [.f..tp...]`
- `src/parser/parse_molecule_bngl.cpp [.f..tp...]`
- `src/parser/parse_observable.cpp [.f..tp...]`
- `src/parser/parse_reaction.cpp [.f..tp...]`
- `src/parser/parse_states.cpp [.f..tp...]`
- `src/parser/populate_reaction_lists.cpp [.f..tp...]`
- `src/parser/read_bonds.cpp [.f..tp...]`
- `src/parser/read_boolean.cpp [.f..tp...]`
- `src/parser/read_internal_coordinates.cpp [.f..tp...]`
- `src/parser/remove_comment.cpp [.f..tp...]`
- `src/parser/write_mol_iface.cpp [.f..tp...]`
- `src/reactions/DDpirr_pfree_ratio_ps.cpp [.f..tp...]`
- `src/reactions/angleSignIsCorrect.cpp [.f..tp...]`
- `src/reactions/areParallel.cpp [.f..tp...]`
- `src/reactions/areSameAngle.cpp [.f..tp...]`
- `src/reactions/associate.cpp [.f..tp...]`
- `src/reactions/associate_ImplicitLipid.cpp [.f..tp...]`
- `src/reactions/associate_ImplicitLipid_box.cpp [.f..tp...]`
- `src/reactions/associate_ImplicitLipid_sphere.cpp [.f..tp...]`
- `src/reactions/associate_sphere.cpp [.f..tp...]`
- `src/reactions/break_interaction.cpp [.f..tp...]`
- `src/reactions/break_interaction_implicitlipid.cpp [.f..tp...]`
- `src/reactions/calc_angular_displacement.cpp [.f..tp...]`
- `src/reactions/calc_one_angular_displacement.cpp [.f..tp...]`
- `src/reactions/calc_pirr.cpp [.f..tp...]`
- `src/reactions/calculate_omega.cpp [.f..tp...]`
- `src/reactions/calculate_phi.cpp [.f..tp...]`
- `src/reactions/check_bases.cpp [.f..tp...]`
- `src/reactions/check_bimolecular_reactions.cpp [.f..tp...]`
- `src/reactions/check_compartment_reaction.cpp [.f..tp...]`
- `src/reactions/check_dissociation.cpp [.f..tp...]`
- `src/reactions/check_dissociation_implicitlipid.cpp [.f..tp...]`
- `src/reactions/check_for_structure_overlap.cpp [.f..tp...]`
- `src/reactions/check_for_structure_overlap_system.cpp [.f..tp...]`
- `src/reactions/check_for_unimolecular_reactions.cpp [.f..tp...]`
- `src/reactions/check_for_unimolstatechange_reactions.cpp [.f..tp...]`
- `src/reactions/check_for_zeroth_order_creation.cpp [.f..tp...]`
- `src/reactions/check_implicit_reactions.cpp [.f..tp...]`
- `src/reactions/check_overlap.cpp [.f..tp...]`
- `src/reactions/check_perform_zeroth_first_order_reactions.cpp [.f..tp...]`
- `src/reactions/com_of_two_tmp_complexes.cpp [.f..tp...]`
- `src/reactions/conservedMags.cpp [.f..tp...]`
- `src/reactions/conservedRigid.cpp [.f..tp...]`
- `src/reactions/correct_structutre.cpp [.f..tp...]`
- `src/reactions/create_DDMatrices.cpp [.f..tp...]`
- `src/reactions/create_arbitrary_vector.cpp [.f..tp...]`
- `src/reactions/create_molecule_and_complex_from_rxn.cpp [.f..tp...]`
- `src/reactions/create_molecule_and_complex_from_transmission_rxn.cpp [.f..tp...]`
- `src/reactions/create_normMatrix.cpp [.f..tp...]`
- `src/reactions/create_pirMatrix.cpp [.f..tp...]`
- `src/reactions/create_survMatrix.cpp [.f..tp...]`
- `src/reactions/determine_2D_bimolecular_reaction_probability.cpp [.f..tp...]`
- `src/reactions/determine_2D_implicitlipid_reaction_probability.cpp [.f..tp...]`
- `src/reactions/determine_3D_bimolecular_reaction_probability.cpp [.f..tp...]`
- `src/reactions/determine_3D_implicitlipid_reaction_probability.cpp [.f..tp...]`
- `src/reactions/determine_entering_compartment_probability.cpp [.f..tp...]`
- `src/reactions/determine_exiting_compartment_probability.cpp [.f..tp...]`
- `src/reactions/determine_if_reaction_occurs.cpp [.f..tp...]`
- `src/reactions/determine_normal.cpp [.f..tp...]`
- `src/reactions/determine_parent_complex.cpp [.f..tp...]`
- `src/reactions/determine_parent_complex_IL.cpp [.f..tp...]`
- `src/reactions/determine_rotation_angles.cpp [.f..tp...]`
- `src/reactions/evaluate_binding_within_complex.cpp [.f..tp...]`
- `src/reactions/find_reaction_rate_state.cpp [.f..tp...]`
- `src/reactions/find_which_reaction.cpp [.f..tp...]`
- `src/reactions/find_which_state_change_reaction.cpp [.f..tp...]`
- `src/reactions/functions_for_spherical_system.cpp [.f..tp...]`
- `src/reactions/functions_implicitlipid.cpp [.f..tp...]`
- `src/reactions/get_distance.cpp [.f..tp...]`
- `src/reactions/get_distance_to_surface.cpp [.f..tp...]`
- `src/reactions/get_geodesic_distance.cpp [.f..tp...]`
- `src/reactions/get_prevNorm.cpp [.f..tp...]`
- `src/reactions/get_prevSurv.cpp [.f..tp...]`
- `src/reactions/hasIntangibles.cpp [.f..tp...]`
- `src/reactions/init_association_events.cpp [.f..tp...]`
- `src/reactions/initialize_molecule_after_reaction.cpp [.f..tp...]`
- `src/reactions/initialize_molecule_after_transmission_reaction.cpp [.f..tp...]`
- `src/reactions/integrator.cpp [.f..tp...]`
- `src/reactions/isReactant.cpp [.f..tp...]`
- `src/reactions/measure_complex_displacement.cpp [.f..tp...]`
- `src/reactions/measure_overlap_free_protein_interfaces.cpp [.f..tp...]`
- `src/reactions/measure_overlap_protein_interfaces.cpp [.f..tp...]`
- `src/reactions/measure_separations_to_identify_possible_reactions.cpp [.f..tp...]`
- `src/reactions/norm_function.cpp [.f..tp...]`
- `src/reactions/omega_rotation.cpp [.f..tp...]`
- `src/reactions/orient_crds_to_template.cpp [.f..tp...]`
- `src/reactions/passocF.cpp [.f..tp...]`
- `src/reactions/passocF_1D.cpp [.f..tp...]`
- `src/reactions/perform_bimolecular_reactions.cpp [.f..tp...]`
- `src/reactions/perform_bimolecular_state_change.cpp [.f..tp...]`
- `src/reactions/perform_bimolecular_state_change_box.cpp [.f..tp...]`
- `src/reactions/perform_bimolecular_state_change_sphere.cpp [.f..tp...]`
- `src/reactions/perform_implicitlipid_state_change.cpp [.f..tp...]`
- `src/reactions/perform_implicitlipid_state_change_box.cpp [.f..tp...]`
- `src/reactions/perform_implicitlipid_state_change_sphere.cpp [.f..tp...]`
- `src/reactions/perform_transmission_reaction.cpp [.f..tp...]`
- `src/reactions/phi_rotation.cpp [.f..tp...]`
- `src/reactions/pir_function.cpp [.f..tp...]`
- `src/reactions/pirr_pfree_ratio_psF.cpp [.f..tp...]`
- `src/reactions/print_association_events.cpp [.f..tp...]`
- `src/reactions/remove_empty_slots.cpp [.f..tp...]`
- `src/reactions/requiresSignFlip.cpp [.f..tp...]`
- `src/reactions/reverse_rotation.cpp [.f..tp...]`
- `src/reactions/rotate.cpp [.f..tp...]`
- `src/reactions/save_mem_orientation.cpp [.f..tp...]`
- `src/reactions/size_lookup.cpp [.f..tp...]`
- `src/reactions/subtract_off_com_position.cpp [.f..tp...]`
- `src/reactions/survival_function.cpp [.f..tp...]`
- `src/reactions/theta_rotation.cpp [.f..tp...]`
- `src/reactions/track_association_events.cpp [.f..tp...]`
- `src/reactions/transform.cpp [.f..tp...]`
- `src/reactions/update_Nboundpairs.cpp [.f..tp...]`
- `src/reactions/update_complex_tmp_com_crds.cpp [.f..tp...]`
- `src/system_setup/areInVicinity.cpp [.f..tp...]`
- `src/system_setup/are_molecules_in_vicinity.cpp [.f..tp...]`
- `src/system_setup/determine_shape_molecule.cpp [.f..tp...]`
- `src/system_setup/generate_coordinates.cpp [.f..tp...]`
- `src/system_setup/generate_coordinates_for_restart.cpp [.f..tp...]`
- `src/system_setup/initialize_complex.cpp [.f..tp...]`
- `src/system_setup/initialize_molecule.cpp [.f..tp...]`
- `src/system_setup/initialize_parameters_for_implicitlipid_and_compartment_model.cpp [.f..tp...]`
- `src/system_setup/initialize_states.cpp [.f..tp...]`
- `src/system_setup/set_exclude_volume_bound.cpp [.f..tp...]`
- `src/system_setup/set_rMaxLimit.cpp [.f..tp...]`
- `src/trajectory_functions/clear_reweight_vecs.cpp [.f..tp...]`
- `src/trajectory_functions/create_complex_propagation_vectors_on_sphere.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_box.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_fiber.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_fiber_box.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_memtest.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_memtest_box.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_memtest_cluster.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_memtest_cluster_box.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_memtest_cluster_sphere.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_memtest_sphere.cpp [.f..tp...]`
- `src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp [.f..tp...]`
- `third_party/nlohmann_json/LICENSE [.f..tp...]`

## Unchanged Files

- `sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir/nerdss`

## Directory Entries From Dry Run

- `deleted directory: obj/trajectory_functions/`
- `deleted directory: obj/system_setup/`
- `deleted directory: obj/reactions/`
- `deleted directory: obj/parser/`
- `deleted directory: obj/math/`
- `deleted directory: obj/io/`
- `deleted directory: obj/error/`
- `deleted directory: obj/classes/`
- `deleted directory: obj/boundary_conditions/`
- `deleted directory: obj/`
- `deleted directory: images/`
- `deleted directory: bin/`
- `deleted directory: assets/`
- `deleted directory: RESTARTS/`
- `deleted directory: PDB/`
- `deleted directory: sample_inputs/VALIDATE_SUITE/hetTrimer/RESTARTS/`
- `deleted directory: sample_inputs/VALIDATE_SUITE/hetTrimer/PDB/`
- `deleted directory: sample_inputs/VALIDATE_SUITE/homoTrimer/RESTARTS/`
- `deleted directory: sample_inputs/VALIDATE_SUITE/homoTrimer/PDB/`
- `deleted directory: sample_inputs/VALIDATE_SUITE/trimer/RESTARTS/`
- `deleted directory: sample_inputs/VALIDATE_SUITE/trimer/PDB/`
- `directory metadata: ./`
- `directory metadata: .cache/`
- `directory metadata: .cache/clangd/`
- `directory metadata: .cache/clangd/index/`
- `directory metadata: .github/`
- `directory metadata: .github/workflows/`
- `directory metadata: EXEs/`
- `directory metadata: MATLAB_functions/`
- `directory metadata: angle_calculation/`
- `directory metadata: docs/`
- `directory metadata: docs/associate_figs/`
- `directory metadata: docs/associate_figs/angles/`
- `directory metadata: docs/examples/`
- `directory metadata: docs/examples/irrev_flat_clat/`
- `directory metadata: docs/examples/zero_create/`
- `directory metadata: docs/html/`
- `directory metadata: docs/html/search/`
- `directory metadata: helper_scripts/`
- `directory metadata: include/`
- `directory metadata: include/boundary_conditions/`
- `directory metadata: include/classes/`
- `directory metadata: include/debug/`
- `directory metadata: include/error/`
- `directory metadata: include/io/`
- `directory metadata: include/math/`
- `directory metadata: include/mpi/`
- `directory metadata: include/parser/`
- `directory metadata: include/reactions/`
- `directory metadata: include/reactions/association/`
- `directory metadata: include/reactions/bimolecular/`
- `directory metadata: include/reactions/implicitlipid/`
- `directory metadata: include/reactions/unimolecular/`
- `directory metadata: include/system_setup/`
- `directory metadata: include/trajectory_functions/`
- `directory metadata: mpi_proto/`
- `directory metadata: python_scripts/`
- `directory metadata: run_code_tests/`
- `directory metadata: run_code_tests/BimolecularTest3D/`
- `directory metadata: run_code_tests/ClathrinAssemblyFlat/`
- `directory metadata: run_code_tests/HomoTimerTest/`
- `directory metadata: run_code_tests/PropagationOnSphereTest/`
- `directory metadata: run_code_tests/PropagationOnSphereTest/ParticleDiffusionOnFlatMembrane/`
- `directory metadata: run_code_tests/PropagationOnSphereTest/SingleParticleBindingToImplicitLipid/`
- `directory metadata: sample_inputs/`
- `directory metadata: sample_inputs/#pucadyil#/`
- `directory metadata: sample_inputs/Gag/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/bimolecular_reversible/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_2D/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/.ipynb_checkpoints/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/clathrin/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/clathrin/.ipynb_checkpoints/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/clock_model/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/closed_homoTrimer/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/create_destroy/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/create_destroy/.ipynb_checkpoints/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/hetTrimer/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/hexamer/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/homoTrimer/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/implicit_lipid/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/mem_localization/SmallBox/SlowDsol/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/michaelis_menten/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/sphere/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/trimer/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/unimolecular_reverse/`
- `directory metadata: sample_inputs/VALIDATE_SUITE/unimolecular_reverse/.ipynb_checkpoints/`
- `directory metadata: sample_inputs/auto_phos/`
- `directory metadata: sample_inputs/clathrin_coat/`
- `directory metadata: sample_inputs/clathrin_coat/cooperative_flat_clat.dir/`
- `directory metadata: sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir/`
- `directory metadata: sample_inputs/clathrin_coat/flat_clat.dir/`
- `directory metadata: sample_inputs/clathrin_coat/puckered_clat.dir/`
- `directory metadata: sample_inputs/clathrin_invitro_invivo/`
- `directory metadata: sample_inputs/clathrin_invitro_invivo/curved_clathrin_solution/`
- `directory metadata: sample_inputs/clathrin_invitro_invivo/invitro/`
- `directory metadata: sample_inputs/clathrin_invitro_invivo/invivo/`
- `directory metadata: sample_inputs/clathrin_pioneer/`
- `directory metadata: sample_inputs/clathrin_pioneer/ap2cla/`
- `directory metadata: sample_inputs/clathrin_pioneer/ap2clafcho/`
- `directory metadata: sample_inputs/clathrin_pioneer/ap2clafchoeps/`
- `directory metadata: sample_inputs/clathrin_pioneer/ap2clafchoeps4site/`
- `directory metadata: sample_inputs/clathrin_pioneer/fel/`
- `directory metadata: sample_inputs/clathrin_pioneer/fela/`
- `directory metadata: sample_inputs/clathrin_pioneer/ffl/`
- `directory metadata: sample_inputs/clathrin_pioneer/figs/`
- `directory metadata: sample_inputs/clathrin_pioneer/macro/`
- `directory metadata: sample_inputs/clathrin_pioneer/test/`
- `directory metadata: sample_inputs/compartment/`
- `directory metadata: sample_inputs/enzyme/`
- `directory metadata: sample_inputs/gag_coat/`
- `directory metadata: sample_inputs/gag_coat/gagpol/`
- `directory metadata: sample_inputs/gag_coat/solution/`
- `directory metadata: sample_inputs/gagsphere/`
- `directory metadata: sample_inputs/genetic_oscillator/`
- `directory metadata: sample_inputs/implicit_lipid/`
- `directory metadata: sample_inputs/membrane_rev_localization/`
- `directory metadata: sample_inputs/membrane_rev_localization/LargeBox/`
- `directory metadata: sample_inputs/membrane_rev_localization/LargeBox/explicit/`
- `directory metadata: sample_inputs/membrane_rev_localization/LargeBox/implicitLipid/`
- `directory metadata: sample_inputs/membrane_rev_localization/SmallBox/`
- `directory metadata: sample_inputs/membrane_rev_localization/SmallBox/FastDsol/`
- `directory metadata: sample_inputs/membrane_rev_localization/SmallBox/FastDsol/EL/`
- `directory metadata: sample_inputs/membrane_rev_localization/SmallBox/FastDsol/IL/`
- `directory metadata: sample_inputs/membrane_rev_localization/SmallBox/SlowDsol/`
- `directory metadata: sample_inputs/sphere/`
- `directory metadata: sample_inputs/testAdd/`
- `directory metadata: sample_inputs/testAdd/box/`
- `directory metadata: sample_inputs/testAdd/box/1/`
- `directory metadata: sample_inputs/testAdd/sphere/`
- `directory metadata: sample_inputs/testAdd/sphere/1/`
- `directory metadata: src/`
- `directory metadata: src/boundary_conditions/`
- `directory metadata: src/classes/`
- `directory metadata: src/debug/`
- `directory metadata: src/error/`
- `directory metadata: src/io/`
- `directory metadata: src/io_mpi/`
- `directory metadata: src/math/`
- `directory metadata: src/mpi/`
- `directory metadata: src/parser/`
- `directory metadata: src/reactions/`
- `directory metadata: src/system_setup/`
- `directory metadata: src/trajectory_functions/`
- `directory metadata: third_party/`
- `directory metadata: third_party/nlohmann_json/`
