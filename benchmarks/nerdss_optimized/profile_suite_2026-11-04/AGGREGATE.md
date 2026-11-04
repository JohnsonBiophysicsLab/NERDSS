# Aggregated profile over 22 cases

Unweighted mean over cases. "top5" = number of cases in which the function is among the five largest by self time.

| function | mean self % | mean inclusive % | max self % (case) | top5 |
| --- | --- | --- | --- | --- |
| `measure_separations_to_identify_possible_reactions` | 11.7 | 37.1 | 25.4 (scale_10k) | 20 |
| `std::vector::operator[] (inline)` | 5.5 | 5.5 | 9.0 (unimol_state_change) | 17 |
| `check_bimolecular_reactions` | 5.0 | 20.8 | 18.4 (hetTrimer) | 10 |
| `Complex::update_properties` | 4.2 | 6.5 | 9.4 (implicit_lipid) | 10 |
| `create_complex_propagation_vectors` | 3.7 | 10.8 | 7.6 (rev_3D) | 7 |
| `gsl_ran_gaussian_ziggurat` | 3.5 | 6.0 | 8.0 (hist_mol_states) | 8 |
| `gsl_ran_wishart` | 2.9 | 2.9 | 7.4 (hist_mol_states) | 2 |
| `std::vector::size (inline)` | 2.8 | 2.8 | 6.3 (hetTrimer) | 3 |
| `__ieee754_exp_fma` | 2.6 | 2.6 | 7.7 (homoTrimer) | 4 |
| `main` | 2.5 | 100.0 | 5.5 (implicit_lipid) | 2 |
| `std::vector::push_back (inline)` | 2.4 | 4.4 | 14.4 (scale_10k) | 4 |
| `do_sin` | 1.8 | 1.8 | 6.4 (hexamer) | 3 |
| `operator- (inline)` | 1.6 | 1.6 | 3.8 (hetTrimer) | 1 |
| `__sincos_fma` | 1.5 | 1.5 | 3.3 (hexamer) | 0 |
| `Complex::propagate` | 1.5 | 15.5 | 3.0 (enzyme) | 0 |
| `SimulVolume::update_memberMolLists` | 1.3 | 6.2 | 2.9 (enzyme) | 0 |
| `QuatRotation::rotate (inline)` | 1.3 | 2.2 | 10.1 (clathrin) | 2 |
| `determine_3D_bimolecular_reaction_probability` | 1.2 | 8.4 | 5.2 (closed_homoTrimer) | 1 |
| `operator+ (inline)` | 1.2 | 1.2 | 3.3 (scale_2k) | 0 |
| `create_pirMatrix` | 1.2 | 1.4 | 16.1 (cluster_gagsphere) | 2 |
| `::sweep_separation_box [clone .constprop.1]` | 1.1 | 6.0 | 4.7 (homoTrimer) | 1 |
| `__erfc` | 1.1 | 1.9 | 4.1 (trimer) | 0 |
| `reflect_traj_complex_rad_rot_box` | 1.0 | 3.5 | 2.1 (implicit_lipid) | 0 |
| `std::_Vector_base::_Vector_impl_data::_M_copy_data (inline)` | 0.9 | 0.9 | 2.7 (cluster_mem_loc) | 0 |
| `Quat::operator* (inline)` | 0.9 | 0.9 | 6.0 (compartment) | 1 |
| `std::vector::_M_erase_at_end (inline)` | 0.9 | 1.0 | 4.1 (scale_10k) | 0 |
| `malloc_consolidate` | 0.9 | 1.3 | 4.8 (rev_2D) | 3 |
| `Vec3D::length (inline)` | 0.9 | 0.9 | 7.6 (sphere) | 1 |
| `SimulVolume::add_member (inline)` | 0.8 | 3.7 | 7.2 (scale_10k) | 1 |
| `__gnu_cxx::__normal_iterator::__normal_iterator (inline)` | 0.7 | 0.7 | 1.7 (hist_mol_states) | 0 |
| `pirr_pfree_ratio_psF` | 0.7 | 3.2 | 3.1 (trimer) | 0 |
| `__GI___exp` | 0.6 | 2.1 | 2.8 (hist_mol_states) | 0 |
| `__erf` | 0.6 | 1.1 | 2.7 (homoTrimer) | 0 |
| `matrix_rotate (inline)` | 0.6 | 0.6 | 2.7 (homoTrimer) | 0 |
| `libc_feholdsetround_sse_ctx` | 0.6 | 0.6 | 1.4 (cluster_gagsphere) | 0 |
| `Vec3D::operator/= (inline)` | 0.6 | 0.6 | 7.1 (sphere) | 1 |
| `nlohmann::json_abi_v3_12_0::basic_json::assert_invariant` | 0.6 | 0.6 | 4.1 (rev_2D) | 1 |
| `nlohmann::json_abi_v3_12_0::basic_json::json_value::destroy` | 0.5 | 4.3 | 3.7 (rev_2D) | 1 |
| `__mpn_divrem` | 0.5 | 0.5 | 3.8 (gagsphere) | 1 |
| `_int_malloc` | 0.5 | 0.6 | 2.6 (cluster_mem_loc) | 0 |
| `check_dissociation_implicitlipid` | 0.5 | 1.2 | 6.7 (implicit_lipid) | 1 |
| `__memcpy_evex_unaligned_erms` | 0.5 | 0.5 | 2.5 (mem_localization) | 0 |
| `reflect_traj_complex_rad_rot` | 0.5 | 4.3 | 1.4 (implicit_lipid) | 0 |
| `do_cos` | 0.5 | 0.5 | 1.3 (hexamer) | 0 |
| `add_3D_rotational_diffusion (inline)` | 0.5 | 0.6 | 1.8 (cluster_gagsphere) | 0 |

## Self % per case, largest functions

| case | `measure_separations_to_ident` | `std::vector::operator[] (inl` | `check_bimolecular_reactions` | `Complex::update_properties` | `create_complex_propagation_v` | `gsl_ran_gaussian_ziggurat` | `gsl_ran_wishart` | `std::vector::size (inline)` | `__ieee754_exp_fma` | `main` | `std::vector::push_back (inli` | `do_sin` | `operator- (inline)` | `__sincos_fma` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| clathrin | 8.2 | 6.0 | 3.1 | 5.6 | 5.2 | 7.4 | 4.1 | 2.3 | 1.5 | 2.0 | 1.4 | 2.8 | 2.7 | 2.8 |
| closed_homoTrimer | 8.4 | 7.7 | 9.6 | 0.9 | 2.6 | 1.5 | 2.8 | 3.3 | 7.4 | 1.1 | 0.2 | 3.3 | 2.4 | 2.2 |
| cluster_gagsphere | 4.3 | 3.6 | 6.8 | 0.7 | 0.4 | 0.0 | 1.4 | 1.4 | 1.8 | 0.4 | 0.7 | 2.9 | 0.4 | 1.1 |
| cluster_mem_loc | 6.1 | 1.8 | 1.5 | 7.1 | 2.8 | 3.3 | 1.1 | 0.7 | 0.2 | 3.7 | 0.9 | 0.0 | 1.9 | 2.0 |
| compartment | 9.9 | 2.2 | 1.0 | 3.5 | 4.5 | 7.4 | 5.5 | 1.7 | 1.2 | 1.7 | 0.5 | 2.0 | 3.0 | 2.2 |
| enzyme | 12.6 | 8.0 | 8.4 | 7.2 | 4.9 | 3.6 | 2.2 | 3.3 | 0.4 | 4.1 | 2.2 | 0.5 | 3.0 | 2.8 |
| gagsphere | 10.7 | 4.5 | 3.1 | 1.7 | 2.5 | 2.8 | 2.6 | 2.0 | 1.8 | 1.4 | 1.4 | 1.5 | 2.4 | 1.4 |
| hetTrimer | 15.8 | 8.6 | 18.4 | 1.7 | 1.5 | 2.0 | 1.7 | 6.3 | 2.6 | 1.2 | 0.3 | 3.5 | 3.8 | 1.9 |
| hexamer | 13.4 | 5.0 | 4.5 | 4.7 | 5.7 | 5.3 | 4.0 | 2.5 | 2.9 | 1.8 | 1.3 | 6.4 | 1.4 | 3.3 |
| hist_mol_states | 11.4 | 6.8 | 7.4 | 5.7 | 5.7 | 8.0 | 7.4 | 3.4 | 6.8 | 4.0 | 0.6 | 0.0 | 1.1 | 0.0 |
| homoTrimer | 8.8 | 7.1 | 9.8 | 1.0 | 2.3 | 1.9 | 2.8 | 3.2 | 7.7 | 1.1 | 0.2 | 4.3 | 1.3 | 1.5 |
| implicit_lipid | 7.1 | 5.9 | 0.0 | 9.4 | 5.5 | 7.2 | 4.6 | 2.8 | 2.8 | 5.5 | 1.3 | 1.6 | 0.8 | 0.7 |
| mem_localization | 8.0 | 2.5 | 1.2 | 6.6 | 2.9 | 2.6 | 1.6 | 1.0 | 0.1 | 3.5 | 0.9 | 0.0 | 1.5 | 1.5 |
| michaelis_menten | 12.2 | 6.0 | 7.4 | 6.4 | 3.7 | 4.4 | 5.6 | 4.6 | 4.2 | 2.5 | 0.5 | 0.0 | 1.3 | 1.4 |
| rev_2D | 3.0 | 0.7 | 0.0 | 0.4 | 1.1 | 0.4 | 0.4 | 0.4 | 0.0 | 0.4 | 0.0 | 4.5 | 0.0 | 1.9 |
| rev_3D | 24.2 | 6.7 | 0.2 | 6.6 | 7.6 | 4.1 | 2.1 | 3.2 | 0.6 | 3.4 | 9.2 | 0.0 | 0.2 | 0.0 |
| rev_3Dto2D | 15.7 | 5.7 | 0.9 | 7.6 | 4.3 | 3.4 | 1.9 | 1.5 | 0.2 | 3.6 | 4.3 | 0.0 | 1.6 | 2.0 |
| scale_10k | 25.4 | 6.9 | 0.8 | 2.9 | 6.8 | 1.3 | 0.9 | 2.1 | 0.2 | 4.2 | 14.4 | 0.0 | 0.2 | 0.0 |
| scale_2k | 24.3 | 6.5 | 0.3 | 5.6 | 6.8 | 4.3 | 2.2 | 3.4 | 0.4 | 3.2 | 9.9 | 0.0 | 0.3 | 0.2 |
| sphere | 2.7 | 2.4 | 0.0 | 4.8 | 2.5 | 4.9 | 3.9 | 1.1 | 1.5 | 2.7 | 1.2 | 2.0 | 1.7 | 1.1 |
| trimer | 7.5 | 8.1 | 16.3 | 0.5 | 0.7 | 0.7 | 1.7 | 5.7 | 7.5 | 0.4 | 0.1 | 4.1 | 3.2 | 1.9 |
| unimol_state_change | 16.9 | 9.0 | 10.1 | 2.8 | 2.3 | 1.5 | 3.7 | 5.3 | 5.5 | 2.1 | 0.4 | 0.0 | 1.4 | 1.6 |

## Inclusive % per case, selected frames

| case | `measure_separations_to_ident` | `check_bimolecular_reactions` | `Complex::propagate` | `create_complex_propagation_v` | `SimulVolume::update_memberMo` | `Complex::update_properties` | `gsl_ran_gaussian_ziggurat` | `QuatRotation::rotate (inline` | `determine_3D_bimolecular_rea` | `reflect_traj_complex_rad_rot` | `SimulVolume::add_member (inl` | `ShellIndex::rebin` | `remove_empty_slots` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| clathrin | 26.9 | 15.1 | 38.9 | 17.3 | 4.6 | 7.7 | 11.8 | 15.6 | 7.2 | 4.3 | 2.3 | 0.0 | 0.0 |
| closed_homoTrimer | 59.5 | 49.5 | 7.0 | 5.3 | 1.0 | 1.6 | 2.6 | 1.2 | 31.7 | 0.4 | 0.4 | 0.0 | 0.0 |
| cluster_gagsphere | 71.1 | 63.2 | 4.3 | 2.9 | 0.7 | 1.1 | 1.4 | 2.9 | 8.6 | 1.4 | 0.7 | 0.0 | 0.0 |
| cluster_mem_loc | 9.9 | 2.4 | 13.0 | 7.5 | 4.3 | 10.4 | 4.3 | 0.0 | 0.3 | 12.0 | 2.2 | 0.0 | 0.0 |
| compartment | 26.6 | 5.0 | 33.7 | 20.1 | 6.2 | 5.0 | 14.1 | 13.4 | 3.2 | 5.2 | 1.2 | 0.0 | 0.0 |
| enzyme | 32.7 | 16.4 | 19.4 | 11.2 | 7.2 | 11.9 | 6.1 | 0.9 | 1.6 | 14.6 | 3.5 | 0.0 | 0.0 |
| gagsphere | 28.8 | 16.4 | 13.3 | 8.0 | 4.5 | 2.5 | 5.1 | 4.2 | 9.1 | 1.1 | 2.5 | 0.0 | 0.0 |
| hetTrimer | 64.4 | 43.9 | 10.7 | 4.7 | 1.0 | 2.7 | 3.1 | 2.2 | 12.6 | 0.7 | 0.5 | 0.0 | 0.0 |
| hexamer | 35.3 | 18.7 | 27.0 | 15.6 | 4.8 | 6.8 | 9.1 | 4.6 | 10.2 | 2.6 | 2.4 | 0.0 | 0.0 |
| hist_mol_states | 31.2 | 13.6 | 14.8 | 18.2 | 5.1 | 11.9 | 12.5 | 0.0 | 3.4 | 5.1 | 1.7 | 0.0 | 0.0 |
| homoTrimer | 58.3 | 47.7 | 7.5 | 5.7 | 1.2 | 1.5 | 3.2 | 1.5 | 31.2 | 0.8 | 0.6 | 0.0 | 0.0 |
| implicit_lipid | 19.9 | 0.0 | 23.4 | 18.0 | 5.5 | 13.8 | 11.3 | 1.8 | 0.0 | 4.9 | 1.9 | 0.0 | 0.0 |
| mem_localization | 12.2 | 3.2 | 12.9 | 7.8 | 3.4 | 10.5 | 4.4 | 0.0 | 0.9 | 11.5 | 1.6 | 0.0 | 0.0 |
| michaelis_menten | 41.4 | 24.3 | 12.1 | 12.0 | 2.6 | 9.6 | 7.8 | 0.0 | 13.5 | 4.7 | 1.0 | 0.0 | 0.0 |
| rev_2D | 33.8 | 30.1 | 1.5 | 1.9 | 0.4 | 0.7 | 0.7 | 0.0 | 0.0 | 1.9 | 0.0 | 0.0 | 0.0 |
| rev_3D | 33.2 | 1.4 | 16.2 | 14.8 | 19.5 | 9.8 | 6.8 | 0.0 | 0.8 | 3.2 | 13.2 | 0.0 | 0.0 |
| rev_3Dto2D | 22.5 | 1.9 | 14.0 | 10.0 | 10.6 | 11.0 | 5.5 | 0.0 | 0.4 | 12.0 | 6.0 | 0.0 | 0.0 |
| scale_10k | 33.1 | 2.1 | 9.4 | 9.5 | 28.2 | 4.1 | 2.5 | 0.0 | 0.6 | 1.5 | 22.4 | 0.0 | 0.0 |
| scale_2k | 34.4 | 1.6 | 15.4 | 14.1 | 20.1 | 8.9 | 6.7 | 0.0 | 0.8 | 3.0 | 13.8 | 0.0 | 0.0 |
| sphere | 7.3 | 0.0 | 37.9 | 27.1 | 3.6 | 7.1 | 8.5 | 0.4 | 0.0 | 1.6 | 1.7 | 0.0 | 0.0 |
| trimer | 76.1 | 65.7 | 3.0 | 1.9 | 0.3 | 0.7 | 0.9 | 0.5 | 32.0 | 0.6 | 0.2 | 0.0 | 0.0 |
| unimol_state_change | 57.4 | 35.5 | 5.8 | 5.1 | 1.9 | 4.0 | 2.6 | 0.0 | 17.6 | 2.5 | 0.9 | 0.1 | 0.0 |
