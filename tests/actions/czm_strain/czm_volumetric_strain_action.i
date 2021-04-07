[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmax = 3
  ymax = 2
  zmax = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 1'
    top_right = '3 2 2'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[AuxVariables]
  [./T_N]
    family = MONOMIAL
    order = CONSTANT
  []
  [./T_S1]
    family = MONOMIAL
    order = CONSTANT
  []
  [./T_S2]
    family = MONOMIAL
    order = CONSTANT
  []
  [./jump_z]
    family = MONOMIAL
    order = CONSTANT
  []
[]




[AuxKernels]
  [./aux_TN]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = T_N
  [../]
  [./aux_TS1]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = T_S1
  [../]
  [./aux_TS2]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = T_S2
  [../]
  [./aux_jump_z]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = displacement_jump_global
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = jump_z
  [../]
[]


[NEMLMechanics]
  kinematics = large
  add_all_output = true
  add_displacements = true
[]

[CohesiveZoneDeer]
   boundary = 'interface'
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "../../neml_test_material.xml"
    model = "powerlaw"
    large_kinematics = true
  [../]
  [./czm]
    type = PureElasticCZMIncremental
    boundary = 'interface'
    E = 1e4
    G = 1e4
  [../]
  [./tensor_rate_mech]
  type = TensorRateMaterial
  rank_two_tensor = mechanical_strain_rotated
  convert_log_strain_to_eng_strain = true
  []
  [./tensor_rate_el]
    type = TensorRateMaterial
    rank_two_tensor = elastic_strain
    convert_log_strain_to_eng_strain = true
  []
  [./tensor_rate_inel]
    type = TensorRateMaterial
    rank_two_tensor = inelastic_strain_rotated
    convert_log_strain_to_eng_strain = true
  []
  [czm_strain]
    type = CZMVolumetricStrain
    boundary = interface
    use_large_kineamtics = true
  []
[]


[Postprocessors]
  [./A0]
    type = InterfaceAreaPostprocessor
    boundary = 'interface'
    use_displaced_mesh = false
    execute_on = 'INITIAL'
  []
  [./V0]
    type = VolumePostprocessor
    use_displaced_mesh = false
    block = '0 1'
    execute_on = 'INITIAL'
  []
  [CZMStrainScaling]
    type = CZMStrainScalingPostprocessor
    A0 = A0
    V0 = V0
  []
[]

[RankTwoTensorIntegralOnDomain]
  [czm_strain]
    rank_two_tensor = 'czm_total_strain czm_normal_strain czm_sliding_strain'
    base_out_names = 'czm_total_strain czm_normal_strain czm_sliding_strain'
    scaling_factor_PP = CZMStrainScaling
    normalize_integral_by_area = true
    boundary = interface
    czm = true
  []
[]

[RankTwoTensorPostprocessorTimeDerivative]
  [czm_strain_rate]
    pp_base_names = 'czm_total_strain czm_normal_strain czm_sliding_strain'
    base_out_names = 'czm_total_strain_rate czm_normal_strain_rate czm_sliding_strain_rate'
  []
[]


[RankTwoTensorIntegralOnDomain]
  [strain_average_average_on_bulk]
    rank_two_tensor = 'elastic_strain_rate inelastic_strain_rotated_rate mechanical_strain_rotated_rate'
    use_displaced_mesh = true
    base_out_names = 'grain_eel_rate grain_einel_rate grain_etot_rate'
    scaling_factor_PP = V0
    block = '0 1'
  []
[]

[RankTwoTensorPostprocessorTimeIntegral]
  [timeintegral1]
    pp_base_names = 'grain_eel_rate grain_einel_rate grain_etot_rate'
    base_out_names = 'grain_eel grain_einel grain_etot'
  []
[]

[RankTwoTensorPostprocessorSum]
  [RVE_strain]
    pp_base_names_1 = 'grain_etot grain_etot_rate'
    pp_base_names_2 = 'czm_total_strain czm_total_strain_rate'
    base_out_names = 'rve_e rve_e_rate'
  []
[]

[Postprocessors]
  [czm_total_strain_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'czm_total_strain_eq'
  []
  [czm_total_strain_rate_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'czm_total_strain_rate_eq'
  []
  [grain_etot_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'grain_etot'
  []
  [grain_einel_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'grain_einel'
  []
  [grain_eel_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'grain_eel'
  []
  [grain_etot_rate_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'grain_etot_rate'
  []
  [grain_einel_rate_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'grain_einel_rate'
  []
  [grain_eel_rate_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'grain_eel_rate'
  []
  [rve_e_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'rve_e'
  []
  [rve_e_rate_eq]
    type = RankTwoTensorInvariantPostprocessor
    invariant = 'EffectiveStrain'
    rank_two_tensor_base_name = 'rve_e_rate'
  []
[]



[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]


[BCs]
  [./fix_x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./fix_y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./fix_z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./move_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = disp_fun
  [../]

  [./rotate_x]
  type = DisplacementAboutAxis
  boundary = 'front back left right top bottom'
  function = '90.'
  angle_units = degrees
  axis_origin = '0. 0. 0.'
  axis_direction = '0. 1. 0.'
  component = 0
  variable = disp_x
  angular_velocity = true
[../]
[./rotate_y]
  type = DisplacementAboutAxis
  boundary = 'front back left right top bottom'
  function = '90.'
  angle_units = degrees
  axis_origin = '0. 0. 0.'
  axis_direction = '0. 1. 0.'
  component = 1
  variable = disp_y
  angular_velocity = true
[../]
[./rotate_z]
  type = DisplacementAboutAxis
  boundary = 'front back left right top bottom'
  function = '90.'
  angle_units = degrees
  axis_origin = '0. 0. 0.'
  axis_direction = '0. 1. 0.'
  component = 2
  variable = disp_z
  angular_velocity = true
[../]
[]

# Constraint System
[Constraints]
  [./x1]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    secondary = 'right'    # boundary
    penalty = 1e6
  [../]
  [./y1]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = 'top'    # boundary
    penalty = 1e6
  [../]
[]

[Controls]
  [./c1]
    type = TimePeriod
    enable_objects = 'BCs::fix_x BCs::fix_y BCs::fix_z BCs::move_top Constraints::x1 Constraints::y1'
    disable_objects = 'BCs::rotate_x BCs::rotate_y BCs::rotate_z'
    start_time = '0'
    end_time = '1'
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1.'
  [../]
  [./dt_fun]
    type = PiecewiseConstant
    x = '0 0.99 2'
    y = '0.01 0.001 0.001'
  []
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  [./TimeStepper]
    type = FunctionDT
    function = dt_fun
  [../]
  end_time = 2
  # num_steps = 2
[]

[Outputs]
  csv = true
  exodus=true
  sync_times = '0.1 0.5 1'
[]
