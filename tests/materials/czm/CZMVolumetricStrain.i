[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
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
[]


[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = small
  add_all_output = true
  add_displacements = true
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm]
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "../../neml_test_material.xml"
    model = "elastic_model"
    large_kinematics = false
  [../]
  [./czm]
    type = PureElasticCZM
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    E = 1e0
    G = 1e0
    interface_thickness = 1
  [../]
  [./czm_voluemtric_strain]
    type = CZMVolumetricStrain
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    large_kinematics = true
  []
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
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

[Controls]
  [./c1]
    type = TimePeriod
    enable_objects = 'BCs::fix_x BCs::fix_y BCs::fix_z BCs::move_top'
    disable_objects = 'BCs::rotate_x BCs::rotate_y BCs::rotate_z'
    start_time = '0'
    end_time = '1'
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 0.1'
  [../]
  [./dt_fun]
    type = PiecewiseConstant
    x = '0 0.99 2'
    y = '0.01 0.001 0.001'
  []
[]

[Postprocessors]
  [./V0]
    type = VolumePostprocessor
    use_displaced_mesh = false
    execute_on = 'INITIAL'
  []
[]

[RankTwoTensorIntegralOnDomain]
 [czm]
   boundary = interface
   rank_two_tensor = 'czm_total_strain_rate czm_normal_strain_rate czm_sliding_strain_rate'
   use_displaced_mesh = false
   base_out_names = 'czm_strain_total_rate czm_strain_normal_rate czm_strain_sliding_rate'
   scaling_factor_PP =V0
 []
[]

[RankTwoTensorPostprocessorTimeIntegral]
  pp_base_names = 'czm_strain_total_rate czm_strain_normal_rate czm_strain_sliding_rate'
  base_out_names = 'czm_strain_total czm_strain_normal czm_strain_sliding'
[]

[Outputs]
  exodus = true
[]
