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
    E = 1e2
    G = 1e2
    interface_thickness = 1
  [../]
  [./czm_voluemtric_strain]
    type = CZMVolumetricStrain
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    use_displaced_mesh = false
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
  dtmin = 1
  dtmax = 1
  end_time = 1
[]

[BCs]
  [./x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./x_top]
    type = DirichletBC
    boundary = front
    variable = disp_x
    value = 0.0
  [../]
  [./y_top]
    type = DirichletBC
    boundary = front
    variable = disp_y
    value = 0.0
  [../]
  [./z_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = disp_fun
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 0.1'
  [../]
[]

[Postprocessors]
  [./A0]
    type = AreaPostprocessor
    boundary = interface
    use_displaced_mesh = false
  []
  [./V0]
    type = VolumePostprocessor
    use_displaced_mesh = false
    execute_on = 'INITIAL'
  []
  [./czm_strain_rate_total_xx]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = czm_total_strain_rate
    index_i = 0
    index_j = 0
    boundary = interface
  [../]
  [./czm_strain_rate_total_yy]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = czm_total_strain_rate
    index_i = 1
    index_j = 1
    boundary = interface
  [../]
  [./czm_strain_rate_total_zz]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = czm_total_strain_rate
    index_i = 2
    index_j = 2
    boundary = interface
    execute_on = 'TIMESTEP_END'
    scaling_factor_PP = V0
  [../]
  [./czm_normal_strain_rate_zz]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = czm_normal_strain_rate
    index_i = 2
    index_j = 2
    boundary = interface
    execute_on = 'TIMESTEP_END'
    scaling_factor_PP = V0
  [../]
  [./czm_sliding_strain_rate_zz]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = czm_sliding_strain_rate
    index_i = 2
    index_j = 2
    boundary = interface
    execute_on = 'TIMESTEP_END'
    scaling_factor_PP = V0
  [../]
  [./czm_sliding_strain_rate_zz2]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = czm_sliding_strain_rate
    index_i = 2
    index_j = 2
    boundary = interface
    execute_on = 'TIMESTEP_END'
  [../]
[]

[Outputs]
  exodus = true
[]
