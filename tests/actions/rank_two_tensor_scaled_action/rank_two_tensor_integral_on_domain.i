# Simple 3D test

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
      [./disp_x]
      [../]
      [./disp_y]
      [../]
      [./disp_z]
      [../]
[]

[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '1 1 1'
  []
  [./interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = new_block
    master_block = 0
    paired_block = 1
    new_boundary = 'interface'
  []
[]

[AuxVariables]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./zstress]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 500'
  [../]
  [./constant]
    type = ConstantFunction
    value = 1.0
  [../]
[]

[BCs]
  [./leftx]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./boty]
    type = DirichletBC
    preset = true
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./backz]
    type = DirichletBC
    preset = true
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./pull_z]
    type = FunctionNeumannBC
    boundary = front
    variable = disp_z
    function = zstress
    use_displaced_mesh = true
  [../]
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "elastic_model"
    large_kinematics = true
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
[]

[Postprocessors]
  [./area]
    type = AreaPostprocessor
    boundary = interface
    use_displaced_mesh = true
  [../]
  [./area0]
    type = AreaPostprocessor
    boundary = interface
    use_displaced_mesh = false
  [../]
  [./volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
  [../]
  [./volume0]
    type = VolumePostprocessor
    use_displaced_mesh = false
  [../]
[]

[RankTwoTensorIntegralOnDomain]
  [integral1]
    rank_two_tensor = 'stress'
    use_displaced_mesh = true
    base_out_names = 'stress_int1'
  []
  [integral2]
    rank_two_tensor = 'stress'
    use_displaced_mesh = true
    base_out_names = 'stress_int2'
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
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Outputs]
  exodus = false
  csv = true
[]
