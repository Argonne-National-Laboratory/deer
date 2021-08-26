# Simple 3D test

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  large_kinematics = true
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
  nz = 1
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.5
  zmax = 0.5
  []
[]

[NEMLMechanics]
  kinematics = large
  formulation = total
  add_all_output = true
  output_internal_strain = true
[]


[Functions]
  [./angles]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0 1.5707963'
  [../]

  [./stretch]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0.1 0.1'
  [../]

  [./move_y]
    type = ParsedFunction
    value = 'y*cos(theta) - z * (1 + a)*sin(theta) - y'
    vars = 'a theta'
    vals = 'stretch angles'
  [../]

  [./move_z]
    type = ParsedFunction
    value = 'y*sin(theta) + z*(1+a)*cos(theta) - z'
    vars = 'a theta'
    vals = 'stretch angles'
  [../]
[]

[BCs]
  [./fix]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = left
    variable = disp_x
  [../]

  [./front_y]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_y
    function = move_y
    preset = true
  [../]

  [./back_y]
    type = FunctionDirichletBC
    boundary = back
    variable = disp_y
    function = move_y
    preset = true
  [../]

  [./front_z]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = move_z
    preset = true
  [../]

  [./back_z]
    type = FunctionDirichletBC
    boundary = back
    variable = disp_z
    function = move_z
    preset = true
  [../]

[]


[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "elastic_model"
  [../]
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
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = 2.0
  [./TimeStepper]
    type = FunctionDT
    time_t = '0 1 2'
    time_dt = '0.1 0.001 0.001'
    interpolate = False
  [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
