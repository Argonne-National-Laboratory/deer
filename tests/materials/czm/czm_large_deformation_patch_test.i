# Patch test for cohesive zone modeling to check convergence
[Mesh]
  [./msh]
  type = FileMeshGenerator
  file = patch_mesh.e
  []
  [./transform]
  type = TransformGenerator
  input = msh
  transform = TRANSLATE
  vector_value = '-0.5 -0.5 -0.5'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = transform
  []
  [./add_surfaces]
    type = SideSetsFromNormalsGenerator
    input = split
    normals = '0  0  1
               0  1  0
               1  0  0
               0  0 -1
               0 -1  0
              -1  0  0'
    fixed_normal = true
    new_boundary = 'front top right back bottom left'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
  formulation = total
[]

[Functions]
  [./angles]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0 0.15707963'
  [../]

  [./stretch]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0.01 0.01'
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

  [./dt_fun]
    type = PiecewiseConstant
    x = '0 2'
    y = '0.25 0.25'
  []
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


[CohesiveZoneDeer]
   boundary = 'interface'
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "mat.xml"
    model = "elastic_model"
    large_kinematics = true
  [../]
  [./czm_mat]
    boundary = 'interface'
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

  solve_type = 'NEWTON'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  start_time = 0.0
  end_time = 2
  [./TimeStepper]
    type = FunctionDT
    function = dt_fun
  [../]
[]

[Postprocessors]
  [./nonlin]
    type = NumNonlinearIterations
  [../]
[]


[Outputs]
  exodus = false
  csv = true
[]
