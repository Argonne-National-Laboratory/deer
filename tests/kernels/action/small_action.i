# Simple 3D test

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 4
[]

[NEMLMechanics]
  displacements = "disp_x disp_y disp_z"
  kinematics = small
  add_all_output = true
[]

[BCs]
  [./leftx]
    type = PresetBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./lefty]
    type = PresetBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./leftz]
    type = PresetBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./pull_z]
    type = PresetBC
    boundary = front
    variable = disp_z
    value = 0.1
  [../]
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "elastic_model"
    large_kinematics = false
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  
  dt = 1
  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
