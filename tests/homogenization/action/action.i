[GlobalParams]
  displacements = "disp_x disp_y disp_z"
  large_kinematics = true
[]

[Mesh]
  [./base]
    type = FileMeshGenerator
    file = '3d.exo'
  [../]

  [./sidesets]
    type = SideSetsFromNormalsGenerator
    input = base
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0
                0 0 -1
                0 0 1'
    fixed_normal = true
    new_boundary = 'left right bottom top back front'
  [../]
[]

[NEMLMechanics]
  kinematics = large
  add_all_output = true
  formulation = total
  homogenize = true
  constraint_types = 'stress strain strain stress stress strain stress stress stress'
  targets = 'stress11 zero zero zero zero zero zero zero zero'
[]

[Functions]
  [./stress11]
    type = ParsedFunction
    value = '8.0e2*t'
  [../]
  [./zero]
    type = ConstantFunction
    value = 0
  [../]
[]

[BCs]
  [./Periodic]
    [./x]
      variable = disp_x
      auto_direction = 'x y z'
    [../]
    [./y]
      variable = disp_y
      auto_direction = 'x y z'
    [../]
    [./z]
      variable = disp_z
      auto_direction = 'x y z'
    [../]
  [../]

  [./fix1_x]
    type = DirichletBC
    boundary = "fix_all"
    variable = disp_x
    value = 0
  [../]
  [./fix1_y]
    type = DirichletBC
    boundary = "fix_all"
    variable = disp_y
    value = 0
  [../]
  [./fix1_z]
    type = DirichletBC
    boundary = "fix_all"
    variable = disp_z
    value = 0
  [../]

  [./fix2_x]
    type = DirichletBC
    boundary = "fix_xy"
    variable = disp_x
    value = 0
  [../]
  [./fix2_y]
    type = DirichletBC
    boundary = "fix_xy"
    variable = disp_y
    value = 0
  [../]

  [./fix3_z]
    type = DirichletBC
    boundary = "fix_z"
    variable = disp_z
    value = 0
  [../]
[]

[Materials]
  [./stress1]
    type = ComputeNEMLStressUpdate
    database = "materials.xml"
    model = "mat1"
    block = '1'
  [../]
  [./stress2]
    type = ComputeNEMLStressUpdate
    database = "materials.xml"
    model = "mat2"
    block = '2'
  [../]
  [./stress3]
    type = ComputeNEMLStressUpdate
    database = "materials.xml"
    model = "mat3"
    block = '3'
  [../]
  [./stress4]
    type = ComputeNEMLStressUpdate
    database = "materials.xml"
    model = "mat4"
    block = '4'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./s11]
    type = ElementAverageValue
    variable = stress_x-x
    execute_on = 'initial timestep_end'
  [../]
  [./s21]
    type = ElementAverageValue
    variable = stress_y-z
    execute_on = 'initial timestep_end'
  [../]
  [./s31]
    type = ElementAverageValue
    variable = stress_x-z
    execute_on = 'initial timestep_end'
  [../]
  [./s12]
    type = ElementAverageValue
    variable = stress_x-y
    execute_on = 'initial timestep_end'
  [../]
  [./s22]
    type = ElementAverageValue
    variable = stress_y-y
    execute_on = 'initial timestep_end'
  [../]
  [./s32]
    type = ElementAverageValue
    variable = stress_y-z
    execute_on = 'initial timestep_end'
  [../]
  [./s13]
    type = ElementAverageValue
    variable = stress_x-z
    execute_on = 'initial timestep_end'
  [../]
  [./s23]
    type = ElementAverageValue
    variable = stress_y-z
    execute_on = 'initial timestep_end'
  [../]
  [./s33]
    type = ElementAverageValue
    variable = stress_z-z
    execute_on = 'initial timestep_end'
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
  nl_max_its = 10
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 0.2
  dtmin = 0.2
  end_time = 1.0
[]

[Outputs]
  exodus = false
  csv = true
[]
