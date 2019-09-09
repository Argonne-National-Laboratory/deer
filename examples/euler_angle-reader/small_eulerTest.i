# Simple 3D test

[GlobalParams]
  displacements = "disp_x disp_y disp_z"
[]

[Mesh]
  type = FileMesh
  file = "small.exo"
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]
[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = 2gr_Euler.tex
  [../]
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
    boundary = back
    variable = disp_y
    value = 0.0
  [../]
  [./leftz]
    type = PresetBC
    boundary = bottom
    variable = disp_z
    value = 0.0
  [../]
  [./pull_z]
    type = FunctionPresetBC
    boundary = top
    variable = disp_z
    function = pfn
  [../]
[]

[Functions]
  [./pfn]
    type = PiecewiseLinear
    x = '0    10'
    y = '0.00 0.1'
  [../]
[]

[AuxVariables]
  [./orientation_q1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./orientation_q2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./orientation_q3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./orientation_q4]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [q1]
    type = MaterialStdVectorAux
    property  = orientation_q
    index = 0
    variable  = orientation_q1
  [../]
  [q2]
    type = MaterialStdVectorAux
    property  = orientation_q
    index = 1
    variable  = orientation_q2
  [../]
  [q3]
    type = MaterialStdVectorAux
    property  = orientation_q
    index = 2
    variable  = orientation_q3
  [../]
  [q4]
    type = MaterialStdVectorAux
    property  = orientation_q
    index = 3
    variable  = orientation_q4
  [../]
[]

[Materials]
  [./strain]
    type = ComputeNEMLStrain
    large_kinematics = true
  [../]
  [./stress1]
    type = ComputeNEMLCPOutput
    database = "test.xml"
    model = "grain_1"
    large_kinematics = true
    euler_angle_provider = euler_angle_file
    # block = 1
    # grain_id = 0
  [../]
[]
[Kernels]
  [./sdx]
      type = StressDivergenceNEML
      variable = disp_x
      component = 0
      use_displaced_mesh = true
  [../]
  [./sdy]
      type = StressDivergenceNEML
      variable = disp_y
      component = 1
      use_displaced_mesh = true
  [../]
  [./sdz]
      type = StressDivergenceNEML
      variable = disp_z
      component = 2
      use_displaced_mesh = true
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

  nl_max_its = 7
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0
  dt = 1
  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
