[Mesh]
  type = FileMesh
  file = 'n10.exo'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
      [./disp_x]
            order = second
      [../]
      [./disp_y]
            order = second
      [../]
      [./disp_z]
            order = second
      [../]
[]

[Functions]
  [./pfn]
    type = PiecewiseLinear
    x = '0    100'
    y = '0.00 0.01'
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

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_10_rand.tex
  [../]
[]
[BCs]
  [./left]
     type = DirichletBC
     preset = true
     variable = disp_x
     boundary = left
     value = 0.0
  [../]

  [./bottom]
    type = DirichletBC
     preset = true
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./back]
    type = DirichletBC
     preset = true
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./front]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = pfn
    preset = true
  [../]
[]

[Materials]
  [./strain]
    type = ComputeNEMLStrain
    large_kinematics = true
  [../]
  [./stress]
    type = ComputeNEMLCPOutput
    database = "test.xml"
    model = "grain_1"
    large_kinematics = true
    euler_angle_provider = euler_angle_file
    # grain_id = 0
    # block = 1
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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  nl_max_its = 20

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8

  end_time = 10.0
  dtmin = 0.1
  dt = 1.0
[]

[Outputs]
  exodus = true
[]
