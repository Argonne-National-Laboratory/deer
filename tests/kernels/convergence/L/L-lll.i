[Mesh]
  type = FileMesh
  file = 'L.exo'
[]

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

[Functions]
  [./pfn]
    type = PiecewiseLinear
    x = '0    1    2'
    y = '0.00 0.3 0.5'
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

[BCs]
  [./left]
     type = DirichletBC
     preset = true
     variable = disp_x
     boundary = fix
     value = 0.0
  [../]

  [./bottom]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = fix
    value = 0.0
  [../]

  [./back]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = fix
    value = 0.0
  [../]

  [./front]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = pull
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

[Executioner]
  type = Transient

  solve_type = 'newton'

  petsc_options_iname = -pc_type
  petsc_options_value = lu

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8

  end_time = 1.0
  dtmin = 0.5
  dt = 0.5
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
