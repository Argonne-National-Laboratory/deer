[Mesh]
  type = FileMesh
  file = 'cube.exo'
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
    y = '0.00 250 500'
  [../]
[]

[Kernels]
  [./sdx]
      type = StressDivergenceNEML
      variable = disp_x
      component = 0
      use_displaced_mesh = false
  [../]
  [./sdy]
      type = StressDivergenceNEML
      variable = disp_y
      component = 1
      use_displaced_mesh = false
  [../]
  [./sdz]
      type = StressDivergenceNEML
      variable = disp_z
      component = 2
      use_displaced_mesh = false
  [../]
[]

[AuxVariables]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
    execute_on = timestep_end
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    execute_on = timestep_end
  [../]

  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./strain_xz]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_xz
    index_i = 0
    index_j = 2
    execute_on = timestep_end
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_yz
    index_i = 1
    index_j = 2
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./left]
     type = PresetBC
     variable = disp_x
     boundary = left
     value = 0.0
  [../]

  [./bottom]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./back]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./front]
    type = FunctionNeumannBC
    variable = disp_z
    boundary = front
    function = pfn
  [../]
[]

[Materials]
  [./strain]
    type = ComputeNEMLStrain
    use_displaced_mesh = false
  [../]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "elastic_model"
    use_displaced_mesh = false
  [../]
[]

[Postprocessors]
  [./exx]
    type = ElementAverageValue
    variable = 'strain_xx'
    execute_on = nonlinear
  [../]

  [./eyy]
    type = ElementAverageValue
    variable = 'strain_yy'
    execute_on = nonlinear
  [../]

  [./ezz]
    type = ElementAverageValue
    variable = 'strain_zz'
    execute_on = nonlinear
  [../]

  [./exy]
    type = ElementAverageValue
    variable = 'strain_xy'
    execute_on = nonlinear
  [../]

  [./exz]
    type = ElementAverageValue
    variable = 'strain_xz'
    execute_on = nonlinear
  [../]

  [./eyz]
    type = ElementAverageValue
    variable = 'strain_yz'
    execute_on = nonlinear
  [../]

  [./sxx]
    type = ElementAverageValue
    variable = 'stress_xx'
    execute_on = nonlinear
  [../]

  [./syy]
    type = ElementAverageValue
    variable = 'stress_yy'
    execute_on = nonlinear
  [../]

  [./szz]
    type = ElementAverageValue
    variable = 'stress_zz'
    execute_on = nonlinear
  [../]

  [./sxy]
    type = ElementAverageValue
    variable = 'stress_xy'
    execute_on = nonlinear
  [../]

  [./sxz]
    type = ElementAverageValue
    variable = 'stress_xz'
    execute_on = nonlinear
  [../]

  [./syz]
    type = ElementAverageValue
    variable = 'stress_yz'
    execute_on = nonlinear
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
  dt = 1

  solve_type = 'newton'

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_max_its = 10

  end_time = 2
  dtmin = 1.0
[]

[Outputs]
  exodus = false
  [./list]
      type = CSV
      execute_on = 'nonlinear'
  [../]
[]
