# 2D test with stress controls

[GlobalParams]
  displacements = 'disp_x disp_y'
  constraint_types = 'stress stress stress'
  ndim = 2
  large_kinematics = false
  macro_gradient = hvar
[]

[Mesh]
  [./base]
    type = FileMeshGenerator
    file = '2d.exo'
  [../]

  [./sidesets]
    type = SideSetsFromNormalsGenerator
    input = base
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0'
    fixed_normal = true
    new_boundary = 'left right bottom top'
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
 [./hvar]
    family = SCALAR
    order = THIRD
  [../]
[]

[AuxVariables]
  [./sxx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./syy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./sxy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./exx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./eyy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./exy]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./sxx]
    type = RankTwoAux
    variable = sxx
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
  [./syy]
    type = RankTwoAux
    variable = syy
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]
  [./sxy]
    type = RankTwoAux
    variable = sxy
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
  [../]
  [./exx]
    type = RankTwoAux
    variable = exx
    rank_two_tensor = mechanical_strain_unrotated
    index_i = 0
    index_j = 0
  [../]
  [./eyy]
    type = RankTwoAux
    variable = eyy
    rank_two_tensor = mechanical_strain_unrotated
    index_i = 1
    index_j = 1
  [../]
  [./exy]
    type = RankTwoAux
    variable = exy
    rank_two_tensor = mechanical_strain_unrotated
    index_i = 0
    index_j = 1
  [../]
[]

[UserObjects]
  [./integrator]
    type = HomogenizationConstraintIntegral
    targets = 'stress11 stress22 stress12'
    execute_on = 'initial linear'
  [../]
[]

[Kernels]
  [./sdx]
    type = TotalStressDivergenceNEML
    variable = disp_x
    component = 0
  [../]
  [./sdy]
    type = TotalStressDivergenceNEML
    variable = disp_y
    component = 1
  [../]
[]

[ScalarKernels]
  [./enforce]
    type = HomogenizationConstraintScalarKernel
    variable = hvar
    integrator = integrator
  [../]
[]

[Functions]
  [./stress11]
    type = ParsedFunction
    value = '400*t'
  [../]
  [./stress22]
    type = ParsedFunction
    value = '-200*t'
  [../]
  [./stress12]
    type = ParsedFunction
    value = '100*t'
  [../]
[]

[BCs]
  [./Periodic]
    [./x]
      variable = disp_x
      auto_direction = 'x y'
    [../]
    [./y]
      variable = disp_y
      auto_direction = 'x y'
    [../]
  [../]

  [./fix1_x]
    type = DirichletBC
    boundary = "fix1"
    variable = disp_x
    value = 0
  [../]
  [./fix1_y]
    type = DirichletBC
    boundary = "fix1"
    variable = disp_y
    value = 0
  [../]

  [./fix2_y]
    type = DirichletBC
    boundary = "fix2"
    variable = disp_y
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
  [./strain]
    type = ComputeNEMLStrain
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./sxx]
    type = ElementAverageValue
    variable = sxx
    execute_on = 'initial timestep_end'
  [../]
  [./syy]
    type = ElementAverageValue
    variable = syy
    execute_on = 'initial timestep_end'
  [../]
  [./sxy]
    type = ElementAverageValue
    variable = sxy
    execute_on = 'initial timestep_end'
  [../]
  [./exx]
    type = ElementAverageValue
    variable = exx
    execute_on = 'initial timestep_end'
  [../]
  [./eyy]
    type = ElementAverageValue
    variable = eyy
    execute_on = 'initial timestep_end'
  [../]
  [./exy]
    type = ElementAverageValue
    variable = exy
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
