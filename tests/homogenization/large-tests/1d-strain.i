# 2D test with just strain control

[GlobalParams]
  displacements = 'disp_x'
  constraint_types = 'strain'
  ndim = 1
  large_kinematics = true
  macro_gradient = hvar
[]

[Mesh]
  type = FileMesh
  file = '1d.exo'
[]

[MeshModifiers]
  [./left]
    type = SideSetsFromPoints
    new_boundary = 'left'
    points = '-1 0 0'
  [../]

  [./right]
    type = SideSetsFromPoints
    new_boundary = 'right'
    points = '7 0 0'
  [../]
[]

[Variables]
  [./disp_x]
  [../]
 [./hvar]
    family = SCALAR
    order = FIRST
  [../]
[]

[AuxVariables]
  [./s11]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s21]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s31]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s12]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s22]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s32]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s13]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s23]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./s33]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./F11]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F21]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F31]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F12]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F22]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F32]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F13]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F23]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./F33]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./s11]
    type = RankTwoAux
    variable = s11
    rank_two_tensor = PK1
    index_i = 0
    index_j = 0
  [../]
  [./s21]
    type = RankTwoAux
    variable = s21
    rank_two_tensor = PK1
    index_i = 1
    index_j = 0
  [../]
  [./s31]
    type = RankTwoAux
    variable = s31
    rank_two_tensor = PK1
    index_i = 2
    index_j = 0
  [../]
  [./s12]
    type = RankTwoAux
    variable = s12
    rank_two_tensor = PK1
    index_i = 0
    index_j = 1
  [../]
  [./s22]
    type = RankTwoAux
    variable = s22
    rank_two_tensor = PK1
    index_i = 1
    index_j = 1
  [../]
  [./s32]
    type = RankTwoAux
    variable = s32
    rank_two_tensor = PK1
    index_i = 2
    index_j = 1
  [../]
  [./s13]
    type = RankTwoAux
    variable = s13
    rank_two_tensor = PK1
    index_i = 0
    index_j = 2
  [../]
  [./s23]
    type = RankTwoAux
    variable = s23
    rank_two_tensor = PK1
    index_i = 1
    index_j = 2
  [../]
  [./s33]
    type = RankTwoAux
    variable = s33
    rank_two_tensor = PK1
    index_i = 2
    index_j = 2
  [../]

  [./F11]
    type = RankTwoAux
    variable = F11
    rank_two_tensor = def_grad
    index_i = 0
    index_j = 0
  [../]
  [./F21]
    type = RankTwoAux
    variable = F21
    rank_two_tensor = def_grad
    index_i = 1
    index_j = 0
  [../]
  [./F31]
    type = RankTwoAux
    variable = F31
    rank_two_tensor = def_grad
    index_i = 2
    index_j = 0
  [../]
  [./F12]
    type = RankTwoAux
    variable = F12
    rank_two_tensor = def_grad
    index_i = 0
    index_j = 1
  [../]
  [./F22]
    type = RankTwoAux
    variable = F22
    rank_two_tensor = def_grad
    index_i = 1
    index_j = 1
  [../]
  [./F32]
    type = RankTwoAux
    variable = F32
    rank_two_tensor = def_grad
    index_i = 2
    index_j = 1
  [../]
  [./F13]
    type = RankTwoAux
    variable = F13
    rank_two_tensor = def_grad
    index_i = 0
    index_j = 2
  [../]
  [./F23]
    type = RankTwoAux
    variable = F23
    rank_two_tensor = def_grad
    index_i = 1
    index_j = 2
  [../]
  [./F33]
    type = RankTwoAux
    variable = F33
    rank_two_tensor = def_grad
    index_i = 2
    index_j = 2
  [../]
[]

[UserObjects]
  [./integrator]
    type = HomogenizationConstraintIntegral
    targets = 'strain11'
    execute_on = 'initial linear'
  [../]
[]

[Kernels]
  [./sdx]
    type = TotalStressDivergenceNEML
    variable = disp_x
    component = 0
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
  [./strain11]
    type = ParsedFunction
    value = '4.0e-1*t'
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = disp_x
      auto_direction = 'x'
    [../]
  [../]

  [./centerfix_x]
    type = DirichletBC
    boundary = "fixme"
    variable = disp_x
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
  [./s11]
    type = ElementAverageValue
    variable = s11
    execute_on = 'initial timestep_end'
  [../]
  [./s21]
    type = ElementAverageValue
    variable = s21
    execute_on = 'initial timestep_end'
  [../]
  [./s31]
    type = ElementAverageValue
    variable = s31
    execute_on = 'initial timestep_end'
  [../]
  [./s12]
    type = ElementAverageValue
    variable = s12
    execute_on = 'initial timestep_end'
  [../]
  [./s22]
    type = ElementAverageValue
    variable = s22
    execute_on = 'initial timestep_end'
  [../]
  [./s32]
    type = ElementAverageValue
    variable = s32
    execute_on = 'initial timestep_end'
  [../]
  [./s13]
    type = ElementAverageValue
    variable = s13
    execute_on = 'initial timestep_end'
  [../]
  [./s23]
    type = ElementAverageValue
    variable = s23
    execute_on = 'initial timestep_end'
  [../]
  [./s33]
    type = ElementAverageValue
    variable = s33
    execute_on = 'initial timestep_end'
  [../]

  [./F11]
    type = ElementAverageValue
    variable = F11
    execute_on = 'initial timestep_end'
  [../]
  [./F21]
    type = ElementAverageValue
    variable = F21
    execute_on = 'initial timestep_end'
  [../]
  [./F31]
    type = ElementAverageValue
    variable = F31
    execute_on = 'initial timestep_end'
  [../]
  [./F12]
    type = ElementAverageValue
    variable = F12
    execute_on = 'initial timestep_end'
  [../]
  [./F22]
    type = ElementAverageValue
    variable = F22
    execute_on = 'initial timestep_end'
  [../]
  [./F32]
    type = ElementAverageValue
    variable = F32
    execute_on = 'initial timestep_end'
  [../]
  [./F13]
    type = ElementAverageValue
    variable = F13
    execute_on = 'initial timestep_end'
  [../]
  [./F23]
    type = ElementAverageValue
    variable = F23
    execute_on = 'initial timestep_end'
  [../]
  [./F33]
    type = ElementAverageValue
    variable = F33
    execute_on = 'initial timestep_end'
  [../]
  [./nonlin]
    type = NumNonlinearIterations
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
