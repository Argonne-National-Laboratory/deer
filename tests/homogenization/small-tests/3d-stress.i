# 2D test with just strain control

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  constraint_types = 'stress stress stress stress stress stress'
  ndim = 3
  large_kinematics = false
  macro_gradient = hvar
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

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
 [./hvar]
    family = SCALAR
    order = SIXTH
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
  [./szz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./syz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./sxz]
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
  [./ezz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./eyz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./exz]
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

  [./zz]
    type = RankTwoAux
    variable = szz
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
  [./syz]
    type = RankTwoAux
    variable = syz
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
  [../]
  [./sxz]
    type = RankTwoAux
    variable = sxz
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
  [../]

  [./exx]
    type = RankTwoAux
    variable = exx
    rank_two_tensor = mechanical_strain
    index_i = 0
    index_j = 0
  [../]
  [./eyy]
    type = RankTwoAux
    variable = eyy
    rank_two_tensor = mechanical_strain
    index_i = 1
    index_j = 1
  [../]
  [./exy]
    type = RankTwoAux
    variable = exy
    rank_two_tensor = mechanical_strain
    index_i = 0
    index_j = 1
  [../]

  [./ezz]
    type = RankTwoAux
    variable = ezz
    rank_two_tensor = mechanical_strain
    index_i = 2
    index_j = 2
  [../]
  [./eyz]
    type = RankTwoAux
    variable = eyz
    rank_two_tensor = mechanical_strain
    index_i = 1
    index_j = 2
  [../]
  [./exz]
    type = RankTwoAux
    variable = exz
    rank_two_tensor = mechanical_strain
    index_i = 0
    index_j = 2
  [../]
[]

[UserObjects]
  [./integrator]
    type = HomogenizationConstraintIntegral
    targets = 'stress11 stress22 stress33 stress23 stress13 stress12'
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
  [./sdz]
    type = TotalStressDivergenceNEML
    variable = disp_z
    component = 2
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
    value = '4.0e2*t'
  [../]
  [./stress22]
    type = ParsedFunction
    value = '-2.0e2*t'
  [../]
  [./stress33]
    type = ParsedFunction
    value = '8.0e2*t'
  [../]
  [./stress23]
    type = ParsedFunction
    value = '2.0e2*t'
  [../]
  [./stress13]
    type = ParsedFunction
    value = '-7.0e2*t'
  [../]
  [./stress12]
    type = ParsedFunction
    value = '1.0e2*t'
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

  [./szz]
    type = ElementAverageValue
    variable = szz
    execute_on = 'initial timestep_end'
  [../]
  [./syz]
    type = ElementAverageValue
    variable = syz
    execute_on = 'initial timestep_end'
  [../]
  [./sxz]
    type = ElementAverageValue
    variable = sxz
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

  [./ezz]
    type = ElementAverageValue
    variable = ezz
    execute_on = 'initial timestep_end'
  [../]
  [./eyz]
    type = ElementAverageValue
    variable = eyz
    execute_on = 'initial timestep_end'
  [../]
  [./exz]
    type = ElementAverageValue
    variable = exz
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
