# 1D strain controlled test

[GlobalParams]
  displacements = 'disp_x'
  constraint_types = 'strain'
  ndim = 1
  large_kinematics = false
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
  [./sxx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./exx]
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
  [./exx]
    type = RankTwoAux
    variable = exx
    rank_two_tensor = mechanical_strain
    index_i = 0
    index_j = 0
  [../]
[]

[UserObjects]
  [./integrator]
    type = HomogenizationConstraintIntegral
    targets = 'func_strain'
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
  [./func_stress]
    type = ParsedFunction
    value = '100*t'
  [../]
  [./func_strain]
    type = ParsedFunction
    value = '1.0e-2*t'
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
  [./sxx]
    type = ElementAverageValue
    variable = sxx
    execute_on = 'initial timestep_end'
  [../]
  [./exx]
    type = ElementAverageValue
    variable = exx
    execute_on = 'initial timestep_end'
  [../]

  [./nonlin]
    type = NumNonlinearIterations
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = default

  automatic_scaling = true

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  start_time = 0.0
  dt = 0.2
  dtmin = 0.2
  end_time = 1.0
[]

[Outputs]
  exodus = false
  csv = true
[]
