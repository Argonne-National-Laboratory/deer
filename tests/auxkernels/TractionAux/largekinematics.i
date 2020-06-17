# Check that PK1 traction is compute correctly and that traction rotation is performedd correctly

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

[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.5
  zmax = 0.5
  []
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

[AuxVariables]
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
  [./traction_N]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./tPK1_N]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./traction_S1]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./traction_S2]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./traction_X]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./traction_Y]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./traction_Z]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Functions]
  [./angles]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0 1.5707963'
  [../]

  [./stretch]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0.1 0.1'
  [../]
  [./dt_fun]
    type = PiecewiseConstant
    x = '0 0.9 2'
    y = '0.1 0.001 0.001'
  []

  [./move_y]
    type = ParsedFunction
    value = 'y*cos(theta) - z * (1 + a)*sin(theta) - y'
    vars = 'a theta'
    vals = 'stretch angles'
  [../]

  [./move_z]
    type = ParsedFunction
    value = 'y*sin(theta) + z*(1+a)*cos(theta) - z'
    vars = 'a theta'
    vals = 'stretch angles'
  [../]

  [./TPK1_check]
    type = ParsedFunction
    value = 'a/a0*tn - tPK1N'
    vars = 'a a0 tn tPK1N'
    vals = 'A A0 t_N tPK1_N'
  [../]
[]

[BCs]
  [./fix]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = left
    variable = disp_x
  [../]

  [./front_y]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_y
    function = move_y
    preset = true
  [../]

  [./back_y]
    type = FunctionDirichletBC
    boundary = back
    variable = disp_y
    function = move_y
    preset = true
  [../]

  [./front_z]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = move_z
    preset = true
  [../]

  [./back_z]
    type = FunctionDirichletBC
    boundary = back
    variable = disp_z
    function = move_z
    preset = true
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
    boundary = back
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
    boundary = back
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
    boundary = back
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
    boundary = back
  [../]
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
    execute_on = timestep_end
    boundary = back
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    execute_on = timestep_end
    boundary = back
  [../]
  [./traction_N]
    type = TractionAux
    scalar_type = 'normal'
    variable = traction_N
    property = stress
    boundary = back
  [../]
  [./tPK1_N]
    type = TractionAux
    scalar_type = 'normal'
    variable = tPK1_N
    property = stress
    boundary = back
    PK1 = true
  [../]
  [./traction_S1]
    type = TractionAux
    scalar_type = 'shear1'
    variable = traction_S1
    property = stress
    boundary = back
  [../]
  [./traction_S2]
    type = TractionAux
    scalar_type = 'shear2'
    variable = traction_S2
    property = stress
    boundary = back
  [../]
  [./traction_X]
    type = TractionAux
    scalar_type = 'X'
    variable = traction_X
    property = stress
    boundary = back
  [../]
  [./traction_Y]
    type = TractionAux
    scalar_type = 'Y'
    variable = traction_Y
    property = stress
    boundary = back
  [../]
  [./traction_Z]
    type = TractionAux
    scalar_type = 'Z'
    variable = traction_Z
    property = stress
    boundary = back
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

[Postprocessors]
  [./sxx]
    type = SideAverageValue
    variable = stress_xx
    boundary = back
  [../]
  [./syy]
    type = SideAverageValue
    variable = stress_yy
    boundary = back
  [../]
  [./szz]
    type = SideAverageValue
    variable = stress_zz
    boundary = back
  [../]
  [./syz]
    type = SideAverageValue
    variable = stress_yz
    boundary = back
  [../]
  [./sxz]
    type = SideAverageValue
    variable = stress_xz
    boundary = back
  [../]
  [./sxy]
    type = SideAverageValue
    variable = stress_xy
    boundary = back
  [../]
  [./t_N]
    type = SideAverageValue
    variable = traction_N
    boundary = back
  [../]
  [./t_S1]
    type = SideAverageValue
    variable = traction_S1
    boundary = back
  [../]
  [./t_S2]
    type = SideAverageValue
    variable = traction_S2
    boundary = back
  [../]
  [./t_X]
    type = SideAverageValue
    variable = traction_X
    boundary = back
  [../]
  [./t_Y]
    type = SideAverageValue
    variable = traction_Y
    boundary = back
  [../]
  [./t_Z]
    type = SideAverageValue
    variable = traction_Z
    boundary = back
  [../]
  [./tPK1_N]
    type = SideAverageValue
    variable = tPK1_N
    boundary = back
  [../]
  [./A]
    type = AreaPostprocessor
    boundary = back
    use_displaced_mesh = true
  [../]
  [./A0]
    type = AreaPostprocessor
    boundary = back
    use_displaced_mesh = false
  [../]
  [./TPK1_check]
    type = FunctionValuePostprocessor
    function = TPK1_check
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
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = 2.0
  [./TimeStepper]
    type = FunctionDT
    function = dt_fun
  [../]
[]

[Outputs]
  csv = true
[]
