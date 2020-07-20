# Test under stress cotrolled condition BC. Load reversal is very quick

[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -1
  zmax = 1
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-0.5 -0.5 0'
    top_right = '0.5 0.5 0.5'
  []
  [./scale]
  type = TransformGenerator
  input = new_block
  transform = SCALE
  vector_value ='0.06 0.06 0.06'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = scale
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./applied_load_x]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 -200 -200'
  [../]
  [./applied_load_y]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 100 100'
  [../]
  [./applied_load_z]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0'
  [../]
[]
[BCs]
    [./x]
      type = DirichletBC
      boundary = back
      variable = disp_x
      value = 0.0
    [../]
    [./y]
      type = DirichletBC
      boundary = back
      variable = disp_y
      value = 0.0
    [../]
    [./z]
      type = DirichletBC
      boundary = back
      variable = disp_z
      value = 0.0
    [../]
    [./x1]
      type = FunctionNeumannBC
      boundary = right
      function = applied_load_x
      variable = disp_x
    [../]
    [./y1]
      type = FunctionNeumannBC
      boundary = top
      function = applied_load_y
      variable = disp_y
    [../]
    [./z1]
      type = FunctionNeumannBC
      boundary = front
      function = applied_load_z
      variable = disp_z
    [../]
[]

[AuxVariables]
  [./tczm_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_Z]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_N]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_S1]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_S2]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_Z]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_N]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_S1]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_S2]
    family = MONOMIAL
    order = CONSTANT
  []

  [./a]
    family = MONOMIAL
    order = CONSTANT
  []
  [./b]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [./tczm_X]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction_deformed
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = tczm_X
  []
  [./tczm_Y]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction_deformed
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = tczm_Y
  []
  [./tczm_Z]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction_deformed
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = tczm_Z
  []
  [./tczm_N]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = tczm_N
  []
  [./tczm_S1]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = tczm_S1
  []
  [./tczm_S2]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = tczm_S2
  []
  [./TPK1czm_X]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_X
  []
  [./TPK1czm_Y]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_Y
  []
  [./TPK1czm_Z]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_Z
  []
  [./TPK1czm_N]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction_natural
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_N
  []
  [./TPK1czm_S1]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction_natural
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_S1
  []
  [./TPK1czm_S2]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction_natural
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_S2
  []
  [./a]
    type = MaterialRealAux
    boundary = 'interface'
    property = average_cavity_radii
    execute_on = 'TIMESTEP_END'
    variable = a
  []
  [./b]
    type = MaterialRealAux
    boundary = 'interface'
    property = average_cavity_spacing
    execute_on = 'TIMESTEP_END'
    variable = b
  []
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
  formulation = total
[]


[CohesiveZoneDeer]
   boundary = 'interface'
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "mat.xml"
    model = "creep_and_hardening"
    large_kinematics = true
  [../]
  [./czm_mat]
    type = GBCavitation
    boundary = 'interface'
    max_time_cut = 4
    D_failure = 0.9
    max_nonlinear_iter = 20
    minimum_allowed_stiffness = 1
    minimum_allowed_residual_life = 10
    nucleation_on = true
    growth_on = true
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  # Executioner
  type = Transient

  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_max_its = 2
  nl_max_its = 10
  start_time = 0.0
  [./TimeStepper]
   type =IterationAdaptiveDT
    optimal_iterations = 10
    dt =  0.1
  []
  dtmin = 1e-4
  end_time = 11
  dtmax = 50
[]

[Postprocessors]
  [a]
    type = SideAverageValue
    variable = a
    boundary = interface
  []
  [b]
    type = SideAverageValue
    variable = b
    boundary = interface
  []
  [t_X]
    type = SideAverageValue
    variable = tczm_X
    boundary = interface
  []
  [t_Y]
    type = SideAverageValue
    variable = tczm_Y
    boundary = interface
  []
  [t_Z]
    type = SideAverageValue
    variable = tczm_Z
    boundary = interface
  []
  [tpk1_X]
    type = SideAverageValue
    variable = TPK1czm_X
    boundary = interface
  []
  [tpk1_Y]
    type = SideAverageValue
    variable = TPK1czm_Y
    boundary = interface
  []
  [tpk1_Z]
    type = SideAverageValue
    variable = TPK1czm_Z
    boundary = interface
  []
[]

[Outputs]
  exodus = false
  csv = true
  sync_times = '0 0.1 10 10.0001 11'
[]
