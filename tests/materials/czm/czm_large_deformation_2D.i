[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    nx = 2
    ny = 2
    dim = 2
  []
  [./subdomain_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    block_id = 1
    top_right = '0.5 1 0'
  []
  [./subdomain_2]
    type = SubdomainBoundingBoxGenerator
    input = subdomain_1
    bottom_left = '0.5 0 0'
    block_id = 2
    top_right = '1 1 0'
  []
  [./breakmesh]
    input = subdomain_2
    type = BreakMeshByBlockGenerator
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[NEMLMechanics]
  kinematics = large
  add_all_output = true
  add_displacements = true
[]

[CohesiveZoneDeer]
   boundary = 'interface'
[]

[AuxVariables]
  [./TPK1_solid_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1_solid_Y]
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
  [./t_solid_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./t_solid_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_Y]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [./T_PK1_solid_x]
    type = TractionAux
    scalar_type = 'X'
    variable = TPK1_solid_X
    property = stress
    boundary = 'interface'
    PK1 = true
  []
  [./T_PK1_solid_y]
    type = TractionAux
    scalar_type = 'Y'
    variable = TPK1_solid_Y
    property = stress
    boundary = 'interface'
    PK1 = true
  []
  [./T_cauchy_solid_x]
    type = TractionAux
    scalar_type = 'X'
    variable = t_solid_X
    property = stress
    boundary = 'interface'
  []
  [./T_cauchy_solid_y]
    type = TractionAux
    scalar_type = 'Y'
    variable = t_solid_Y
    property = stress
    boundary = 'interface'
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
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./right_y]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '0.2*t'
  [../]
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "mat.xml"
    model = "elastic_model"
    large_kinematics = true
  [../]
  [./czm_mat]
    type = PureElasticCZMIncremental
    E = 1e7
    G = 1e6
    boundary = 'interface'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  solve_type = NEWTON
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt =  0.05
  end_time = 5
  dtmin = 0.05
  line_search = none
[]

[Outputs]
  csv = true
  exodus = false
[]

[Postprocessors]
  [./sxx]
    type = SideAverageValue
    variable = stress_x-x
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./syy]
    type = SideAverageValue
    variable = stress_y-y
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./sxy]
    type = SideAverageValue
    variable = stress_x-y
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./disp_x]
    type = SideAverageValue
    variable = disp_x
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'right'
  [../]
  [./disp_y]
    type = SideAverageValue
    variable = disp_y
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'right'
  [../]
  [./TPK1_solid_X]
    type = SideAverageValue
    variable = TPK1_solid_X
    boundary = interface
  [../]
  [./TPK1_solid_Y]
    type = SideAverageValue
    variable = TPK1_solid_Y
    boundary = interface
  [../]
  [./t_solid_X]
    type = SideAverageValue
    variable = t_solid_X
    boundary = interface
  [../]
  [./t_solid_Y]
    type = SideAverageValue
    variable = t_solid_Y
    boundary = interface
  [../]
  [./TPK1czm_X]
    type = SideAverageValue
    variable = TPK1czm_X
    boundary = interface
  [../]
  [./TPK1czm_Y]
    type = SideAverageValue
    variable = TPK1czm_Y
    boundary = interface
  [../]
  [./tczm_X]
    type = SideAverageValue
    variable = tczm_X
    boundary = interface
  [../]
  [./tczm_Y]
    type = SideAverageValue
    variable = tczm_Y
    boundary = interface
  [../]
[]
