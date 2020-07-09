#
# Stretch + rotation test
#
# This test is designed to compute a uniaxial stress and then follow that
# stress as the mesh is rotated 90 degrees.
#
# The mesh is composed of two blocks, each with a single element.

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
    y = '0 0 0'
  [../]
  [./applied_load_y]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0'
  [../]
  [./applied_load_z]
    type = PiecewiseLinear
    x = '0 0.1 300 300.02  1e7'
    y = '0 180 180 -180 -180'
  [../]
  [./dt_fun]
    type = PiecewiseConstant
    x = '0 0.99 2'
    y = '0.01 0.001 0.001'
  []
[]
[BCs]
    [./x]
      type = DirichletBC
      boundary = left
      variable = disp_x
      value = 0.0
    [../]
    [./y]
      type = DirichletBC
      boundary = bottom
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
  [./t_solid_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./t_solid_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./t_solid_Z]
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
  [./TPK1_solid_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1_solid_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1_solid_Z]
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
  [./e]
    family = MONOMIAL
    order = CONSTANT
  []
  [./edot]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
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
  [./T_cauchy_solid_z]
    type = TractionAux
    scalar_type = 'Z'
    variable = t_solid_Z
    property = stress
    boundary = 'interface'
  []
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
  [./T_PK1_solid_z]
    type = TractionAux
    scalar_type = 'Z'
    variable = TPK1_solid_Z
    property = stress
    boundary = 'interface'
    PK1 = true
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
  [./e]
    type = MaterialRealAux
    boundary = 'interface'
    property = strain_eq_interface
    execute_on = 'TIMESTEP_END'
    variable = e
  []
  [./edot]
    type = MaterialRealAux
    boundary = 'interface'
    property = strain_rate_eq_interface
    execute_on = 'TIMESTEP_END'
    variable = edot
  []
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
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
    scale_variables = false
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_z-z
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_y-y
  [../]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_x-x
  [../]
  [./stress_xy]
    type = ElementAverageValue
    variable = stress_x-y
  [../]
  [./stress_xz]
    type = ElementAverageValue
    variable = stress_x-z
  [../]
  [./stress_yz]
    type = ElementAverageValue
    variable = stress_y-z
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
  [./t_solid_Z]
    type = SideAverageValue
    variable = t_solid_Z
    boundary = interface
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
  [./TPK1_solid_Z]
    type = SideAverageValue
    variable = TPK1_solid_Z
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
  [./tczm_Z]
    type = SideAverageValue
    variable = tczm_Z
    boundary = interface
  [../]
  [./tczm_N]
    type = SideAverageValue
    variable = tczm_N
    boundary = interface
  [../]
  [./tczm_S1]
    type = SideAverageValue
    variable = tczm_S1
    boundary = interface
  [../]
  [./tczm_S2]
    type = SideAverageValue
    variable = tczm_S2
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
  [./TPK1czm_Z]
    type = SideAverageValue
    variable = TPK1czm_Z
    boundary = interface
  [../]
  [./TPK1czm_N]
    type = SideAverageValue
    variable = TPK1czm_N
    boundary = interface
  [../]
  [./TPK1czm_S1]
    type = SideAverageValue
    variable = TPK1czm_S1
    boundary = interface
  [../]
  [./TPK1czm_S2]
    type = SideAverageValue
    variable = TPK1czm_S2
    boundary = interface
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
  nl_abs_tol = 1e-6
  l_max_its = 2
  nl_max_its = 10
  start_time = 0.0
  [./TimeStepper]
   type =IterationAdaptiveDT
    dt =  0.1
  []
  end_time =1200
  dtmax = 1000
[]

[Outputs]
  csv = true
  exodus = true
[]
