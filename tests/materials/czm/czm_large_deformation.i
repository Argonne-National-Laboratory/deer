#
# Rotation Test
#
# This test is designed to compute a uniaxial stress and then follow that
# stress as the mesh is rotated 90 degrees.
#
# The mesh is composed of one block with a single element.  The nodal
# displacements in the x and y directions are prescribed.  Poisson's
# ratio is zero.
#

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
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./stretch]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1'
  [../]
  [./dt_fun]
    type = PiecewiseConstant
    x = '0 0.9 2'
    y = '0.01 0.001 0.001'
  []
[]
[BCs]
  [./fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = left
    variable = disp_x
  [../]
  [./fix_y]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = bottom
    variable = disp_y
  [../]
  [./fix_z]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = back
    variable = disp_z
  [../]
  [./back_z]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = stretch
    preset = true
  [../]

  [./rotate_x]
    type = DisplacementAboutAxis
    boundary = 'front back'
    function = '90.'
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 1. 0.'
    component = 0
    variable = disp_x
    angular_velocity = true
  [../]
  [./rotate_y]
    type = DisplacementAboutAxis
    boundary = 'front back'
    function = '90.'
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 1. 0.'
    component = 1
    variable = disp_y
    angular_velocity = true
  [../]
  [./rotate_z]
    type = DisplacementAboutAxis
    boundary = 'front back'
    function = '90.'
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 1. 0.'
    component = 2
    variable = disp_z
    angular_velocity = true
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
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
[]


[CohesiveZoneDeer]
   boundary = 'interface'
   displacements = 'disp_x disp_y disp_z'
[]

[Controls]
  [./c1]
    type = TimePeriod
    enable_objects = 'BCs::fix_x BCs::fix_y BCs::fix_z BCs::back_z'
    disable_objects = 'BCs::rotate_x BCs::rotate_y BCs::rotate_z'
    start_time = '0'
    end_time = '1'
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
  nl_rel_tol = 1e-30
  nl_abs_tol = 1e-8
  l_max_its = 20
  start_time = 0.0
  [./TimeStepper]
    type = FunctionDT
    function = dt_fun
  [../]
  end_time = 2
[]

[Outputs]
  csv = true
  exodus = true
[]
