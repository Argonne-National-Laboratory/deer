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
    y = '0 0 0'
  [../]
  [./applied_load_y]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0'
  [../]
  [./applied_load_z]
    type = PiecewiseLinear
    x = '0 0.1 1e6'
    y = '0 100 100'
  [../]
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


[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm_ik]
    boundary = 'interface'
    strain = FINITE
    generate_output='traction_x traction_y traction_z jump_x jump_y jump_z normal_traction tangent_traction interface_jump_x interface_jump_y interface_jump_z normal_jump tangent_jump pk1_traction_x pk1_traction_y pk1_traction_z'
  [../]
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "mat.xml"
    model = "creep_and_hardening"
    large_kinematics = true
  [../]
  [./czm_mat]
    type = GBCavitationNewCZM
    boundary = 'interface'
    max_time_cut = 4
    D_failure = 0.9
    max_nonlinear_iter = 20
    minimum_allowed_stiffness = 1
    minimum_allowed_residual_life = 10
    nucleation_on = true
    growth_on = true
    vdot_method = 2
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
  end_time = 1e6
  dtmax = 1e4
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
[]

[Outputs]
  exodus = true
  csv = true
  # sync_times = '0 0.1 150 150.0001 300 300.02 410 410.1'
[]
