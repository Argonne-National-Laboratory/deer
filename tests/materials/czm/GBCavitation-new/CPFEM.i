a0 = 5e-5
b0 = 6e-2
D_GB = 5e-15
E = 1.5e5
G = 5.837e4
w = 0.0113842
eta_s = 1e6
T0 = 100
FN = 1.77e6
Nc = 90
Tx = 150
Ty = 0

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  large_kinematics = true
  volumetric_locking_correction = true
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 8
    ny = 4
    nz = 4
    xmax = 2
  []
  [left]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '1 1 1'
  []
  [right]
    type = SubdomainBoundingBoxGenerator
    input = left
    block_id = 2
    bottom_left = '1 0 0'
    top_right = '2 1 1'
  []
  [breakmesh]
    type = BreakMeshByBlockGenerator
    input = right
  []
  use_displaced_mesh = false
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = FINITE
        new_system = true
        add_variables = true
        formulation = TOTAL
        volumetric_locking_correction = false
      []
    []
    [CohesiveZoneMaster]
      [interface]
        boundary = 'interface'
        strain = FINITE
        generate_output = 'traction_x traction_y traction_z'
      []
    []
  []
[]

[Functions]
  [load_x]
    type = PiecewiseLinear
    x = '0 3600'
    y = '0 ${Tx}'
  []
  [load_y]
    type = PiecewiseLinear
    x = '0 3600'
    y = '0 ${Ty}'
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  [force_x]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = 'right'
    function = 'load_x'
  []
  [force_y]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = 'top'
    function = 'load_y'
  []
[]

[Constraints]
  [ev_x]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    secondary = 'right'
    penalty = 1e8
  []
  [ev_y]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = 'top'
    penalty = 1e8
  []
  [ev_z]
    type = EqualValueBoundaryConstraint
    variable = disp_z
    secondary = 'front'
    penalty = 1e8
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = 'test.xml'
    model = 'test'
  []
  [GB_props]
    type = GenericConstantMaterial
    prop_names = 'a0 b0 D_GB E G w eta_s T0 FN Nc'
    prop_values = '${a0} ${b0} ${D_GB} ${E} ${G} ${w} ${eta_s} ${T0} ${FN} ${Nc}'
    boundary = 'interface'
  []
  [GB]
    type = GrainBoundaryCavitation
    a0 = a0
    b0 = b0
    psi = 70
    n = 5
    P = 10
    gamma = 2
    eps = 1e-6
    boundary = 'interface'
  []
  [stress_xx]
    type = RankTwoCartesianComponent
    property_name = stress_xx
    rank_two_tensor = cauchy_stress
    index_i = 0
    index_j = 0
  []
  [strain_xx]
    type = RankTwoCartesianComponent
    property_name = strain_xx
    rank_two_tensor = mechanical_strain
    index_i = 0
    index_j = 0
  []
  [creep_strain_xx]
    type = RankTwoCartesianComponent
    property_name = creep_strain_xx
    rank_two_tensor = inelastic_strain
    index_i = 0
    index_j = 0
  []
[]

[Postprocessors]
  [a]
    type = SideAverageMaterialProperty
    property = average_cavity_radius
    boundary = interface
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [b]
    type = SideAverageMaterialProperty
    property = average_cavity_half_spacing
    boundary = interface
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [D]
    type = SideAverageMaterialProperty
    property = damage
    boundary = interface
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [avg_disp_x]
    type = SideAverageValue
    variable = disp_x
    boundary = right
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [strain]
    type = ParsedPostprocessor
    pp_names = 'avg_disp_x'
    function = 'avg_disp_x / 2'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [delta_strain]
    type = ChangeOverTimePostprocessor
    postprocessor = strain
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [dt]
    type = TimestepSize
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [strain_rate]
    type = ParsedPostprocessor
    pp_names = 'delta_strain dt'
    function = 'delta_strain / dt'
    execute_on = 'TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  line_search = none
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 12
  dt = 10
  end_time = 2e7
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    iteration_window = 2
    linear_iteration_ratio = 1000
    dt = 10
    growth_factor = 1.5
    cutback_factor = 0.5
    cutback_factor_at_failure = 0.2
  []
  [Predictor]
    type = SimplePredictor
    scale = 1
    skip_after_failed_timestep = true
  []
  dtmax = 1e4
[]

[Outputs]
  file_base = 'Tx_${Tx}_Ty_${Ty}'
  exodus = true
  csv = true
  print_linear_residuals = false
[]
