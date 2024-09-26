# Test under stress cotrolled condition BC. Load reversal is very quick

[Mesh]
  [msh]
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
  [new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-0.5 -0.5 0'
    top_right = '0.5 0.5 0.5'
  []
  [scale]
    type = TransformGenerator
    input = new_block
    transform = SCALE
    vector_value = '0.06 0.06 0.06'
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = scale
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [applied_load_x]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0'
  []
  [applied_load_y]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0'
  []
  [applied_load_z]
    type = PiecewiseLinear
    x = '0 0.1 300 300.02 410 410.1'
    y = '0 100 100 -200 -200 0'
  []
[]
[BCs]
  [x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  []
  [z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  []
  [x1]
    type = FunctionNeumannBC
    boundary = right
    function = applied_load_x
    variable = disp_x
  []
  [y1]
    type = FunctionNeumannBC
    boundary = top
    function = applied_load_y
    variable = disp_y
  []
  [z1]
    type = FunctionNeumannBC
    boundary = front
    function = applied_load_z
    variable = disp_z
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = FINITE
        add_variables = true
        new_system = true
        formulation = TOTAL
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[Physics/SolidMechanics/CohesiveZone]
  [czm]
    strain = FINITE
    boundary = 'interface'
    generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "mat.xml"
    model = "creep_and_hardening"
    large_kinematics = true
  []
  [czm_mat]
    type = GBCavitation
    boundary = 'interface'
    max_time_cut = 4
    D_failure = 0.9
    max_nonlinear_iter = 20
    minimum_allowed_stiffness = 1
    minimum_allowed_residual_life = 10
    nucleation_on = true
    growth_on = true
    vdot_method = 2
    output_properties = 'a b'
    outputs = exodus
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
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
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    dt = 0.1
  []
  dtmin = 1e-4
  end_time = 410.1
  dtmax = 50
[]

[Outputs]
  exodus = true
  sync_times = '0 0.1 150 150.0001 198.9001 300 300.02 301.1199 303.3197 307.7193 316.5185 334.1169 369.3137 410 410.1'
[]
