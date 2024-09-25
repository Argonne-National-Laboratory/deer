# Test under stress cotrolled condition BC. Load reversal is very quick

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 3
    xmin = -0.5
    xmax = 0.5
    ymin = -0.5
    ymax = 0.5
    zmin = 0
    zmax = 3
  []
  [new_block1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-0.5 -0.5 1'
    top_right = '0.5 0.5 2'
  []
  [new_block2]
    type = SubdomainBoundingBoxGenerator
    input = new_block1
    block_id = 2
    bottom_left = '-0.5 -0.5 2'
    top_right = '0.5 0.5 3'
  []
  [scale]
    type = TransformGenerator
    input = new_block2
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
    x = '0 0.1 300'
    y = '0 1 1'
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

# [AuxVariables]
#   [./a]
#     family = MONOMIAL
#     order = CONSTANT
#   []
#   [./b]
#     family = MONOMIAL
#     order = CONSTANT
#   []
# []

[UserObjects]
  [cavitation_properties]
    type = GBCavitationBoundaryPropertyUO
    file_name = 'cavitation_properties_test.txt'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    boundary = 'interface'
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
    nucleation_on = false
    growth_on = false
    vdot_method = 2
    GBCavitationBoundaryPropertyUO = 'cavitation_properties'
    output_properties = 'a0 a b n_exponent beta_exponent GB_diffusivity GB_sliding_viscosity NI FN GB_young_modulus GB_shear_modulus GB_thickness sigma_0 h'
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
  dt = 0.1
  dtmin = 0.1
  end_time = 0.1
  dtmax = 50
[]

[Outputs]
  exodus = true
  sync_times = '0 0.1'
[]
