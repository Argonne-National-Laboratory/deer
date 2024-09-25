[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
  [new_block]
    type = SubdomainBoundingBoxGenerator
    input = generated_mesh
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '1 1 1'
  []
  [boundary]
    type = SideSetsBetweenSubdomainsGenerator
    input = new_block
    primary_block = '0 1'
    paired_block = 1
    new_boundary = 'interface'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        new_system = true
        formulation = UPDATED
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "../test_materials.xml"
    model = "elastic_model"
    large_kinematics = false
  []
  [tensor_rate]
    type = TensorRateMaterial
    rank_two_tensor = cauchy_stress
  []
[]

[Functions]
  [topfunc_x]
    type = PiecewiseLinear
    x = '0 2'
    y = '0 10'
  []
  [topfunc_y]
    type = PiecewiseLinear
    x = '0 2'
    y = '0 20'
  []
  [topfunc_z]
    type = PiecewiseLinear
    x = '0 2'
    y = '0 -30'
  []
[]

[BCs]
  [x_0]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [y_0]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [z_0]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  [x_1]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = 'right'
    function = topfunc_x
  []
  [y_1]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = 'top'
    function = topfunc_y
  []
  [z_1]
    type = FunctionNeumannBC
    variable = disp_z
    boundary = 'front'
    function = topfunc_z
  []

[]

[AuxVariables]
  [sdot_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [sdot_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [sdot_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [sdot_xx]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress_rate
    variable = sdot_xx
    index_i = 0
    index_j = 0
  []
  [sdot_yy]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress_rate
    variable = sdot_yy
    index_i = 1
    index_j = 1
  []
  [sdot_zz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress_rate
    variable = sdot_zz
    index_i = 2
    index_j = 2
  []
[]

[Postprocessors]
  [sdot_xx]
    type = ElementAverageValue
    variable = sdot_xx
  []
  [sdot_yy]
    type = ElementAverageValue
    variable = sdot_yy
  []
  [sdot_zz]
    type = ElementAverageValue
    variable = sdot_zz
  []
  [s_xx]
    type = ElementAverageValue
    variable = cauchy_stress_xx
  []
  [s_yy]
    type = ElementAverageValue
    variable = cauchy_stress_yy
  []
  [s_zz]
    type = ElementAverageValue
    variable = cauchy_stress_zz
  []
  [sderivative_xx]
    type = TimeDerivativePostprocessor
    postprocessor = s_xx
  []
  [sderivative_yy]
    type = TimeDerivativePostprocessor
    postprocessor = s_yy
  []
  [sderivative_zz]
    type = TimeDerivativePostprocessor
    postprocessor = s_zz
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 4
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dtmin = 1
  dt = 1
  end_time = 2.0
[]
[Outputs]
  csv = true
[]
