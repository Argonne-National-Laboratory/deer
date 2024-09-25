# Simple 3D test

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 2
  []
  [new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '1 1 1'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = new_block
    primary_block = 0
    paired_block = 1
    new_boundary = 'interface'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = FINITE
        add_variables = true
        new_system = true
        formulation = UPDATED
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[AuxVariables]
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [zstress]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 500'
  []
  [constant]
    type = ConstantFunction
    value = 1.0
  []
  [szz_V0]
    type = ParsedFunction
    symbol_names = 'sd V0 V'
    symbol_values = 's_def volume0 volume'
    expression = 'sd *V / V0'
  []
  [szz_A0]
    type = ParsedFunction
    symbol_names = 'sd A A0'
    symbol_values = 's_def_interface area area0'
    expression = 'sd*A/A0'
  []
[]

[BCs]
  [leftx]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [boty]
    type = DirichletBC
    preset = true
    boundary = bottom
    variable = disp_y
    value = 0.0
  []
  [backz]
    type = DirichletBC
    preset = true
    boundary = back
    variable = disp_z
    value = 0.0
  []
  [pull_z]
    type = FunctionNeumannBC
    boundary = front
    variable = disp_z
    function = zstress
    use_displaced_mesh = true
  []
[]

[AuxKernels]
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "test.xml"
    model = "elastic_model"
    large_kinematics = true
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  #volume related functions
  [s_def]
    type = ElementAverageValue
    variable = stress_zz
    use_displaced_mesh = true
  []
  [s_def_V0]
    type = FunctionValuePostprocessor
    function = szz_V0
  []
  [volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
  []
  [volume0]
    type = VolumePostprocessor
    use_displaced_mesh = false
  []
  [s_def_new]
    type = MaterialTensorIntegralScaled
    rank_two_tensor = cauchy_stress
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
    scaling_factor_PP = volume
  []
  [s_undef_new]
    type = MaterialTensorIntegralScaled
    rank_two_tensor = cauchy_stress
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
    scaling_factor_PP = volume0
  []

  #interface related functions
  [s_def_interface]
    type = SideAverageValue
    variable = stress_zz
    boundary = interface
    use_displaced_mesh = true
  []
  [s_def_A0]
    type = FunctionValuePostprocessor
    function = szz_A0
  []
  [area]
    type = AreaPostprocessor
    boundary = interface
    use_displaced_mesh = true
  []
  [area0]
    type = AreaPostprocessor
    boundary = interface
    use_displaced_mesh = false
    execute_on = 'INITIAL'
  []
  [s_def_interface_new]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = cauchy_stress
    index_i = 2
    index_j = 2
    boundary = interface
    execute_on = 'TIMESTEP_END'
    scaling_factor_PP = area
    use_displaced_mesh = true
  []
  [s_undef_interface_new]
    type = MaterialTensorIntegralInterfaceScaled
    rank_two_tensor = cauchy_stress
    index_i = 2
    index_j = 2
    boundary = interface
    execute_on = 'TIMESTEP_END'
    scaling_factor_PP = area0
    use_displaced_mesh = true
  []
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Outputs]
  exodus = false
  csv = true
[]
