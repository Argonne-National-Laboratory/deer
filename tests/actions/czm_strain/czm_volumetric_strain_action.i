bulk_volume_PP = 'czm_strain_V0'

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 2
    xmin = -1
    xmax = 1
    ymin = -1.5
    ymax = 1.5
    zmin = -1
    zmax = 1
  []
  [new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-1 -1.5 0'
    top_right = '1 1.5 1'
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
  [add_side_sets]
    input = split
    type = SideSetsFromNormalsGenerator
    normals = '0 -1  0
               0  1  0
               -1 0  0
               1  0  0
               0  0 -1
               0  0  1'
    fixed_normal = true
    new_boundary = 'y0 y1 x0 x1 z0 z1'
  []
  [x0_bottom]
    type = ExtraNodesetGenerator
    input = add_side_sets
    new_boundary = 'x0_bottom'
    nodes = '0 3 4 7'
  []
  [x1_bottom]
    type = ExtraNodesetGenerator
    input = x0_bottom
    new_boundary = 'x1_bottom'
    nodes = '1 2 5 6'
  []
  [y0_bottom]
    type = ExtraNodesetGenerator
    input = x1_bottom
    new_boundary = 'y0_bottom'
    nodes = '0 1 4 5'
  []
  [y1_bottom]
    type = ExtraNodesetGenerator
    input = y0_bottom
    new_boundary = 'y1_bottom'
    nodes = '2 3 6 7'
  []
  [z1_bottom]
    type = ExtraNodesetGenerator
    input = y1_bottom
    new_boundary = 'z1_bottom'
    nodes = '4 5 6 7'
  []
  [x0_top]
    type = ExtraNodesetGenerator
    input = z1_bottom
    new_boundary = 'x0_top'
    nodes = '8 11 12 15'
  []
  [x1_top]
    type = ExtraNodesetGenerator
    input = x0_top
    new_boundary = 'x1_top'
    nodes = '9 10 13 14'
  []
  [y0_top]
    type = ExtraNodesetGenerator
    input = x1_top
    new_boundary = 'y0_top'
    nodes = '8 9 12 13'
  []
  [y1_top]
    type = ExtraNodesetGenerator
    input = y0_top
    new_boundary = 'y1_top'
    nodes = '10 11 14 15'
  []
  [z0_top]
    type = ExtraNodesetGenerator
    input = y1_top
    new_boundary = 'z0_top'
    nodes = ' 12 13 14 15'
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
        volumetric_locking_correction = false
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[Physics/SolidMechanics/CohesiveZone]
  [czm]
    boundary = 'interface'
    strain = FINITE
  []
[]

[Postprocessors]
  [czm_strain_V0]
    type = VolumePostprocessor
    block = '0 1'
    execute_on = INITIAL
    use_displaced_mesh = false
  []
  [czm_total_strain_xx]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 0
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_total_strain
    strain = FINITE
  []
  [czm_total_strain_xy]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 1
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_total_strain
    strain = FINITE
  []
  [czm_total_strain_xz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_total_strain
    strain = FINITE
  []
  [czm_total_strain_yy]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 1
    index_j = 1
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_total_strain
    strain = FINITE
  []
  [czm_total_strain_yz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 1
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_total_strain
    strain = FINITE
  []
  [czm_total_strain_zz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 2
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_total_strain
    strain = FINITE
  []
  [czm_normal_strain_xx]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 0
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_normal_strain
    strain = FINITE
  []
  [czm_normal_strain_xy]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 1
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_normal_strain
    strain = FINITE
  []
  [czm_normal_strain_xz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_normal_strain
    strain = FINITE
  []
  [czm_normal_strain_yy]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 1
    index_j = 1
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_normal_strain
    strain = FINITE
  []
  [czm_normal_strain_yz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 1
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_normal_strain
    strain = FINITE
  []
  [czm_normal_strain_zz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 2
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_normal_strain
    strain = FINITE
  []
  [czm_sliding_strain_xx]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 0
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_sliding_strain
    strain = FINITE
  []
  [czm_sliding_strain_xy]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 1
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_sliding_strain
    strain = FINITE
  []
  [czm_sliding_strain_xz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 0
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_sliding_strain
    strain = FINITE
  []
  [czm_sliding_strain_yy]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 1
    index_j = 1
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_sliding_strain
    strain = FINITE
  []
  [czm_sliding_strain_yz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 1
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_sliding_strain
    strain = FINITE
  []
  [czm_sliding_strain_zz]
    type = CZMStrainComponent
    boundary = interface
    execute_on = TIMESTEP_END
    index_i = 2
    index_j = 2
    initial_bulk_volume_pp = ${bulk_volume_PP}
    rank_two_tensor = czm_sliding_strain
    strain = FINITE
  []
  [czm_total_strain_rate_xx]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_total_strain_xx
  []
  [czm_total_strain_rate_xy]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_total_strain_xy
  []
  [czm_total_strain_rate_xz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_total_strain_xz
  []
  [czm_total_strain_rate_yy]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_total_strain_yy
  []
  [czm_total_strain_rate_yz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_total_strain_yz
  []
  [czm_total_strain_rate_zz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_total_strain_zz
  []
  [czm_normal_strain_rate_xx]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_normal_strain_xx
  []
  [czm_normal_strain_rate_xy]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_normal_strain_xy
  []
  [czm_normal_strain_rate_xz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_normal_strain_xz
  []
  [czm_normal_strain_rate_yy]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_normal_strain_yy
  []
  [czm_normal_strain_rate_yz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_normal_strain_yz
  []
  [czm_normal_strain_rate_zz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_normal_strain_zz
  []
  [czm_sliding_strain_rate_xx]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_sliding_strain_xx
  []
  [czm_sliding_strain_rate_xy]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_sliding_strain_xy
  []
  [czm_sliding_strain_rate_xz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_sliding_strain_xz
  []
  [czm_sliding_strain_rate_yy]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_sliding_strain_yy
  []
  [czm_sliding_strain_rate_yz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_sliding_strain_yz
  []
  [czm_sliding_strain_rate_zz]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_sliding_strain_zz
  []
  [czm_total_strain_eq]
    type = RankTwoTensorInvariantPostprocessor
    execute_on = TIMESTEP_END
    invariant = EffectiveStrain
    rank_two_tensor_base_name = czm_total_strain
  []
  [czm_normal_strain_eq]
    type = RankTwoTensorInvariantPostprocessor
    execute_on = TIMESTEP_END
    invariant = EffectiveStrain
    rank_two_tensor_base_name = czm_normal_strain
  []
  [czm_sliding_strain_eq]
    type = RankTwoTensorInvariantPostprocessor
    execute_on = TIMESTEP_END
    invariant = EffectiveStrain
    rank_two_tensor_base_name = czm_sliding_strain
  []
  [czm_total_strain_eq_rate]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_total_strain_eq
  []
  [czm_normal_strain_eq_rate]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_normal_strain_eq
  []
  [czm_sliding_strain_eq_rate]
    type = TimeDerivativePostprocessor
    execute_on = TIMESTEP_END
    postprocessor = czm_sliding_strain_eq
  []
[]

[Postprocessors]
  [V0]
    type = VolumePostprocessor
    use_displaced_mesh = false
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "../../neml_test_material.xml"
    model = "elastic_model"
    large_kinematics = false
  []
  [czm]
    type = PureElasticCZM
    boundary = 'interface'
    E = 1e0
    G = 1e0
    interface_thickness = 1
  []
  [czm_strain]
    type = CZMVolumetricStrain
    boundary = 'interface'
    strain = FINITE
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
  nl_abs_tol = 1e-6
  dt = 0.25
  end_time = 2
[]

[BCs]
  [fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = x0
    variable = disp_x
  []
  [fix_y]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = y0
    variable = disp_y
  []
  [fix_z]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = 'z0 z1_bottom'
    variable = disp_z
  []
  [z_stretch_top]
    type = FunctionDirichletBC
    preset = true
    function = t
    boundary = 'z1 z0_top'
    variable = disp_z
  []
  [y_stretch]
    type = FunctionDirichletBC
    preset = true
    function = -1.5*t
    boundary = y1
    variable = disp_y
  []
  [x_stretch]
    type = FunctionDirichletBC
    preset = true
    function = -1.*t
    boundary = x1
    variable = disp_x
  []
  [rotate_x]
    type = DisplacementAboutAxis
    boundary = 'y0 y1 x0 x1 z0 z1'
    function = '90.'
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 1. 0.'
    component = 0
    variable = disp_x
    angular_velocity = true
  []
  [rotate_y]
    type = DisplacementAboutAxis
    boundary = 'y0 y1 x0 x1 z0 z1'
    function = '90.'
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 1. 0.'
    component = 1
    variable = disp_y
    angular_velocity = true
  []
  [rotate_z]
    type = DisplacementAboutAxis
    boundary = 'y0 y1 x0 x1 z0 z1'
    function = '90.'
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 1. 0.'
    component = 2
    variable = disp_z
    angular_velocity = true
  []
[]

[Controls]
  [c1]
    type = TimePeriod
    enable_objects = 'BCs::fix_x BCs::fix_y BCs::fix_z BCs::z_stretch_top BCs::y_stretch BCs::x_stretch'
    disable_objects = 'BCs::rotate_x BCs::rotate_y BCs::rotate_z'
    start_time = '0'
    end_time = '1'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
