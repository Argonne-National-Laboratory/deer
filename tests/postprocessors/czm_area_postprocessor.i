[Mesh]
  [./msh]
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
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-1 -1.5 0'
    top_right = '1 1.5 1'
  []
  [./split]
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
  [./x0_bottom]
    type = ExtraNodesetGenerator
    input = add_side_sets
    new_boundary = 'x0_bottom'
    nodes = '0 3 4 7'
  []
  [./x1_bottom]
    type = ExtraNodesetGenerator
    input = x0_bottom
    new_boundary = 'x1_bottom'
    nodes = '1 2 5 6'
  []
  [./y0_bottom]
    type = ExtraNodesetGenerator
    input = x1_bottom
    new_boundary = 'y0_bottom'
    nodes = '0 1 4 5'
  []
  [./y1_bottom]
    type = ExtraNodesetGenerator
    input = y0_bottom
    new_boundary = 'y1_bottom'
    nodes = '2 3 6 7'
  []
  [./z1_bottom]
    type = ExtraNodesetGenerator
    input = y1_bottom
    new_boundary = 'z1_bottom'
    nodes = '4 5 6 7'
  []
  [./x0_top]
    type = ExtraNodesetGenerator
    input = z1_bottom
    new_boundary = 'x0_top'
    nodes = '8 11 12 15'
  []
  [./x1_top]
    type = ExtraNodesetGenerator
    input = x0_top
    new_boundary = 'x1_top'
    nodes = '9 10 13 14'
  []
  [./y0_top]
    type = ExtraNodesetGenerator
    input = x1_top
    new_boundary = 'y0_top'
    nodes = '8 9 12 13'
  []
  [./y1_top]
    type = ExtraNodesetGenerator
    input = y0_top
    new_boundary = 'y1_top'
    nodes = '10 11 14 15'
  []
  [./z0_top]
    type = ExtraNodesetGenerator
    input = y1_top
    new_boundary = 'z0_top'
    nodes = ' 12 13 14 15'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[NEMLMechanics]
  kinematics = small
  add_all_output = false
  add_displacements = true
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm]
    boundary = 'interface'
    strain = FINITE
  [../]
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "../neml_test_material.xml"
    model = "elastic_model"
    large_kinematics = false
  [../]
  [./czm]
    type = PureElasticCZM
    boundary = 'interface'
    E = 1e0
    G = 1e0
    interface_thickness = 1
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
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
  [./fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = x0
    variable = disp_x
  [../]
  [./fix_y]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = y0
    variable = disp_y
  [../]
  [./fix_z]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = 'z0 z1_bottom'
    variable = disp_z
  [../]
  [./z_stretch_top]
    type = FunctionDirichletBC
    preset = true
    function = t
    boundary = 'z1 z0_top'
    variable = disp_z
  [../]
  [./y_stretch]
    type = FunctionDirichletBC
    preset = true
    function = -1.5*t
    boundary = y1
    variable = disp_y
  [../]
  [./x_stretch]
    type = FunctionDirichletBC
    preset = true
    function = -1.*t
    boundary = x1
    variable = disp_x
  [../]

  [./rotate_x]
  type = DisplacementAboutAxis
  boundary = 'y0 y1 x0 x1 z0 z1'
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
  boundary = 'y0 y1 x0 x1 z0 z1'
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
  boundary = 'y0 y1 x0 x1 z0 z1'
  function = '90.'
  angle_units = degrees
  axis_origin = '0. 0. 0.'
  axis_direction = '0. 1. 0.'
  component = 2
  variable = disp_z
  angular_velocity = true
[../]
[]

[Controls]
  [./c1]
    type = TimePeriod
    enable_objects = 'BCs::fix_x BCs::fix_y BCs::fix_z BCs::z_stretch_top BCs::y_stretch BCs::x_stretch'
    disable_objects = 'BCs::rotate_x BCs::rotate_y BCs::rotate_z'
    start_time = '0'
    end_time = '1'
  [../]
[]

[Postprocessors]
  [./czm_A]
    type = CZMAreaPostprocessor
    strain = FINITE
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  []
  [./czm_daDA]
    type = CZMAreaRatioPostprocessor
    strain = FINITE
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  []
  [parsed]
    type = ParsedPostprocessor
    function = 'czm_A / czm_daDA'
    pp_names = 'czm_A czm_daDA'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]


[Outputs]
  csv = true
[]
