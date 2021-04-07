[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmax = 3
  ymax = 2
  zmax = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 1'
    top_right = '3 2 2'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]



[NEMLMechanics]
  kinematics = large
  add_all_output = true
  add_displacements = true
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "neml_test_material.xml"
    model = "powerlaw"
    large_kinematics = true
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
  [./TimeStepper]
    type = FunctionDT
    function = dt_fun
  [../]
  end_time = 2
[]

[BCs]
  [./fix_x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./fix_y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./fix_z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./move_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = disp_fun
  [../]

  [./rotate_x]
  type = DisplacementAboutAxisDeer
  boundary = 'front back left right top bottom'
  function = '90.'
  angle_units = degrees
  axis_origin = '0. 0. 0.'
  axis_direction = '0. 1. 0.'
  component = 0
  variable = disp_x
  angular_velocity = true
[../]
[./rotate_y]
  type = DisplacementAboutAxisDeer
  boundary = 'front back left right top bottom'
  function = '90.'
  angle_units = degrees
  axis_origin = '0. 0. 0.'
  axis_direction = '0. 1. 0.'
  component = 1
  variable = disp_y
  angular_velocity = true
[../]
[./rotate_z]
  type = DisplacementAboutAxisDeer
  boundary = 'front back left right top bottom'
  function = '90.'
  angle_units = degrees
  axis_origin = '0. 0. 0.'
  axis_direction = '0. 1. 0.'
  component = 2
  variable = disp_z
  angular_velocity = true
[../]
[]

# Constraint System
[Constraints]
  [./x1]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    secondary = 'right'    # boundary
    penalty = 1e6
  [../]
  [./y1]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = 'top'    # boundary
    penalty = 1e6
  [../]
[]

[Controls]
  [./c1]
    type = TimePeriod
    enable_objects = 'BCs::fix_x BCs::fix_y BCs::fix_z BCs::move_top Constraints::x1 Constraints::y1'
    disable_objects = 'BCs::rotate_x BCs::rotate_y BCs::rotate_z'
    start_time = '0'
    end_time = '1'
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1.'
  [../]
  [./dt_fun]
    type = PiecewiseConstant
    x = '0 0.99 2'
    y = '0.01 0.001 0.001'
  []
[]

[Postprocessors]
  [e_xx_inelastic]
    type = MaterialTensorAverage
    rank_two_tensor = inelastic_strain_rotated
    index_i = 0
    index_j = 0
    use_displaced_mesh = true
  []
  [e_yy_inelastic]
    type = MaterialTensorAverage
    rank_two_tensor = inelastic_strain_rotated
    index_i = 1
    index_j = 1
    use_displaced_mesh = true
  []
  [e_zz_inelastic]
    type = MaterialTensorAverage
    rank_two_tensor = inelastic_strain_rotated
    index_i = 2
    index_j = 2
    use_displaced_mesh = true
  []
  [e_xy_inelastic]
    type = MaterialTensorAverage
    rank_two_tensor = inelastic_strain_rotated
    index_i = 0
    index_j = 1
    use_displaced_mesh = true
  []
  [e_xz_inelastic]
    type = MaterialTensorAverage
    rank_two_tensor = inelastic_strain_rotated
    index_i = 0
    index_j = 2
    use_displaced_mesh = true
  []
  [e_yz_inelastic]
    type = MaterialTensorAverage
    rank_two_tensor = inelastic_strain_rotated
    index_i = 1
    index_j = 2
    use_displaced_mesh = true
  []
  [e_xx_mechanical]
    type = MaterialTensorAverage
    rank_two_tensor = mechanical_strain_rotated
    index_i = 0
    index_j = 0
    use_displaced_mesh = true
  []
  [e_yy_mechanical]
    type = MaterialTensorAverage
    rank_two_tensor = mechanical_strain_rotated
    index_i = 1
    index_j = 1
    use_displaced_mesh = true
  []
  [e_zz_mechanical]
    type = MaterialTensorAverage
    rank_two_tensor = mechanical_strain_rotated
    index_i = 2
    index_j = 2
    use_displaced_mesh = true
  []
  [e_xy_mechanical]
    type = MaterialTensorAverage
    rank_two_tensor = mechanical_strain_rotated
    index_i = 0
    index_j = 1
    use_displaced_mesh = true
  []
  [e_xz_mechanical]
    type = MaterialTensorAverage
    rank_two_tensor = mechanical_strain_rotated
    index_i = 0
    index_j = 2
    use_displaced_mesh = true
  []
  [e_yz_mechanical]
    type = MaterialTensorAverage
    rank_two_tensor = mechanical_strain_rotated
    index_i = 1
    index_j = 2
    use_displaced_mesh = true
  []
  [e_xx_elastic]
    type = MaterialTensorAverage
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    use_displaced_mesh = true
  []
  [e_yy_elastic]
    type = MaterialTensorAverage
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    use_displaced_mesh = true
  []
  [e_zz_elastic]
    type = MaterialTensorAverage
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    use_displaced_mesh = true
  []
  [e_xy_elastic]
    type = MaterialTensorAverage
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    use_displaced_mesh = true
  []
  [e_xz_elastic]
    type = MaterialTensorAverage
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    use_displaced_mesh = true
  []
  [e_yz_elastic]
    type = MaterialTensorAverage
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    use_displaced_mesh = true
  []
[]

[Outputs]
  csv = true
  exodus=true
  sync_times = '0.1 0.5 1 2'
[]
