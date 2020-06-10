[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = generated_mesh
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '1 1 1'
  []
  [./boundary]
    type = SideSetsBetweenSubdomainsGenerator
    input = new_block
    master_block = '0 1'
    paired_block = 1
    new_boundary = 'interface'
  []
[]

[NEMLMechanics]
  displacements = "disp_x disp_y disp_z"
  kinematics = small
  add_all_output = true
  add_displacements = true
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "../../test_materials.xml"
    model = "elastic_model"
    large_kinematics = false
  [../]
  [./tensor_rate]
    type = TensorRateMaterial
    rank_two_tensor = stress
   []
[]

[Functions]
 [./topfunc_x]
   type = PiecewiseLinear
   x = '0 2'
   y = '0 10'
 [../]
 [./topfunc_y]
   type = PiecewiseLinear
   x = '0 2'
   y = '0 20'
 [../]
 [./topfunc_z]
   type = PiecewiseLinear
   x = '0 2'
   y = '0 -30'
 [../]
[]

[BCs]
  [./x_0]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  [../]
  [./y_0]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  [./z_0]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  [../]
  [./x_1]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = 'right'
    function = topfunc_x
  [../]
  [./y_1]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = 'top'
    function = topfunc_y
  [../]
  [./z_1]
    type = FunctionNeumannBC
    variable = disp_z
    boundary = 'front'
    function = topfunc_z
  [../]

[]

[AuxVariables]
  [./sdot_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sdot_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sdot_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./sdot_xx]
    type = RankTwoAux
    rank_two_tensor = stress_rate
    variable = sdot_xx
    index_i = 0
    index_j = 0
  [../]
  [./sdot_yy]
    type = RankTwoAux
    rank_two_tensor = stress_rate
    variable = sdot_yy
    index_i = 1
    index_j = 1
  [../]
  [./sdot_zz]
    type = RankTwoAux
    rank_two_tensor = stress_rate
    variable = sdot_zz
    index_i = 2
    index_j = 2
  [../]
[]

[Postprocessors]
  [./sdot_xx]
    type = ElementAverageValue
    variable = sdot_xx
  [../]
  [./sdot_yy]
    type = ElementAverageValue
    variable = sdot_yy
  [../]
  [./sdot_zz]
    type = ElementAverageValue
    variable = sdot_zz
  [../]
  [./s_xx]
    type = ElementAverageValue
    variable = stress_x-x
  [../]
  [./s_yy]
    type = ElementAverageValue
    variable = stress_y-y
  [../]
  [./s_zz]
    type = ElementAverageValue
    variable = stress_z-z
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
