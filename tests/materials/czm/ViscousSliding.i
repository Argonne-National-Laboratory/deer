[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  zmax = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 1'
    top_right = '1 1 2'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[AuxVariables]
  [./T_N]
    family = MONOMIAL
    order = CONSTANT
  []
  [./T_S1]
    family = MONOMIAL
    order = CONSTANT
  []
  [./T_S2]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[NEMLMechanics]
  displacements = "disp_x disp_y disp_z"
  kinematics = small
  add_all_output = true
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm]
    boundary = 'interface'
    generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
  [../]
[]

[BCs]
  [./x]
    type = DirichletBC
    boundary = back
    variable = disp_x
    value = 0.0
  [../]
  [./y]
    type = DirichletBC
    boundary = back
    variable = disp_y
    value = 0.0
  [../]
  [./z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./z_top]
    type = DirichletBC
    boundary = front
    variable = disp_z
    value = 0.0
  [../]
  [./x_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_x
    function = disp_fun
  [../]
  [./y_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_y
    function = disp_fun_neg
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 1 10 11 20'
    y = '0 1 1 0 0'
  [../]
  [./disp_fun_neg]
    type = PiecewiseLinear
    x = '0 1 10 11 20'
    y = '0 -1.5 -1 -1 0'
  [../]
[]



[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "../../neml_test_material.xml"
    model = "elastic_model"
    large_kinematics = false
  [../]

  [./czm]
    type = ViscousSlidingCZM
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    E = 1e2
    G = 1e2
    shear_viscosity = 1e3
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
  [./TimeStepper]
   type =IterationAdaptiveDT
    dt =  1
  []
  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  dtmin = 0.5
  dtmax = 9
  end_time = 20
[]

[Outputs]
  exodus = true
  sync_times = '0 1 10 11 20'
[]
