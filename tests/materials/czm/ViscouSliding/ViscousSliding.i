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
  [./lower]
     input = split
     type = LowerDBlockFromSidesetGenerator
     new_block_name = 'LD_interface'
     new_block_id = 1000
     sidesets = 6
  [../]
[]

[AuxVariables]
  [./dummy]
    family = MONOMIAL
    order = CONSTANT
  []
  [./interface_damage]
    family = MONOMIAL
    order = CONSTANT
  []
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

[UserObjects]
  [./interface_TN]
    type = Map2LDInterfaceMPRealVectorValue
    mp_name = traction
    component = 0
    ld_block_names = 'LD_interface'
    boundary = 'interface'
    execute_on = 'TIMESTEP_END'
  []
  [./interface_TS1]
    type = Map2LDInterfaceMPRealVectorValue
    mp_name = traction
    component = 1
    ld_block_names = 'LD_interface'
    boundary = 'interface'
    execute_on = 'TIMESTEP_END'
  []
  [./interface_TS2]
    type = Map2LDInterfaceMPRealVectorValue
    mp_name = traction
    component = 2
    ld_block_names = 'LD_interface'
    boundary = 'interface'
    execute_on = 'TIMESTEP_END'
  []
[]


[AuxKernels]
  [./aux_TN]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_TN
    variable = T_N
    execute_on = 'TIMESTEP_END'
  [../]
  [./aux_TS1]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_TS1
    variable = T_S1
    execute_on = 'TIMESTEP_END'
  [../]
  [./aux_TS2]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_TS2
    variable = T_S2
    execute_on = 'TIMESTEP_END'
  [../]
[]


[NEMLMechanics]
  displacements = "disp_x disp_y disp_z"
  kinematics = small
  add_all_output = true
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm]
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
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
    x = '0 1 1000 1001 2000'
    y = '0 1 1 0 0'
  [../]
  [./disp_fun_neg]
    type = PiecewiseLinear
    x = '0 1 1000 1200 2000'
    y = '0 -1.5 -1 -1 0'
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
  dtmax = 500
  end_time = 2000
[]

[Outputs]
  exodus = true
  sync_times = '1 100 200 300 400 500 600 700 800 900 1000 1001 1200 1300 1400 1500 1600 1700 1800 1900 2000'
[]
