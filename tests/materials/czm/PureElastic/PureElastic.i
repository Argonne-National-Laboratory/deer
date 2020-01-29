[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  zmax = 2
  elem_type = HEX20
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
     sidesets = '6'
  [../]
[]

[AuxVariables]
  [./T_N]
    family = MONOMIAL
    order = CONSTANT
  []
  [./stress_TN]
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
  [./bound_prop_const]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[UserObjects]
  [./interface_const]
    type = Map2LDInterfaceMPReal
    mp_name = bound_prop_const
    ld_block_names = 'LD_interface'
    boundary = 'interface'
    execute_on = 'TIMESTEP_END'
  []
  [./interface_TN]
    type = Map2LDInterfaceMPRealVectorValue
    mp_name = traction
    component = 0
    ld_block_names = 'LD_interface'
    boundary = 'interface'
    execute_on = 'TIMESTEP_END'
  []
[]

[AuxKernels]
  [./aux_normal_traction]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    variable = T_N
    component = 0
  [../]
  [./aux_shear_traction_1]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    variable = T_S1
    component = 1
  [../]
  [./aux_shear_traction_2]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    variable = T_S2
    component = 2
  [../]
  [./aux_const_material]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_const
    variable = bound_prop_const
    execute_on = 'TIMESTEP_END'
  [../]
  [./aux_TN]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_TN
    variable = T_N
    execute_on = 'TIMESTEP_END'
  [../]
  [./stress_TN]
    type = TractionAux
    boundary = 'interface'
    variable = stress_TN
    scalar_type = 'normal'
    property = stress
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
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[BCs]
  [./x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  # [./x_top]
  #   type = FunctionDirichletBC
  #   boundary = front
  #   variable = disp_x
  #   function = disp_fun
  # [../]
  # [./y_top]
  #   type = FunctionDirichletBC
  #   boundary = front
  #   variable = disp_y
  #   function = disp_fun
  # [../]
  [./x_top]
    type = DirichletBC
    boundary = front
    variable = disp_x
    value = 0.0
  [../]
  [./y_top]
    type = DirichletBC
    boundary = front
    variable = disp_y
    value = 0.0
  [../]
  [./z_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = disp_fun
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 1  3  4 6'
    y = '0 1 -1 -2 0'
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
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    E = 1e1
    G = 5e0
    penetration_penalty = 100
    interface_thickness = 0.1
  [../]
  [./GonsantBoundaryMaterial]
    type = GenericConstantMaterial
    boundary = 'interface back'
    prop_names = 'bound_prop_const'
    prop_values = 100
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
  dtmin = 1
  dtmax = 100
  end_time = 6
[]

[Outputs]
  exodus = true
  sync_times = '1 2 4 6'
[]
