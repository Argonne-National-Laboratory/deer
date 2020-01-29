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
  [./interface_damage]
    family = MONOMIAL
    order = CONSTANT
  []
  [./T_N]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[UserObjects]
  [./interface_damage]
    type = Map2LDInterfaceMPReal
    mp_name = interface_damage
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
  [./aux_TN]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_TN
    variable = T_N
    execute_on = 'TIMESTEP_END'
  [../]
  [./aux_interface_damage]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_damage
    variable = interface_damage
    execute_on = 'TIMESTEP_END'
  [../]
[]


[NEMLMechanics]
  displacements = "disp_x disp_y disp_z"
  kinematics = small
  add_all_output = true
[]

[InterfaceKernels]
  [./czmx]
    type = CZMInterfaceKernel
    displacements = 'disp_x disp_y disp_z'
    component = 0
    variable = disp_x
    neighbor_var = disp_x
    boundary = 'interface'
  []
  [./czmy]
    type = CZMInterfaceKernel
    displacements = 'disp_x disp_y disp_z'
    component = 1
    variable = disp_y
    neighbor_var = disp_y
    boundary = 'interface'
  []
  [./czmz]
    type = CZMInterfaceKernel
    displacements = 'disp_x disp_y disp_z'
    component = 2
    variable = disp_z
    neighbor_var = disp_z
    boundary = 'interface'
  []
[]

[BCs]
  [./x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./right]
    type = DirichletBC
    boundary = right
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
  [./z_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = disp_z_fun
  [../]
[]

[Functions]
  [./disp_z_fun]
    type = PiecewiseLinear
    x = '0 1   3600'
    y = '0 0.1 0.1'
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
    type = TimeDependentDamage
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    K_n = 1e4
    K_t = 1e4
    max_damage =.99
    stiffness_reduction_factor = 1000.
    residual_life_scaling_factor = 1.
    effective_stress_mp_name = vonmises_interface
    x = '0 1000'
    y = '1e3 1e2'
  [../]
  [./vonmises_interface]
    type = InterfaceEffectiveStress
    effective_stress_type = VonMises
    boundary = interface
    effective_stress_mp_name = vonmises_interface
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

  dtmin = 1
  dtmax = 20
  end_time = 3600
  [TimeStepper]
  type = IterationAdaptiveDT
  optimal_iterations = 10
  dt = 1
[]
[]

[Outputs]
  exodus = true
[]
