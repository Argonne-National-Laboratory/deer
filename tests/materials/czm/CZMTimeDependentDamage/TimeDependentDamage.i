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
    type = PresetBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./right]
    type = PresetBC
    boundary = right
    variable = disp_x
    value = 0.0
  [../]
  [./y]
    type = PresetBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./z]
    type = PresetBC
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
    database = "test.xml"
    model = "elastic_model"
    large_kinematics = false
  [../]

  [./czm]
    type = TimeDependentDamage
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    K_n = 1e4
    K_t = 1e4
    max_damage =.9
    stiffness_reduction_factor = 100.
    residual_life_scaling_factor = 10.
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
