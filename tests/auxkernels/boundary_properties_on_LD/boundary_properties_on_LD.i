[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    nx = 2
    ny = 2
    nz = 2
    dim = 3
  []
  [./subdomain_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    block_id = 1
    top_right = '0.5 1 1'
  []
  [./subdomain_2]
    type = SubdomainBoundingBoxGenerator
    input = subdomain_1
    bottom_left = '0.5 0 0'
    block_id = 2
    top_right = '1 1 1'
  []
  [./breakmesh]
    input = subdomain_2
    type = BreakMeshByBlockGenerator
  [../]
  [./lower]
     input = breakmesh
     type = LowerDBlockFromSidesetGenerator
     new_block_name = 'LD_interface'
     new_block_id = 1000
     sidesets = 6
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Kernels]
  [./dummy_x]
    type = NullKernel
    block = 'LD_interface'
    variable = disp_x
  []
  [./dummy_y]
    type = NullKernel
    block = 'LD_interface'
    variable = disp_y
  []
  [./dummy_z]
    type = NullKernel
    block = 'LD_interface'
    variable = disp_z
  []
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy'
    block = '1 2'
  [../]
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm1]
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[AuxVariables]
  [./T_N]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    preset = false
    boundary = left
    value = 0.0
  [../]
  [./left_y]
    type = DirichletBC
    variable = disp_y
    preset = false
    boundary = left
    value = 0.0
  [../]
  [./left_z]
    type = DirichletBC
    variable = disp_z
    preset = false
    boundary = left
    value = 0.0
  [../]
  [./right_x]
    type = FunctionDirichletBC
    variable = disp_x
    preset = false
    boundary = right
    function = '1*t'
  [../]
  [./right_y]
    type = FunctionDirichletBC
    variable = disp_y
    preset = false
    boundary = right
    function = '0*t'
  [../]
  [./right_z]
    type = FunctionDirichletBC
    variable = disp_z
    preset = false
    boundary = right
    function = '0*t'
  [../]
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
[]
[AuxKernels]
  [./aux_TN]
    type = Boundary2LDAux
    block = 1000
    map2LDelem_uo_name = interface_TN
    variable = T_N
    execute_on = 'TIMESTEP_END'
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2'
  [../]
  [./czm_3dc]
    type = SalehaniIrani3DCTraction
    boundary = 'interface'
    normal_gap_at_maximum_normal_traction = 1
    tangential_gap_at_maximum_shear_traction = 0.5
    maximum_normal_traction = 100
    maximum_shear_traction = 70
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./dummy]
    type =GenericConstantMaterial
    block = 1000
    prop_names = dummy
    prop_values = 0
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  solve_type = NEWTON
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 0.2
  end_time = 5
  dtmin = 0.2
  line_search = none
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
  sync_times = '0 1 2 3 4 5'
[]
