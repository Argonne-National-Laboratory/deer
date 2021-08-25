[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -1
  zmax = 1
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-0.5 -0.5 0'
    top_right = '0.5 0.5 0.5'
  []
  [./scale]
  type = TransformGenerator
  input = new_block
  transform = SCALE
  vector_value ='0.06 0.06 0.06'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = scale
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./applied_load_x]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0'
  [../]
  [./applied_load_y]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0'
  [../]
  [./applied_load_z]
    type = PiecewiseLinear
    x = '0 0.1 300 300.02 1e6'
    y = '0 100 100 -500 -500'
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
    [./x1]
      type = FunctionNeumannBC
      boundary = right
      function = applied_load_x
      variable = disp_x
    [../]
    [./y1]
      type = FunctionNeumannBC
      boundary = top
      function = applied_load_y
      variable = disp_y
    [../]
    [./z1]
      type = FunctionNeumannBC
      boundary = front
      function = applied_load_z
      variable = disp_z
    [../]
[]

[AuxVariables]
  [./a]
    family = MONOMIAL
    order = CONSTANT
  []
  [./b]
    family = MONOMIAL
    order = CONSTANT
  []
  [./e]
    family = MONOMIAL
    order = CONSTANT
  []
  [./edot]
    family = MONOMIAL
    order = CONSTANT
  []
  [./SH_interface]
  family = MONOMIAL
  order = CONSTANT
  []
  [./SVM_interface]
    family = MONOMIAL
    order = CONSTANT
  []
  [./VLdot]
    family = MONOMIAL
    order = CONSTANT
  []
  [./VL1dot]
    family = MONOMIAL
    order = CONSTANT
  []
  [./VL2dot]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [./a]
    type = MaterialRealAux
    boundary = 'interface'
    property = average_cavity_radii
    execute_on = 'TIMESTEP_END'
    variable = a
    check_boundary_restricted = false
  []
  [./b]
    type = MaterialRealAux
    boundary = 'interface'
    property = average_cavity_spacing
    execute_on = 'TIMESTEP_END'
    variable = b
    check_boundary_restricted = false
  []
  [./e]
    type = MaterialRealAux
    boundary = 'interface'
    property = strain_eq_interface
    execute_on = 'TIMESTEP_END'
    variable = e
    check_boundary_restricted = false
  []
  [./edot]
    type = MaterialRealAux
    boundary = 'interface'
    property = strain_rate_eq_interface
    execute_on = 'TIMESTEP_END'
    variable = edot
    check_boundary_restricted = false
  []
  [./SH_interface]
    type = MaterialRealAux
    boundary = 'interface'
    property = stress_H_interface
    execute_on = 'TIMESTEP_END'
    variable = SH_interface
    check_boundary_restricted = false
  []
  [./SVM_interface]
    type = MaterialRealAux
    boundary = 'interface'
    property = stress_vm_interface
    execute_on = 'TIMESTEP_END'
    variable = SVM_interface
    check_boundary_restricted = false
  []
  [./VLdot]
    type = MaterialRealAux
    boundary = 'interface'
    property = VLdot
    execute_on = 'TIMESTEP_END'
    variable = VLdot
    check_boundary_restricted = false
  []
  [./VL1dot]
    type = MaterialRealAux
    boundary = 'interface'
    property = VL1dot
    execute_on = 'TIMESTEP_END'
    variable = VL1dot
    check_boundary_restricted = false
  []
  [./VL2dot]
    type = MaterialRealAux
    boundary = 'interface'
    property = VL2dot
    execute_on = 'TIMESTEP_END'
    variable = VL2dot
    check_boundary_restricted = false
  []
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
  formulation = total
[]


[Modules/TensorMechanics/CohesiveZoneMaster]
  [czm]
    strain = FINITE
    boundary = 'interface'
    generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
  []
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "mat.xml"
    model = "creep_and_hardening"
    large_kinematics = true
  [../]
  [./czm_mat]
    type = GBCavitation
    boundary = 'interface'
    max_time_cut = 10
    D_failure = 0.9
    max_nonlinear_iter = 20
    minimum_allowed_stiffness = 1
    minimum_allowed_residual_life = 10
    nucleation_on = true
    growth_on = true
    nl_residual_abs_tol = 1e-12
    vdot_method = 2
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  # Executioner
  type = Transient

  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type -snes_stol'
  petsc_options_value = 'lu 0'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_max_its = 2
  nl_max_its = 10
  start_time = 0.0
  [./TimeStepper]
   type =IterationAdaptiveDT
    optimal_iterations = 10
    dt =  0.1
  []
  dtmin = 1e-4
  end_time = 1000
  dtmax = 50
[]

[Outputs]
  exodus = true
  sync_times = '0 0.1 300 300.02 410 410.1'
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
[]
