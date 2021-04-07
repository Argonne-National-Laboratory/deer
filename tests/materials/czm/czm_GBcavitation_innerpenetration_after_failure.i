# Test under stress cotrolled condition BC. Load reversal is very quick

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
  [./applied_disp_z]
    type = PiecewiseLinear
    x = '0.0000000000000 1000 5000 10000'
    y = '0 0.06 0 -0.03'
  [../]
[]

[BCs]
  [./x0]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./y0]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./z0]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./z1]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function  = applied_disp_z
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

[AuxVariables]
  [./t_solid_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./t_solid_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./t_solid_Z]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_Z]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_N]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_S1]
    family = MONOMIAL
    order = CONSTANT
  []
  [./tczm_S2]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1_solid_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1_solid_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1_solid_Z]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_X]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_Y]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_Z]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_N]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_S1]
    family = MONOMIAL
    order = CONSTANT
  []
  [./TPK1czm_S2]
    family = MONOMIAL
    order = CONSTANT
  []

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
  [./u_N]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [./T_cauchy_solid_x]
    type = TractionAux
    scalar_type = 'X'
    variable = t_solid_X
    property = stress
    boundary = 'interface'
  []
  [./T_cauchy_solid_y]
    type = TractionAux
    scalar_type = 'Y'
    variable = t_solid_Y
    property = stress
    boundary = 'interface'
  []
  [./T_cauchy_solid_z]
    type = TractionAux
    scalar_type = 'Z'
    variable = t_solid_Z
    property = stress
    boundary = 'interface'
  []
  [./T_PK1_solid_x]
    type = TractionAux
    scalar_type = 'X'
    variable = TPK1_solid_X
    property = stress
    boundary = 'interface'
    PK1 = true
  []
  [./T_PK1_solid_y]
    type = TractionAux
    scalar_type = 'Y'
    variable = TPK1_solid_Y
    property = stress
    boundary = 'interface'
    PK1 = true
  []
  [./T_PK1_solid_z]
    type = TractionAux
    scalar_type = 'Z'
    variable = TPK1_solid_Z
    property = stress
    boundary = 'interface'
    PK1 = true
  []
  [./tczm_X]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction_deformed
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = tczm_X
  []
  [./tczm_Y]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction_deformed
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = tczm_Y
  []
  [./tczm_Z]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction_deformed
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = tczm_Z
  []
  [./tczm_N]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = tczm_N
  []
  [./tczm_S1]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = tczm_S1
  []
  [./tczm_S2]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = traction
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = tczm_S2
  []
  [./TPK1czm_X]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_X
  []
  [./TPK1czm_Y]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_Y
  []
  [./TPK1czm_Z]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_Z
  []
  [./TPK1czm_N]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction_natural
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_N
  []
  [./TPK1czm_S1]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction_natural
    component = 1
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_S1
  []
  [./TPK1czm_S2]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = PK1traction_natural
    component = 2
    execute_on = 'TIMESTEP_END'
    variable = TPK1czm_S2
  []
  [./a]
    type = MaterialRealAux
    boundary = 'interface'
    property = average_cavity_radii
    execute_on = 'TIMESTEP_END'
    variable = a
  []
  [./b]
    type = MaterialRealAux
    boundary = 'interface'
    property = average_cavity_spacing
    execute_on = 'TIMESTEP_END'
    variable = b
  []
  [./e]
    type = MaterialRealAux
    boundary = 'interface'
    property = strain_eq_interface
    execute_on = 'TIMESTEP_END'
    variable = e
  []
  [./edot]
    type = MaterialRealAux
    boundary = 'interface'
    property = strain_rate_eq_interface
    execute_on = 'TIMESTEP_END'
    variable = edot
  []
  [./SH_interface]
    type = MaterialRealAux
    boundary = 'interface'
    property = stress_H_interface
    execute_on = 'TIMESTEP_END'
    variable = SH_interface
  []
  [./SVM_interface]
    type = MaterialRealAux
    boundary = 'interface'
    property = stress_vm_interface
    execute_on = 'TIMESTEP_END'
    variable = SVM_interface
  []
  [./VLdot]
    type = MaterialRealAux
    boundary = 'interface'
    property = VLdot
    execute_on = 'TIMESTEP_END'
    variable = VLdot
  []
  [./VL1dot]
    type = MaterialRealAux
    boundary = 'interface'
    property = VL1dot
    execute_on = 'TIMESTEP_END'
    variable = VL1dot
  []
  [./VL2dot]
    type = MaterialRealAux
    boundary = 'interface'
    property = VL2dot
    execute_on = 'TIMESTEP_END'
    variable = VL2dot
  []
  [./U_N]
    type = MaterialRealVectorValueAux
    boundary = 'interface'
    property = displacement_jump
    component = 0
    execute_on = 'TIMESTEP_END'
    variable = u_N
  []
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
  formulation = total
[]


[CohesiveZoneDeer]
   boundary = 'interface'
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
    use_old_bulk_property = false
    nl_residual_abs_tol = 1e-12
    interface_thickness = 1e-4
    E_penalty_after_failure_minus_thickenss = 1e8
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
  automatic_scaling = true
  compute_scaling_once = false
  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5
  l_max_its = 2
  nl_max_its = 10
  start_time = 0.0
  [./TimeStepper]
   type =IterationAdaptiveDT
    optimal_iterations = 10
    dt =  0.1
  []
  dtmin = 1e-10
  end_time = 10000
  dtmax = 50
[]

[Outputs]
  exodus = true
  print_nonlinear_converged_reason = false
  print_linear_converged_reason = false
  # sync_times = '0 0.1 300 300.02 410 410.1 2000 3000'
[]
