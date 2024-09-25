[Mesh]
  [msh]
    type = FileMeshGenerator
    file = test_10_gg_rcl1pt1_tet10.exo
  []
  [breakmesh]
    # split the mesh while preserving connectivity
    input = msh
    type = BreakMeshByBlockGenerator
  []
  [add_side_sets]
    # adding sidesets to apply boundary conditions
    input = breakmesh
    type = SideSetsFromNormalsGenerator
    normals = '0 -1  0
               0  1  0
               -1 0  0
               1  0  0
               0  0 -1
               0  0  1'
    fixed_normal = true
    new_boundary = 'y0 y1 x0 x1 z0 z1'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = FINITE
        new_system = true
        add_variables = true
        formulation = TOTAL
        volumetric_locking_correction = false
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

#output some material propertiy from the cavitation model
[AuxVariables]
  [a]
    family = MONOMIAL
    order = CONSTANT
  []
  [b]
    family = MONOMIAL
    order = CONSTANT
  []
  [D]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [a]
    type = MaterialRealAux
    boundary = 'interface'
    property = a
    execute_on = 'TIMESTEP_END'
    variable = a
    check_boundary_restricted = false #this is important
  []
  [b]
    type = MaterialRealAux
    boundary = 'interface'
    property = b
    execute_on = 'TIMESTEP_END'
    variable = b
    check_boundary_restricted = false #this is important
  []
  [D]
    type = MaterialRealAux
    boundary = 'interface'
    property = interface_damage
    execute_on = 'TIMESTEP_END'
    variable = D
    check_boundary_restricted = false #this is important
  []
[]

[Physics/SolidMechanics/CohesiveZone]
  [czm_ik]
    # add the proper cohesive interface kernels
    boundary = 'interface'
    strain = FINITE # use finite strins, total lagrangian formulation
    generate_output = 'traction_x traction_y traction_z jump_x jump_y jump_z normal_traction tangent_traction normal_jump tangent_jump' #output traction and jump
  []
[]

[Functions]
  # loading functions for each direction
  [applied_load_x]
    type = PiecewiseLinear
    x = '0 0.1 1e7' #time#
    y = '0 100 100' #PK1 stress in x direction
  []
  [applied_load_y]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 50 50' #PK1 stress in y direction
  []
  [applied_load_z]
    type = PiecewiseLinear
    x = '0 0.1 1e7'
    y = '0 0 0' #PK1 stress in z direction
  []
[]

[BCs]
  [x0]
    type = DirichletBC
    variable = disp_x
    boundary = x0
    value = 0.0
  []
  [y0]
    type = DirichletBC
    variable = disp_y
    boundary = y0
    value = 0.0
  []
  [z0]
    type = DirichletBC
    variable = disp_z
    boundary = z0
    value = 0.0
  []
  [x1]
    type = FunctionNeumannBC
    boundary = x1
    function = applied_load_x
    variable = disp_x
  []
  [y1]
    type = FunctionNeumannBC
    boundary = y1
    function = applied_load_y
    variable = disp_y
  []
  [z1]
    type = FunctionNeumannBC
    boundary = z1
    function = applied_load_z
    variable = disp_z
  []
[]

# Constraint System, used to keep the RVE faces flat thus imposing block periodic BC
[Constraints]
  [x1]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    secondary = 'x1'
    penalty = 1e6
  []
  [y1]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = 'y1'
    penalty = 1e6
  []
  [z1]
    type = EqualValueBoundaryConstraint
    variable = disp_z
    secondary = 'z1'
    penalty = 1e6
  []
[]

[UserObjects]
  [euler_angle_file]
    type = ElementPropertyReadFile
    nprop = 3
    prop_file_name = grn_10_rand.tex
    read_type = block
    nblock = 10
    use_zero_based_block_indexing = false
  []
[]
[Materials]
  [stress]
    # define the bulk material model, euler angles for each grain come from the `euler_angle_file` UserObjects
    type = NEMLCrystalPlasticity
    database = "Gr91.xml"
    model = "Gr91"
    large_kinematics = true
    euler_angle_reader = euler_angle_file
  []
  [ShamNeedleman]
    #setting up the GBcavitation model. Most model paramters have default values that have been calibrated for Grade 91 at 600C. If you want to modify them just list their value in this block. For the complete paramter list see GBCavitation.C
    type = GBCavitation
    boundary = 'interface' #the boundary on which we want to use the Cavitation model, this boundry is generate by BreakMeshByBlockGnerator
    # parameter below are optional
    a0 = 4e-5
    b0 = 5.9e-2
    D_failure = 0.9 # a/b ratio at which the traction decay model kicks in
  []
[]
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  solve_type = NEWTON
  line_search = none
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6
  nl_max_its = 20
  l_tol = 1e-15
  l_max_its = 1
  start_time = 0.0
  dtmin = 1e-4
  dtmax = 1e3
  end_time = 10
  n_max_nonlinear_pingpong = 1
  nl_forced_its = 2

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 2
    cutback_factor = 0.5
    optimal_iterations = 8
    iteration_window = 2
  []
[]

[Outputs]
  sync_times = '0 0.1 1 10'
  exodus = true
[]
