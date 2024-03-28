[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1
    nx = 6
    ny = 6
    nz = 6
    elem_type = HEX20
  []
[]

[Variables]
    [p]
        order = FIRST
        family = LAGRANGE
    []
    [u]
        order = SECOND
        family = LAGRANGE_VEC
    []
[]

[AuxVariables]
    [u_exact]
        order = FIRST
        family = LAGRANGE_VEC
    []
    [p_exact]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxKernels]
    [u_exact]
        type = VectorFunctionAux
        variable = u_exact
        function = u_exact
    []
    [p_exact]
        type = FunctionAux
        variable = p_exact
        function = p_exact
    []
[]

[Functions]
    [u_exact]
        type = ParsedVectorFunction
        expression_x = '2 * x^2 + y^2 + z^2'
        expression_y = '2 * x^2 - 2 * x * y'
        expression_z = '2 * x^2 - 2 * x * z'
    []
    [p_exact]
        type = ParsedFunction
        expression = 'x + y + z - 3.0/2.0'
    []
[]

[ICs]
  [u]
    type = VectorConstantIC
    variable = u
    x_value = 1e-15
    y_value = 1e-15
    z_value = 1e-15
  []
[]

[Kernels]
    [equil]
        type = ADStokesStressDivergence
        variable = u
        pressure = p
    []
    [incompressible]
        type = ADStokesIncompressibility
        variable = p
        velocity = u
    []
    [body]
        type = VectorBodyForce
        variable = u
        function_x = 7
        function_y = 3
        function_z = 3
    []
[]

[Materials]
    [strain]
        type = StokesStrainRate
        velocity = u
    []
    [stress]
        type = StokesLinearViscous
        mu = 1.0
    []
[]

[BCs]
    [fixed]
        type = ADVectorFunctionDirichletBC
        variable = u
        boundary = 'bottom top back front'
        function_x = '2 * x^2 + y^2 + z^2'
        function_y = '2 * x^2 - 2 * x * y'
        function_z = '2 * x^2 - 2 * x * z'
    []
#    [natural]
#        type = ADVectorFunctionNeumannBC
#        variable = u
#        boundary = 'right'
#        function_x = 7
#        function_y = 3
#        function_z = 3
#    []
#    [pressure]
#        type = ADFunctionDirichletBC
#        boundary = 'left'
#        variable = p
#        function = p_exact
#    []
[]

[Preconditioning]
  [FSP]
    type = FSP
    topsplit = 'up'
    [up]
      splitting = 'u p'
      splitting_type  = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_type -ksp_pc_side'
      petsc_options_value = 'full                            self                              300                fgmres    right'
    []
    [u]
        vars = 'u'
        petsc_options_iname = '-pc_type'
        petsc_options_value = 'lu' 
    []
    [p]
        vars = 'p'
        petsc_options = '-pc_lsc_scale_diag'
        petsc_options_iname = '-ksp_type -ksp_gmres_restart -pc_type -ksp_pc_side -lsc_pc_type -lsc_ksp_type -lsc_ksp_pc_side -lsc_ksp_gmres_restart'
        petsc_options_value = 'fgmres    300                lsc      right        lu            gmres        right            300'
    []
  []
[]

[Executioner]
  type = Steady

  solve_type = 'newton'
  line_search = none

  l_max_its = 10
  l_tol = 1e-8
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
[]

[Postprocessors]
    [u_error]
        type = ElementVectorL2Error
        function = u_exact
        variable = u
    []
    [p_error]
        type = ElementL2Error
        function = p_exact
        variable = p
    []
[]

[Outputs]
    exodus = true
[]
