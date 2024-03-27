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
    nx = 5
    ny = 5
    nz = 5
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
    [stress_xx]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_yy]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_zz]
        family = MONOMIAL
        order = CONSTANT
    []
[]

[AuxKernels]
    [stress_xx]
        type = ADMaterialRankTwoTensorAux
        variable = stress_xx
        property = 'stress'
        i = 0
        j = 0
    []
    [stress_yy]
        type = ADMaterialRankTwoTensorAux
        variable = stress_yy
        property = 'stress'
        i = 1
        j = 1
    []
    [stress_zz]
        type = ADMaterialRankTwoTensorAux
        variable = stress_zz
        property = 'stress'
        i = 2
        j = 2
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
    [no_slip]
        type = ADVectorFunctionDirichletBC
        variable = u
        boundary = 'bottom top'
    []
    [inflow]
        type = ADVectorFunctionDirichletBC
        variable = u
        boundary = 'right'
        function_x = '-sin(pi * y)'
    []
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
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
[]

[Outputs]
    exodus = true
[]
