[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = SMALL
        add_variables = true
        new_system = true
        formulation = UPDATED
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy '
                          'cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx '
                          'mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy '
                          'mechanical_strain_xz mechanical_strain_yz'
        eigenstrain_names = 'thermal'
      []
    []
  []
[]

[BCs]
  [x]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [right]
    type = DirichletBC
    preset = true
    boundary = right
    variable = disp_x
    value = 0.0
  []
  [y]
    type = DirichletBC
    preset = true
    boundary = bottom
    variable = disp_y
    value = 0.0
  []
  [z]
    type = DirichletBC
    preset = true
    boundary = back
    variable = disp_z
    value = 0.0
  []
[]

[AuxVariables]
  [T]
    initial_condition = 0
  []
  [T0]
    initial_condition = 0
  []
[]

[AuxKernels]
  [temp]
    type = FunctionAux
    variable = T
    function = temperature
  []
[]

[Functions]
  [temperature]
    type = ParsedFunction
    value = 100.0*t
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "test.xml"
    model = "elastic_model"
  []
  [eigenstrain_bases]
    type = GenericConstantVectorMaterial
    prop_names = 'n0 n1 n2'
    prop_values = '1 0 0 0 1 0 0 0 1'
  []
  [thermal_strain]
    type = ComputeAnisotropicThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = T0
    eigenstrain_name = thermal
    basis_0 = n0
    basis_1 = n1
    basis_2 = n2
    CTE_0 = 1e-3
    CTE_1 = 2e-3
    CTE_2 = 3e-3
  []
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
  nl_abs_tol = 1e-10

  dt = 1
  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
