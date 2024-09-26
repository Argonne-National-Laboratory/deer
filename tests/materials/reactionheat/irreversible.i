[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        new_system = true
        formulation = UPDATED
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[Functions]
  [temperature]
    type = PiecewiseLinear
    x = '0 50 200'
    y = '0 51 0'
  []
[]

[AuxVariables]
  [temperature]
  []
  [reaction_heat]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [temperature]
    type = FunctionAux
    function = temperature
    variable = temperature
  []
  [reaction_heat]
    type = MaterialRealAux
    variable = reaction_heat
    property = "reaction_power"
  []
[]

[BCs]
  [leftx]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [lefty]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_y
    value = 0.0
  []
  [leftz]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_z
    value = 0.0
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "../../test_materials.xml"
    model = "elastic_model"
    large_kinematics = false
  []
  [reaction_heat]
    type = EmpiricalReactionHeat
    total_heat = 100.0
    start_temperature = 50.0
    time_constant = 10.0
    temperature = temperature
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
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

  start_time = 0.0
  dt = 10.0
  end_time = 200.0
[]

[Postprocessors]
  [temperature]
    type = ElementAverageValue
    variable = temperature
  []
  [power]
    type = ElementAverageValue
    variable = reaction_heat
  []
  [energy]
    type = TimeIntegratedPostprocessor
    value = power
  []
[]

[Outputs]
  csv = true
[]
