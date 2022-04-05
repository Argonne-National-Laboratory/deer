[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 20
    zmax = 10
  []
[]

[Variables]
  [./temperature]
  [../]
[]

[Kernels]
  [./hc]
    type = HeatConduction
    variable = temperature
  [../]
  [./ht]
    type = HeatCapacityConductionTimeDerivative
    variable = temperature
  [../]
  [./hs]
    type = AdiabaticHeating
    variable = temperature
    heat = reaction_heat
  [../]
[]

[BCs]
  [heat_up]
    type = NeumannBC
    boundary = 'front'
    value = 100.0
    variable = temperature
  [../]
[]

[Materials]
  [./thermal]
    type = GenericConstantMaterial
    prop_names = 'density heat_capacity thermal_conductivity'
    prop_values = '1.0 12.0 4.0'
  [../]
  [./reaction_heat]
    type = EmpiricalReactionHeat
    total_heat = 10000.0
    start_temperature = 100.0
    time_constant = 200.0
    temperature = temperature
  [../]
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[AuxVariables]
  [./reaction_heat]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./reaction_heat]
    type = MaterialRealAux
    variable = reaction_heat
    property = "reaction_power"
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-12
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  dt = 10.0
  end_time = 200.0
[]

[Outputs]
  exodus = true
[]
