[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = generated_mesh
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '1 1 1'
  []
  [./boundary]
    type = SideSetsBetweenSubdomainsGenerator
    input = new_block
    primary_block = '0'
    paired_block = 1
    new_boundary = 'interface'
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
                          'cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy '
                          'mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[Materials]
  [./stress]
    type = CauchyStressFromNEML
    database = "../../test_materials.xml"
    model = "elastic_model"
  [../]
  [./vonmises]
    type = EffectiveStressMaterial
    effective_stress_mp_name = vonmises
    effective_stress_type = VonMises
  [../]
  [./s1]
    type = EffectiveStressMaterial
    effective_stress_mp_name = s1
    effective_stress_type = MaxPrincipal
  [../]
  [./tresca]
    type = EffectiveStressMaterial
    effective_stress_mp_name = tresca
    effective_stress_type = Tresca
  [../]
  [./hydrostatic]
    type = EffectiveStressMaterial
    effective_stress_mp_name = hydrostatic
    effective_stress_type = Hydrostatic
  [../]
  [./huddleston]
    type = EffectiveStressMaterial
    effective_stress_mp_name = huddleston
    params_vector = '0.03'
    effective_stress_type = Huddleston
  [../]
  [./hayhurst]
    type = EffectiveStressMaterial
    effective_stress_mp_name = hayhurst
    params_vector = '0.2 0.3'
    effective_stress_type = Hayhurst
  [../]
  [./rccmrx_mises]
    type = EffectiveStressMaterial
    effective_stress_mp_name = rccmrx_mises
    params_vector = '0.2'
    effective_stress_type = RCCMRXMises
  [../]
  [./rccmrx_tresca]
    type = EffectiveStressMaterial
    effective_stress_mp_name = rccmrx_tresca
    params_vector = '0.2'
    effective_stress_type = RCCMRXTresca
  [../]
  [./max_s1_mises]
    type = EffectiveStressMaterial
    effective_stress_mp_name = max_s1_mises
    effective_stress_type = maxS1AndMises
  [../]
[]

[Functions]
 [./topfunc_x]
   type = PiecewiseLinear
   x = '0 1'
   y = '0 10'
 [../]
 [./topfunc_y]
   type = PiecewiseLinear
   x = '0 1'
   y = '0 20'
 [../]
 [./topfunc_z]
   type = PiecewiseLinear
   x = '0 1'
   y = '0 -30'
 [../]
[]

[BCs]
  [./x_0]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  [../]
  [./y_0]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  [./z_0]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  [../]
  [./x_1]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = 'right'
    function = topfunc_x
  [../]
  [./y_1]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = 'top'
    function = topfunc_y
  [../]
  [./z_1]
    type = FunctionNeumannBC
    variable = disp_z
    boundary = 'front'
    function = topfunc_z
  [../]

[]

[AuxVariables]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tresca]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hydrostatic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./huddleston]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hayhurst]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rccmrx_mises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rccmrx_tresca]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./max_s1_mises]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./vonmises]
    type = MaterialRealAux
    variable = vonmises
    property = vonmises
  [../]
  [./s1]
    type = MaterialRealAux
    variable = s1
    property = s1
  [../]
  [./tresca]
    type = MaterialRealAux
    variable = tresca
    property = tresca
  [../]
  [./hydrostatic]
    type = MaterialRealAux
    variable = hydrostatic
    property = hydrostatic
  [../]
  [./huddleston]
    type = MaterialRealAux
    variable = huddleston
    property = huddleston
  [../]
  [./hayhurst]
    type = MaterialRealAux
    variable = hayhurst
    property = hayhurst
  [../]
  [./rccmrx_mises]
    type = MaterialRealAux
    variable = rccmrx_mises
    property = rccmrx_mises
  [../]
  [./rccmrx_tresca]
    type = MaterialRealAux
    variable = rccmrx_tresca
    property = rccmrx_tresca
  [../]
  [./max_s1_mises]
    type = MaterialRealAux
    variable = max_s1_mises
    property = max_s1_mises
  [../]
[]

[Postprocessors]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./s1]
    type = ElementAverageValue
    variable = s1
  [../]
  [./tresca]
    type = ElementAverageValue
    variable = tresca
  [../]
  [./hydrostatic]
    type = ElementAverageValue
    variable = hydrostatic
  [../]
  [./huddleston]
    type = ElementAverageValue
    variable = huddleston
  [../]
  [./hayhurst]
    type = ElementAverageValue
    variable = hayhurst
  [../]
  [./rccmrx_mises]
    type = ElementAverageValue
    variable = rccmrx_mises
  [../]
  [./rccmrx_tresca]
    type = ElementAverageValue
    variable = rccmrx_tresca
  [../]
  [./max_s1_mises]
    type = ElementAverageValue
    variable = max_s1_mises
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

  solve_type = NEWTON
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 4
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dtmin = 1
  dt = 1
  end_time = 1.0
[]
[Outputs]
  exodus = true
[]
