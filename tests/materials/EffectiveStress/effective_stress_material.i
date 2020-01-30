[Mesh]
  displacements = 'disp_x disp_y disp_z'
  [generated_mesh]
    type = GeneratedMeshGenerator
    elem_type = HEX8
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
  []
  [node]
    type = ExtraNodesetGenerator
    coord = '0.0 0.0 0.0'
    new_boundary = 6
    input = generated_mesh
  []
  [snode]
    type = ExtraNodesetGenerator
    coord = '1.0 0.0 0.0'
    new_boundary = 7
    input = node
  []
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
 [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]



[Materials]
  [./fplastic]
    type = FiniteStrainPlasticMaterial
    block = 0
    yield_stress='0. 445. 0.05 610. 0.1 680. 0.38 810. 0.95 920. 2. 950.'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    C_ijkl = '2.827e5 1.21e5 1.21e5 2.827e5 1.21e5 2.827e5 0.808e5 0.808e5 0.808e5'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
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
 [./topfunc]
   type = ParsedFunction
   value = 't'
 [../]
[]

[BCs]
  [./bottom3]
    type = DirichletBC
    variable = disp_z
    boundary = 0
    value = 0.0
    preset = true
  [../]
  [./top]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 5
    function = topfunc
    preset = true
  [../]
  [./corner1]
    type = DirichletBC
    variable = disp_x
    boundary = 6
    value = 0.0
    preset = true
  [../]
  [./corner2]
    type = DirichletBC
    variable = disp_y
    boundary = 6
    value = 0.0
    preset = true
  [../]
  [./corner3]
    type = DirichletBC
    variable = disp_z
    boundary = 6
    value = 0.0
    preset = true
  [../]
  [./side1]
    type = DirichletBC
    variable = disp_y
    boundary = 7
    value = 0.0
    preset = true
  [../]
  [./side2]
    type = DirichletBC
    variable = disp_z
    boundary = 7
    value = 0.0
    preset = true
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

[Executioner]

  type = Transient

  dt=0.1
  dtmin=0.1
  dtmax=1
  end_time=1.0

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
