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
[]

[Functions]
 [./topfunc]
   type = ParsedFunction
   value = 't'
 [../]
[]

[BCs]
  [./bottom3]
    type = PresetBC
    variable = disp_z
    boundary = 0
    value = 0.0
  [../]
  [./top]
    type = FunctionPresetBC
    variable = disp_z
    boundary = 5
    function = topfunc
  [../]
  [./corner1]
    type = PresetBC
    variable = disp_x
    boundary = 6
    value = 0.0
  [../]
  [./corner2]
    type = PresetBC
    variable = disp_y
    boundary = 6
    value = 0.0
  [../]
  [./corner3]
    type = PresetBC
    variable = disp_z
    boundary = 6
    value = 0.0
  [../]
  [./side1]
    type = PresetBC
    variable = disp_y
    boundary = 7
    value = 0.0
  [../]
  [./side2]
    type = PresetBC
    variable = disp_z
    boundary = 7
    value = 0.0
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
    type = EffectiveStressAux
    variable = vonmises
    effective_stress_type = VonMises
  [../]
  [./s1]
    type = EffectiveStressAux
    variable = s1
    effective_stress_type = MaxPrincipal
  [../]
  [./tresca]
    type = EffectiveStressAux
    variable = tresca
    effective_stress_type = Tresca
  [../]
  [./hydrostatic]
    type = EffectiveStressAux
    variable = hydrostatic
    effective_stress_type = Hydrostatic
  [../]
  [./huddleston]
    type = EffectiveStressAux
    variable = huddleston
    params_vector = '0.03'
    effective_stress_type = Huddleston
  [../]
  [./hayhurst]
    type = EffectiveStressAux
    variable = hayhurst
    params_vector = '0.2 0.3'
    effective_stress_type = Hayhurst
  [../]
  [./rccmrx_mises]
    type = EffectiveStressAux
    variable = rccmrx_mises
    params_vector = '0.2'
    effective_stress_type = RCCMRXMises
  [../]
  [./rccmrx_tresca]
    type = EffectiveStressAux
    variable = rccmrx_tresca
    params_vector = '0.2'
    effective_stress_type = RCCMRXTresca
  [../]
  [./max_s1_mises]
    type = EffectiveStressAux
    variable = max_s1_mises
    effective_stress_type = maxS1AndMises
  [../]
[]

[Postprocessors]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./s1]
    type = ElementAverageValue
    variable = s1
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./tresca]
    type = ElementAverageValue
    variable = tresca
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./hydrostatic]
    type = ElementAverageValue
    variable = hydrostatic
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./huddleston]
    type = ElementAverageValue
    variable = huddleston
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./hayhurst]
    type = ElementAverageValue
    variable = hayhurst
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./rccmrx_mises]
    type = ElementAverageValue
    variable = rccmrx_mises
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./rccmrx_tresca]
    type = ElementAverageValue
    variable = rccmrx_tresca
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./max_s1_mises]
    type = ElementAverageValue
    variable = max_s1_mises
    block = 'ANY_BLOCK_ID 0'
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
