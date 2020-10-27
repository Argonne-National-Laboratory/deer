[Mesh]
  ## This file is where you split your mesh
  ## execution must NOT be --mes-only as we are saving element integers for later use
  [msh]
    type = FileMeshGenerator
    file = RVE_test_mesh.e
  []
  [breakmesh]
    input = msh
    type = BreakMeshByBlockGenerator
    # block = '4 5'
    # add_transition_interface = true
    write_fake_neighbor_list_to_file = true
    fake_neighbor_list_file_name = 'fake_neighbors_test_bmbb.csv'
  []
  [add_side_sets]
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

[AuxVariables]
  [bmbb_element_id]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [set_material_id]
    type = ElemExtraIDAux
    variable = bmbb_element_id
    extra_id_name = bmbb_element_id
  []
[]

[Problem]
  kernel_coverage_check = false
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
