#pragma once

#include "MeshGenerator.h"
#include "libmesh/elem.h"

/*
 * A mesh generator to split a mesh by breaking all element-element interfaces in the
 * specified subdomains
 */
class ExplodeMeshGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  ExplodeMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  std::unordered_map<dof_id_type, std::vector<dof_id_type>>
  buildSubdomainRestrictedNodeToElemMap(std::unique_ptr<MeshBase> & mesh,
                                        const std::vector<SubdomainID> & subdomains) const;

  void duplicateNodes(
      std::unique_ptr<MeshBase> & mesh,
      const std::unordered_map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map) const;

  void duplicateNode(std::unique_ptr<MeshBase> & mesh, Elem * elem, const Node * node) const;

  void createInterface(
      MeshBase & mesh,
      const std::unordered_map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map) const;

  /// The mesh to modify
  std::unique_ptr<MeshBase> & _input;

  // The subdomains to explode
  const std::vector<SubdomainID> & _subdomains;

  // The name of the new boundary
  const BoundaryName _interface_name;
};
