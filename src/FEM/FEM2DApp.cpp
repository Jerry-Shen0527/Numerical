#include<FEM/FEM2DApp.hpp>

void StaticFEM2D::FillMatrix()
{
	std::vector<Eigen::Triplet<Float>> triplets;
	for (int i = 0; i < mat_size; ++i)
	{
		auto related_vec = RelatedFuncIdx(i);

		for (int j : related_vec)
		{
			triplets.emplace_back(i, j, GradientInnerProduct(i, j));
			triplets.emplace_back(i, j, SelfInnerProduct(i, j));
			triplets.emplace_back(i, j, GradientSelfInnerProduct(i, j));
		}
	}
	matrix_ = Eigen::SparseMatrix<Float>(mat_size, mat_size);
	matrix_.setFromTriplets(triplets.begin(), triplets.end());
}

void StaticFEM2D::FillRhs()
{
	rhs = Vector(mat_size);
	for (int i = 0; i < mat_size; ++i)
	{
		rhs(i) = RHSInnerProduct(i);
	}
}

std::vector<int> StaticFEM2DAppP1::IdxToMesh(int idx, std::vector<int>& shapeFuncId)
{
	std::vector<int> ret;
	auto v = mesh_.vertex_handle(idx);

	auto vf = mesh_.vf_ccw_range(v);

	for (auto smart_face_handle : vf)
	{
		ret.push_back(smart_face_handle.idx());

		auto vertices = smart_face_handle.vertices();
		auto indices = vertices.to_vector([](OpenMesh::SmartVertexHandle handle) {return handle.idx(); });

		std::sort(indices.begin(), indices.end());

		shapeFuncId.push_back(std::distance(indices.begin(), std::find(indices.begin(), indices.end(), idx)));
	}

	return ret;
}

bool StaticFEM2DAppP1::MeshToIdx(int mesh_idx, int shapefun_idx, int& idx)
{
	auto face_handle = mesh_.face_handle(mesh_idx);

	auto face_vertices = mesh_.fv_ccw_range(face_handle);

	auto indices = face_vertices.to_vector([](OpenMesh::SmartVertexHandle handle) {return handle.idx(); });

	std::sort(indices.begin(), indices.end());
	int non_zero_id = indices[shapefun_idx];

	for (auto face_vertex : face_vertices)
	{
		if (face_vertex.idx() == non_zero_id && face_vertex.is_boundary())
		{
			return false;
		}
	}

	return true;
}