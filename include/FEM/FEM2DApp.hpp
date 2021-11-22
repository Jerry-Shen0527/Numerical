#pragma once
#include <set>

#include"FEM.hpp"
#define _USE_MATH_DEFINES
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include "Domain.hpp"
#include "NIntegrate/2DIntegrate.hpp"
#include "NIntegrate/2DIntegrate.hpp"
#include "NIntegrate/Integrate.hpp"

class StaticFEM2D :public StaticFEM
{
protected:
	void FillMatrix() override;
	void FillRhs() override;

public:
	virtual Float Value(int i, Eigen::Vector3d vector) = 0;
};

class StaticFEM2DApp :public StaticFEM2D
{
public:
	using FEM2DMesh = OpenMesh::PolyMesh_ArrayKernelT<>;
	TriangleDomain BuildTriangleElement(int face_idx)
	{
		auto mesh_triangle = mesh_.ArrayKernel::face_handle(face_idx);
		auto face_vertices = mesh_.fv_ccw_range(mesh_triangle);
		auto indices = face_vertices.to_vector([](OpenMesh::SmartVertexHandle handle) {return handle.idx(); });

		std::sort(indices.begin(), indices.end());

		Eigen::Vector3d ps[3];
		for (int i = 0; i < 3; ++i)
		{
			auto point = mesh_.point(mesh_.vertex_handle(indices[i]));
			for (int j = 0; j < 3; ++j)
			{
				ps[i][j] = point[j];
			}
		}

		TriangleDomain triangle(ps[0], ps[1], ps[2]);

		return triangle;
	}

	StaticFEM2DApp(const std::function<Float(Eigen::Vector3d)>& rhs_func, const std::function<Float(Eigen::Vector3d)>& d_func, const std::function<Float(Eigen::Vector3d)>& c_func) : RHS_func(rhs_func), d_func(d_func),
		c_func(c_func) {}

	Float Value(int i, Eigen::Vector3d x) override
	{
		if (mat_size == 0)
		{
			return 0;
		}
		Float ret = 0;

		auto triangle = BuildTriangleElement(i);

		for (int j = 0; j < ShapeFunctions.size(); ++j)
		{
			int idx;
			MeshToIdx(i, j, idx);

			ret += triangle.remap(ShapeFunctions[j])(x) * rst(idx);
		}

		return ret;
	}

	FEM2DMesh& GetMesh() { return mesh_; }

protected:

	Float GradientInnerProduct(int i, int j) override
	{
		std::vector<int> i_id, j_id;
		auto i_mesh = IdxToMesh(i, i_id);
		auto j_mesh = IdxToMesh(j, j_id);

		Float ret = 0;

		for (int a = 0; a < i_mesh.size(); ++a)
		{
			for (int b = 0; b < j_mesh.size(); ++b)
			{
				if (i_mesh[a] == j_mesh[b])
				{
					int face_idx = i_mesh[a];

					auto triangle = BuildTriangleElement(face_idx);

					auto shape_func1 = triangle.remap(ShapeFunctions[i_id[a]]);
					auto shape_func2 = triangle.remap(ShapeFunctions[j_id[b]]);

					NDifferential<3> gradient;

					ret += WeightedL2InnerProduct2D<3>(gradient(shape_func1), gradient(shape_func2), d_func, triangle);
				}
			}
		}

		return ret;
	}

	Float SelfInnerProduct(int i, int j) override
	{
		std::vector<int> i_id, j_id;
		auto i_mesh = IdxToMesh(i, i_id);
		auto j_mesh = IdxToMesh(j, j_id);

		Float ret = 0;

		for (int a = 0; a < i_mesh.size(); ++a)
		{
			for (int b = 0; b < j_mesh.size(); ++b)
			{
				if (i_mesh[a] == j_mesh[b])
				{
					int face_idx = i_mesh[a];
					auto triangle = BuildTriangleElement(face_idx);

					auto shape_func1 = triangle.remap(ShapeFunctions[i_id[a]]);
					auto shape_func2 = triangle.remap(ShapeFunctions[j_id[b]]);

					ret += WeightedL2InnerProduct2D(shape_func1, shape_func2, c_func, triangle);
				}
			}
		}

		return ret;
	}

	Float GradientSelfInnerProduct(int i, int j) override
	{
		return 0;
	}

	Float RHSInnerProduct(int i) override
	{
		std::vector<int> func_id;
		auto i_mesh = IdxToMesh(i, func_id);

		Float ret = 0;

		for (int a = 0; a < func_id.size(); ++a)
		{
			auto triangle = BuildTriangleElement(i_mesh[a]);
			ret += L2InnerProduct2D(triangle.remap(ShapeFunctions[func_id[a]]), RHS_func, triangle);
		}

		return ret;
	}

	std::vector<int> RelatedFuncIdx(int idx) override
	{
		std::vector<int> ret;
		std::vector<int> foo_id;

		auto MeshIds = IdxToMesh(idx, foo_id);

		std::set<int> set_ret;

		for (auto mesh_id : MeshIds)
		{
			for (int i = 0; i < ShapeFunctions.size(); ++i)
			{
				int idx;
				if (MeshToIdx(mesh_id, i, idx))
				{
					set_ret.emplace(idx);
				}
			}
		}
		ret.assign(set_ret.begin(), set_ret.end());
		return ret;
	}

	virtual std::vector<int> IdxToMesh(int idx, std::vector<int>& shapeFuncId) = 0;
	virtual bool MeshToIdx(int mesh_idx, int shapefun_idx, int& idx) = 0;

	std::vector<std::function<Float(Eigen::Vector2d)>> ShapeFunctions;

protected:
	FEM2DMesh mesh_;

	std::function<Float(Eigen::Vector3d)> RHS_func;
	std::function<Float(Eigen::Vector3d)> d_func;
	std::function<Float(Eigen::Vector3d)> c_func;
};

class StaticFEM2DAppP1 :public StaticFEM2DApp
{
public:
	StaticFEM2DAppP1(const std::function<Float(Eigen::Vector3d)>& rhs_func, const std::function<Float(Eigen::Vector3d)>& d_func, const std::function<Float(Eigen::Vector3d)>& c_func) : StaticFEM2DApp(rhs_func, d_func, c_func)
	{
		TriangleDomain domain;
		ShapeFunctions =
		{
			domain.scale
			(
				[](Eigen::Vector3d vec)
				{
					assert(vec.norm() < 1E10);
					return 1 - vec.x() - vec.y();
				}
			),
			domain.scale([](Eigen::Vector3d vec) { return vec.x(); }),
			domain.scale([](Eigen::Vector3d vec) { return vec.y(); }),
		};
	}

protected:
	std::vector<int> IdxToMesh(int idx, std::vector<int>& shapeFuncId) override;
	bool MeshToIdx(int mesh_idx, int shapefun_idx, int& idx) override;

	void SetMatSize() override
	{
		int total = mesh_.n_vertices();

		mat_size = total;
	}
};

class StaticFEM2DAppP2 :public StaticFEM2DApp
{
public:
	StaticFEM2DAppP2(const std::function<Float(Eigen::Vector3d)>& rhs_func, const std::function<Float(Eigen::Vector3d)>& d_func, const std::function<Float(Eigen::Vector3d)>& c_func) : StaticFEM2DApp(rhs_func, d_func, c_func)
	{
		TriangleDomain domain;
		ShapeFunctions =
		{
			domain.scale
			(
				[](Eigen::Vector3d vec)
				{
					assert(vec.norm() < 1E10);
					return 2 * (1 - vec.x() - vec.y()) * (1.0 / 2.0 - vec.x() - vec.y());
				}
			),
			domain.scale([](Eigen::Vector3d vec) { return 2.0 * vec.x() * (vec.x() - 1.0 / 2); }),
			domain.scale([](Eigen::Vector3d vec) { return 2.0 * vec.y() * (vec.y() - 1.0 / 2); }),

			domain.scale
			(
				[](Eigen::Vector3d vec)
				{
					assert(vec.norm() < 1E10);
					return 4.0 * vec.x() * vec.y();
				}
			),
			domain.scale([](Eigen::Vector3d vec) { return 4.0 * vec.x() * (1 - vec.x() - vec.y()); }),
			domain.scale([](Eigen::Vector3d vec) { return 4.0 * (1 - vec.x() - vec.y()) * vec.y(); }),
		};
	}

protected:
	std::vector<int> IdxToMesh(int idx, std::vector<int>& shapeFuncId) override
	{
		std::vector<int> ret;
		if (idx < mesh_.n_vertices())
		{
			auto v = mesh_.vertex_handle(idx);
			if (mesh_.is_boundary(v))
			{
				return std::vector<int>();
			}

			auto vf = mesh_.vf_range(v);

			for (auto smart_face_handle : vf)
			{
				ret.push_back(smart_face_handle.idx());

				auto vertices = smart_face_handle.vertices();
				auto indices = vertices.to_vector([](OpenMesh::SmartVertexHandle handle) {return handle.idx(); });

				std::sort(indices.begin(), indices.end());

				auto dis = std::distance(indices.begin(), std::find(indices.begin(), indices.end(), idx));

				shapeFuncId.push_back(dis);
			}
		}
		else
		{
			idx = idx - mesh_.n_vertices();
			auto e = mesh_.edge_handle(idx);
			if (mesh_.is_boundary(e))
			{
				return std::vector<int>();
			}

			std::vector<OpenMesh::HalfedgeHandle> he_s(2);
			he_s[0] = mesh_.halfedge_handle(e, 0);
			he_s[1] = mesh_.halfedge_handle(e, 1);

			//auto vf = mesh_.vf_range(v);
			std::vector<OpenMesh::FaceHandle> handles(2);

			handles[0] = mesh_.face_handle(he_s[0]);
			handles[1] = mesh_.face_handle(he_s[1]);

			for (auto face_handle : handles)
			{
				ret.push_back(face_handle.idx());

				auto edges = mesh_.fe_range(face_handle);

				auto edge_indices = edges.to_vector([&](OpenMesh::SmartEdgeHandle handle)
					{
						auto fv = mesh_.fv_range(face_handle);
						int v_idx = 0;
						for (auto smart_vertex_handle : fv)
						{
							if (handle.v0() != smart_vertex_handle && handle.v1() != smart_vertex_handle)
							{
								v_idx = smart_vertex_handle.idx();
								break;
							}
						}

						return std::make_pair(handle.idx(), v_idx);
					});

				std::sort(edge_indices.begin(), edge_indices.end(), [](std::pair<int, int>pair, std::pair<int, int>pair2) {return pair.second < pair2.second; });

				auto dis = std::distance(edge_indices.begin(), std::find_if(edge_indices.begin(), edge_indices.end(), [idx](std::pair<int, int>pair) {return idx == pair.first; }));

				shapeFuncId.push_back(3 + dis);
			}
		}

		return ret;
	}

	bool MeshToIdx(int mesh_idx, int shapefun_idx, int& idx) override
	{
		auto face_handle = mesh_.face_handle(mesh_idx);

		if (shapefun_idx < 3)
		{
			auto face_vertices = mesh_.fv_ccw_range(face_handle);

			auto indices = face_vertices.to_vector([](OpenMesh::SmartVertexHandle handle) {return handle.idx(); });

			std::sort(indices.begin(), indices.end());
			int non_zero_id = indices[shapefun_idx];

			for (auto face_vertex : face_vertices)
			{
				if (face_vertex.idx() == non_zero_id)
				{
					idx = face_vertex.idx();
				}
			}
		}
		else
		{
			auto face_vertices = mesh_.fe_range(face_handle);

			auto indices = face_vertices.to_vector([](OpenMesh::SmartEdgeHandle handle) {return handle.idx(); });

			std::sort(indices.begin(), indices.end());
			int non_zero_id = indices[shapefun_idx - 3];

			for (auto face_edge : face_vertices)
			{
				if (face_edge.idx() == non_zero_id)
				{
					idx = face_edge.idx();
				}
			}
		}

		return true;
	}

	void SetMatSize() override
	{
		int total = mesh_.n_vertices() + mesh_.n_edges();

		mat_size = total;
	}
};