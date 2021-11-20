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