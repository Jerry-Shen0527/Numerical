#pragma once
#include <set>

#include"FEM.hpp"
#define _USE_MATH_DEFINES
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

class StaticFEM2D :public StaticFEM
{
protected:
	void FillMatrix() override;
	void FillRhs() override;

public:
	virtual Float Value(int i, Eigen::Vector2f vector) = 0;
};

class StaticFEM2DApp :public StaticFEM2D
{
public:
	using FEM2DMesh = OpenMesh::PolyMesh_ArrayKernelT<>;

	StaticFEM2DApp(const FEM2DMesh& mesh) :mesh_(mesh) {}

	Float Value(int i, Eigen::Vector2f vector) override;
protected:
	void SetMatSize() override;
	Float GradientInnerProduct(int i, int j) override;
	Float SelfInnerProduct(int i, int j) override;
	Float GradientSelfInnerProduct(int i, int j) override;
	Float RHSInnerProduct(int i) override;

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

	std::vector<std::function<Float(Eigen::Vector2f)>> ShapeFunctions;
	std::vector<std::function<Eigen::Vector2f(Eigen::Vector2f)>> ShapeFunctionGradients;

protected:
	FEM2DMesh mesh_;
};

class StaticFEM2DAppPoly :public StaticFEM2DApp
{
	StaticFEM2DAppPoly(const FEM2DMesh& mesh) :StaticFEM2DApp(mesh) {}

protected:
	std::vector<int> IdxToMesh(int idx, std::vector<int>& shapeFuncId) override;
	bool MeshToIdx(int mesh_idx, int shapefun_idx, int& idx) override;
};