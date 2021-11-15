#pragma once
#include"FEM.hpp"

#include "Geometry/Mesh/HEMesh.hpp"

class StaticFEM2D :public StaticFEM
{
protected:
	void SetMatSize() override;
	void FillMatrix() override;
	void FillRhs() override;
	Float GradientInnerProduct(int i, int j) override;
	Float SelfInnerProduct(int i, int j) override;
	Float GradientSelfInnerProduct(int i, int j) override;
	Float RHSInnerProduct(int i) override;
	std::vector<int> RelatedFuncIdx(int idx) override;
};

class FEM2DApp :public StaticFEM2D
{

};