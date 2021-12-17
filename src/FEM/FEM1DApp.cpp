#include <FEM/FEM1DApp.hpp>

Float StaticFEM1DApp::GradientSelfInnerProduct(int i, int j)
{
	std::vector<int> i_id, j_id;
	auto i_mesh = IdxToMesh(i, i_id);
	auto j_mesh = IdxToMesh(j, j_id);

	Float ret = 0;

#pragma omp parallel for
	for (int a = 0; a < i_mesh.size(); ++a)
	{
		for (int b = 0; b < j_mesh.size(); ++b)
		{
			if (i_mesh[a] == j_mesh[b])
			{
				auto sub_interval = interval.SubInterval(i_mesh[a]);
				ret += WeightedL2InnerProduct(sub_interval.remap(ShapeFunctions[i_id[a]]),
					sub_interval.remap(ShapeFunctionGradients[j_id[b]],
						1.0 / sub_interval.length()), b_func, sub_interval);
			}
		}
	}

	return ret;
}

Float StaticFEM1DApp::GradientInnerProduct(int i, int j)
{
	std::vector<int> i_id, j_id;
	auto i_mesh = IdxToMesh(i, i_id);
	auto j_mesh = IdxToMesh(j, j_id);

	Float ret = 0;

#pragma omp parallel for
	for (int a = 0; a < i_mesh.size(); ++a)
	{
		for (int b = 0; b < j_mesh.size(); ++b)
		{
			if (i_mesh[a] == j_mesh[b])
			{
				auto sub_interval = interval.SubInterval(i_mesh[a]);
				ret += WeightedL2InnerProduct(
					sub_interval.remap(ShapeFunctionGradients[i_id[a]], 1.0 / sub_interval.length()),
					sub_interval.remap(ShapeFunctionGradients[j_id[b]], 1.0 / sub_interval.length()), d_func,
					sub_interval);
			}
		}
	}

	return ret;
}

Float StaticFEM1DApp::SelfInnerProduct(int i, int j)
{
	std::vector<int> i_id, j_id;
	auto i_mesh = IdxToMesh(i, i_id);
	auto j_mesh = IdxToMesh(j, j_id);

	Float ret = 0;
#pragma omp parallel for
	for (int a = 0; a < i_mesh.size(); ++a)
	{
		for (int b = 0; b < j_mesh.size(); ++b)
		{
			if (i_mesh[a] == j_mesh[b])
			{
				auto sub_interval = interval.SubInterval(i_mesh[a]);
				ret += WeightedL2InnerProduct(sub_interval.remap(ShapeFunctions[i_id[a]]),
					sub_interval.remap(ShapeFunctions[j_id[b]]), c_func, sub_interval);
			}
		}
	}

	return ret;
}

Float StaticFEM1DApp::RHSInnerProduct(int i)
{
	std::vector<int> func_id;
	auto i_mesh = IdxToMesh(i, func_id);

	Float ret = 0;

	for (int a = 0; a < func_id.size(); ++a)
	{
		auto sub_interval = interval.SubInterval(i_mesh[a]);
		ret += L2InnerProduct(sub_interval.remap(ShapeFunctions[func_id[a]]), RHS_func, sub_interval);
	}

	return ret;
}

std::vector<int> StaticFEM1DApp::RelatedFuncIdx(int idx)
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

Float StaticFEM1DApp::Value(Float x)
{
	if (mat_size == 0)
	{
		return 0;
	}
	Float ret = 0;
	for (int i = 0; i < interval.GetPartitionCount(); ++i)
	{
		auto sub_interval = interval.SubInterval(i);
		if (sub_interval.Inside(x))
		{
			for (int j = 0; j < ShapeFunctions.size(); ++j)
			{
				int idx;
				if (MeshToIdx(i, j, idx))
				{
					ret += sub_interval.remap(ShapeFunctions[j])(x) * rst(idx);
				}
			}
		}
	}
	return ret;
}

std::function<Float(Float)> LagrangianBase(int N, int i)
{
	std::vector<Point2d> points(N + 1);
	Float h = 1.0 / N;
	for (int i = 0; i <= N; ++i)
	{
		points[i] = Point2d(i * h, 0);
	}
	points[i] = Point2d(i * h, 1.0);
	return LagrangianPolynomial(points);
}

std::function<Float(Float)> LagrangianBaseDerivative(int N, int i)
{
	return [=](Float x)
	{
		Float ret = 0;
		for (int missing = 0; missing <= N; ++missing)
		{
			if (missing != i)
			{
				std::vector<Point2d> points;
				Float h = 1.0 / N;

				for (int j = 0; j <= N; ++j)
					if (j != missing)
						if (j == i)
							points.emplace_back(j * h, 1.0);
						else
							points.emplace_back(j * h, 0.0);

				ret += LagrangianPolynomial(points)(x) / (h * (i - missing));
			}
		}
		return ret;
	};
}