#include <iostream>
#include<drjit/array.h>

#include "drjit/jit.h"
#include <drjit/autodiff.h>

#include "drjit/tensor.h"
#include "drjit/math.h"
#include "drjit/random.h"

namespace dr = drjit;

template <typename Value_, size_t Size_>
struct Vector : dr::StaticArrayImpl<Value_, Size_, false, Vector<Value_, Size_>> {
    using Base = dr::StaticArrayImpl<Value_, Size_, false, Vector<Value_, Size_>>;

    /// Helper alias used to implement type promotion rules
    template <typename T> using ReplaceValue = Vector<T, Size_>;

    using ArrayType = Vector;
    using MaskType = dr::Mask<Value_, Size_>;

    DRJIT_ARRAY_IMPORT(Vector, Base)
};

using Float =  dr::DiffArray< dr::CUDAArray<float>>;

template <typename Value,
          typename T = std::conditional_t<dr::is_static_array_v<Value>,
                                          dr::value_t<Value>, Value>>
using DynamicBuffer =
std::conditional_t<dr::is_dynamic_array_v<T>, T, dr::DynamicArray<T>>;
//using Float = float;


using TensorXf = dr::Tensor<DynamicBuffer<Float>>;

size_t shape[] = { 20000,20000 };

int main()
{
	using namespace drjit;

	jit_init();
    dr::PCG32<TensorXf::Array> rand;

	TensorXf u(100, 2, shape);
    TensorXf v(50, 2, shape);
    auto& arr = v.array();
    arr.set_entry(0, 10);




    //

    enable_grad(u);
    enable_grad(v);


    auto loss =sum((u*v).array());
    backward(loss);
    //


	//for (int i = 0; i < 10; ++i)
	//{
 //       d += v;
 //       d *= v;
 //       d -= v;
	//}
    auto grad = u.array().grad_();
    std::cout << grad << std::endl;
    //auto a = u.;
    //auto b = v[0].grad_();

    //std::cout << u * v << std::endl;
    //std::cout << a << std::endl;
    //std::cout << b << std::endl;

    
	//auto s= std::string(Float::graphviz_());
	//std::cout << s << std::endl;



}
