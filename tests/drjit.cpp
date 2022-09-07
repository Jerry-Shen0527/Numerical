#include <iostream>
#include<drjit/array.h>

#include "drjit/jit.h"

int main()
{
	using namespace drjit;

	jit_init();

	LLVMArray<float> arr = LLVMArray<float>(3);

	std::cout << arr << std::endl;

}
