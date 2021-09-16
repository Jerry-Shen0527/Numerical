#pragma once

#include "Visualizer.h"

class Visualizer2D :public Visualizer
{
public:
	Visualizer2D();
	~Visualizer2D();
	void Init() override;

private:

	void InitWindow();
};

Visualizer2D::Visualizer2D()
{
}

Visualizer2D::~Visualizer2D()
{
}

void Visualizer2D::Init()
{
}
