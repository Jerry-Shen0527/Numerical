//#include <iostream>
//#include <numeric>
//#include <type.hpp>
//
//#include <Eigen/Eigen>
//
//#include "FEM/FEM1DApp.hpp"
//#include "FEM/FEM2DApp.hpp"
//#include "imgui/implot.h"
//#include "Visualization/Visualizer.h"
//
//#include <numeric>
//
//using Linear = PolynomialFEMApp<1>;
//using Quadratic = PolynomialFEMApp<2>;
//
//class FEM2DVisualizer :public Visualizer
//{
//protected:
//
//	StaticFEM2DAppP1 app;
//protected:
//
//	void evaluate()
//	{
//		app.evaluate();
//	}
//
//	void Control_UI();
//	void draw(bool* p_open) override;
//
//	bool updated = true;
//
//	void error(std::vector<float>& ref, std::vector<float>& eval, Float& L_1, Float& L_2, Float& L_inf)
//	{
//		assert(ref.size() == eval.size());
//
//		std::vector<Float> minus(ref.size());
//
//		for (int i = 0; i < ref.size(); ++i)
//		{
//			minus[i] = abs(ref[i] - eval[i]);
//		}
//
//		L_inf = *std::max_element(minus.begin(), minus.end(), [](Float a, Float b) {return a < b; });
//		L_2 = sqrt(std::accumulate(minus.begin(), minus.end(), static_cast<Float>(0), [](Float r, Float a) {return r + a * a; }) / Float(ref.size()));
//		L_1 = std::accumulate(minus.begin(), minus.end(), static_cast<Float>(0), [](Float r, Float a) {return r + a; }) / Float(ref.size());
//	}
//
//public:
//	void CalcAccurateRst()
//	{
//	}
//
//	FEM2DVisualizer(const StaticFEM2DApp::FEM2DMesh& mesh) :app(mesh, [](Eigen::Vector3d vec) {return Float(0); }, [](Eigen::Vector3d vec) {return Float(1); }, [](Eigen::Vector3d vec) {return Float(0); }) {
//		CalcAccurateRst();
//		//segemnt = 16;
//		//Float L1, L2, L_inf;
//		//do
//		//{
//		//	evaluate();
//
//		//	error(precise_val, linear_val, L1, L2, L_inf);
//
//		//	using std::cout;
//		//	using std::endl;
//
//		//	//if (segemnt == 16)
//		//	//{
//		//	//	cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
//		//	//}
//		//	//else
//		//	//	cout << segemnt << '&' << L1 << '&' <<- log2(L1 / linear_L1.back()) << '&' << L2 << '&' << -log2(L2 / linear_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / linear_Linf.back()) << "\\\\" << endl;
//
//		//	pointcount.push_back(segemnt);
//		//	linear_L1.push_back(L1);
//		//	linear_L2.push_back(L2);
//		//	linear_Linf.push_back(L_inf);
//		//	error(precise_val, quadratic_val, L1, L2, L_inf);
//
//		//	if (segemnt == 16)
//		//	{
//		//		cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
//		//	}
//		//	else
//		//		cout << segemnt << '&' << L1 << '&' << -log2(L1 / quadratic_L1.back()) << '&' << L2 << '&' << -log2(L2 / quadratic_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / quadratic_Linf.back()) << "\\\\" << endl;
//
//		//	quadratic_L1.push_back(L1);
//		//	quadratic_L2.push_back(L2);
//		//	quadratic_Linf.push_back(L_inf);
//
//		//	segemnt *= 2;
//		//} while (segemnt != 4096);
//		//segemnt = 16;
//
//		evaluate();
//	}
//	const size_t Length = 10001;
//
//	std::vector<float> xs = std::vector<float>(Length);
//	std::vector<float> precise_val = std::vector<float>(Length);
//	std::vector<float> quadratic_val = std::vector<float>(Length);
//	std::vector<float> linear_val = std::vector<float>(Length);
//	std::vector<float> quadratic_diff = std::vector<float>(Length);
//	std::vector<float> linear_diff = std::vector<float>(Length);
//
//	std::vector<float> pointcount;
//	std::vector<float> linear_L1;
//	std::vector<float> linear_L2;
//	std::vector<float> linear_Linf;
//
//	std::vector<float> quadratic_L1;
//	std::vector<float> quadratic_L2;
//	std::vector<float> quadratic_Linf;
//};
//
//static inline ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs) { return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y); }
//
//void FEM2DVisualizer::Control_UI()
//{
//	//if (ImGui::SliderInt("Number of segments", &segemnt, 2, 200))
//	//{
//	//	segemnt = segemnt < 2 ? 2 : segemnt;
//	//	evaluate();
//	//}
//	evaluate();
//	CalcAccurateRst();
//}
//
//void FEM2DVisualizer::draw(bool* p_open)
//{
//}

/*
* Copyright (c) 2014-2021, NVIDIA CORPORATION. All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*/

#include "FEM/FEM2DApp.hpp"

#include <donut/app/ApplicationBase.h>
#include <donut/engine/ShaderFactory.h>
#include <donut/engine/TextureCache.h>
#include <donut/engine/CommonRenderPasses.h>
#include <donut/app/DeviceManager.h>
#include <donut/core/log.h>
#include <donut/core/vfs/VFS.h>
#include <nvrhi/utils.h>

#include "donut/app/Camera.h"
#include "donut/app/imgui_console.h"
#include "donut/app/imgui_renderer.h"
#include "donut/app/UserInterfaceUtils.h"
#include "donut/engine/ConsoleInterpreter.h"
#include "donut/engine/Scene.h"
#include "donut/engine/SceneGraph.h"
#include "donut/engine/View.h"
#include "donut/render/SkyPass.h"

#define _USE_MATH_DEFINES
#undef dim
#include <OpenMesh/Core/IO/MeshIO.hh>

using namespace donut;
using namespace donut::math;
using namespace donut::app;
using namespace donut::vfs;
using namespace donut::engine;

static const char* g_WindowTitle = "Donut Example: Vertex Buffer";

struct Vertex
{
	math::float3 position;
};

struct UIData
{
	bool                                ShowUI = true;
	bool                                ShowConsole = false;
	bool                                UseDeferredShading = true;
	bool                                Stereo = false;
	bool                                EnableSsao = true;

	render::SkyParameters                       SkyParams;

	bool                                EnableVsync = true;
	bool                                ShaderReoladRequested = false;
	bool                                EnableProceduralSky = true;
	bool                                EnableBloom = true;
	float                               BloomSigma = 32.f;
	float                               BloomAlpha = 0.05f;
	bool                                EnableTranslucency = true;
	bool                                EnableMaterialEvents = false;
	bool                                EnableShadows = true;
	float                               AmbientIntensity = 1.0f;
	bool                                EnableLightProbe = true;
	float                               LightProbeDiffuseScale = 1.f;
	float                               LightProbeSpecularScale = 1.f;
	float                               CsmExponent = 4.f;
	bool                                DisplayShadowMap = false;
	bool                                UseThirdPersonCamera = true;
	bool                                EnableAnimations = false;
	std::shared_ptr<Material>           SelectedMaterial;
	std::shared_ptr<SceneGraphNode>     SelectedNode;
	std::string                         ScreenshotFileName;
	std::shared_ptr<SceneCamera>        ActiveSceneCamera;
};

auto rhs_func = [](Eigen::Vector3d vector) {
	Float x = vector.x();
	Float y = vector.y();

	return  -10 * (-2 * cos(1 - y) * cos(y) * sin(1 - x) * sin(x) - 2 * cos(1 - x) * cos(x) * sin(1 - y) * sin(y) -
		4 * sin(1 - x) * sin(x) * sin(1 - y) * sin(y)); };

class FEM2DVisualizer : public app::IRenderPass
{
	std::vector<Vertex> g_Vertices;

	std::vector< uint32_t> g_Indices;

private:
	nvrhi::ShaderHandle m_VertexShader;
	nvrhi::ShaderHandle m_PixelShader;
	nvrhi::BufferHandle m_ConstantBuffer;
	nvrhi::BufferHandle m_VertexBuffer;
	nvrhi::BufferHandle m_IndexBuffer;
	//nvrhi::TextureHandle m_Texture;
	nvrhi::InputLayoutHandle m_InputLayout;
	nvrhi::BindingLayoutHandle m_BindingLayout;
	nvrhi::BindingSetHandle m_BindingSet;
	nvrhi::GraphicsPipelineHandle m_Pipeline;
	nvrhi::CommandListHandle m_CommandList;
	float m_Rotation = 0.f;

public:
	using IRenderPass::IRenderPass;

	FEM2DVisualizer(DeviceManager* deviceManager, UIData& ui) :IRenderPass(deviceManager), m_ui(ui), fem_app(rhs_func, [](Eigen::Vector3d vector) {return 1.0; }, [](Eigen::Vector3d vector) {return 0.0; })
	{
	}

	bool Init()
	{
		auto nativeFS = std::make_shared<vfs::NativeFileSystem>();

		std::filesystem::path frameworkShaderPath = app::GetDirectoryWithExecutable() / "shaders/framework" / app::GetShaderTypeName(GetDevice()->getGraphicsAPI());
		std::filesystem::path appShaderPath = app::GetDirectoryWithExecutable() / "shaders/vertex_buffer" / app::GetShaderTypeName(GetDevice()->getGraphicsAPI());

		m_RootFs = std::make_shared<vfs::RootFileSystem>();
		m_RootFs->mount("/shaders/donut", frameworkShaderPath);
		m_RootFs->mount("/shaders/app", appShaderPath);

		m_SceneFilesAvailable = { "gd0.obj","gd1.obj","gd2.obj","gd3.obj","gd4.obj" };
		SetCurrentSceneName("gd0.obj");

		m_ShaderFactory = std::make_shared<engine::ShaderFactory>(GetDevice(), m_RootFs, "/shaders");
		m_VertexShader = m_ShaderFactory->CreateShader("app/shaders.hlsl", "main_vs", nullptr, nvrhi::ShaderType::Vertex);
		m_PixelShader = m_ShaderFactory->CreateShader("app/shaders.hlsl", "main_ps", nullptr, nvrhi::ShaderType::Pixel);

		if (!m_VertexShader || !m_PixelShader)
		{
			return false;
		}

		m_ConstantBuffer = GetDevice()->createBuffer(nvrhi::utils::CreateVolatileConstantBufferDesc(sizeof(math::float4x4), "ConstantBuffer", engine::c_MaxRenderPassConstantBufferVersions));

		nvrhi::VertexAttributeDesc attributes[] = {
			nvrhi::VertexAttributeDesc()
				.setName("POSITION")
				.setFormat(nvrhi::Format::RGB32_FLOAT)
				.setOffset(offsetof(Vertex, position))
				.setElementStride(sizeof(Vertex)),
		};
		m_InputLayout = GetDevice()->createInputLayout(attributes, uint32_t(std::size(attributes)), m_VertexShader);

		engine::CommonRenderPasses commonPasses(GetDevice(), m_ShaderFactory);
		engine::TextureCache textureCache(GetDevice(), nativeFS, nullptr);

		m_CommandList = GetDevice()->createCommandList();
		m_CommandList->open();

		nvrhi::BufferDesc vertexBufferDesc;
		vertexBufferDesc.byteSize = sizeof(Vertex) * g_Vertices.size();
		vertexBufferDesc.isVertexBuffer = true;
		vertexBufferDesc.debugName = "VertexBuffer";
		vertexBufferDesc.initialState = nvrhi::ResourceStates::CopyDest;
		m_VertexBuffer = GetDevice()->createBuffer(vertexBufferDesc);

		m_CommandList->beginTrackingBufferState(m_VertexBuffer, nvrhi::ResourceStates::CopyDest);
		m_CommandList->writeBuffer(m_VertexBuffer, &g_Vertices[0], sizeof(Vertex) * g_Vertices.size());
		m_CommandList->setPermanentBufferState(m_VertexBuffer, nvrhi::ResourceStates::VertexBuffer);

		nvrhi::BufferDesc indexBufferDesc;
		indexBufferDesc.byteSize = g_Indices.size() * sizeof(uint32_t);
		indexBufferDesc.isIndexBuffer = true;
		indexBufferDesc.debugName = "IndexBuffer";
		indexBufferDesc.initialState = nvrhi::ResourceStates::CopyDest;
		m_IndexBuffer = GetDevice()->createBuffer(indexBufferDesc);

		m_CommandList->beginTrackingBufferState(m_IndexBuffer, nvrhi::ResourceStates::CopyDest);
		m_CommandList->writeBuffer(m_IndexBuffer, &g_Indices[0], g_Indices.size() * sizeof(uint32_t));
		m_CommandList->setPermanentBufferState(m_IndexBuffer, nvrhi::ResourceStates::IndexBuffer);

		m_CommandList->close();
		GetDevice()->executeCommandList(m_CommandList);

		nvrhi::BindingSetDesc bindingSetDesc;
		bindingSetDesc.bindings = {
			nvrhi::BindingSetItem::ConstantBuffer(0, m_ConstantBuffer),
		};

		if (!nvrhi::utils::CreateBindingSetAndLayout(GetDevice(), nvrhi::ShaderType::All, 0, bindingSetDesc, m_BindingLayout, m_BindingSet))
		{
			log::error("Couldn't create the binding set or layout");
			return false;
		}

		return true;
	}

	void Animate(float fElapsedTimeSeconds) override
	{
		if (!m_ui.ActiveSceneCamera)
			GetActiveCamera().Animate(fElapsedTimeSeconds);
		GetDeviceManager()->SetInformativeWindowTitle(g_WindowTitle);
	}

	void BackBufferResizing() override
	{
		m_Pipeline = nullptr;
	}

	void CopyActiveCameraToFirstPerson()
	{
		if (m_ui.ActiveSceneCamera)
		{
			dm::affine3 viewToWorld = m_ui.ActiveSceneCamera->GetViewToWorldMatrix();
			dm::float3 cameraPos = viewToWorld.m_translation;
			m_FirstPersonCamera.LookAt(cameraPos, cameraPos + viewToWorld.m_linear.row2, viewToWorld.m_linear.row1);
		}
		else if (m_ui.UseThirdPersonCamera)
		{
			m_FirstPersonCamera.LookAt(m_ThirdPersonCamera.GetPosition(), m_ThirdPersonCamera.GetPosition() + m_ThirdPersonCamera.GetDir(), m_ThirdPersonCamera.GetUp());
		}
	}

	BaseCamera& GetActiveCamera() const
	{
		return m_ui.UseThirdPersonCamera ? (BaseCamera&)m_ThirdPersonCamera : (BaseCamera&)m_FirstPersonCamera;
	}

	std::shared_ptr<IView>              m_View;
	std::shared_ptr<IView>              m_ViewPrevious;
	float                               m_CameraVerticalFov = 60.f;

	bool SetupView(nvrhi::IFramebuffer* framebuffer)
	{
		const nvrhi::FramebufferInfo& fbinfo = framebuffer->getFramebufferInfo();

		float2 renderTargetSize = float2(fbinfo.width, fbinfo.height);

		std::shared_ptr<StereoPlanarView> stereoView = std::dynamic_pointer_cast<StereoPlanarView, IView>(m_View);
		std::shared_ptr<PlanarView> planarView = std::dynamic_pointer_cast<PlanarView, IView>(m_View);

		dm::affine3 viewMatrix;
		float verticalFov = dm::radians(m_CameraVerticalFov);
		float zNear = 0.01f;
		if (m_ui.ActiveSceneCamera)
		{
			auto perspectiveCamera = std::dynamic_pointer_cast<PerspectiveCamera>(m_ui.ActiveSceneCamera);
			if (perspectiveCamera)
			{
				zNear = perspectiveCamera->zNear;
				verticalFov = perspectiveCamera->verticalFov;
			}

			viewMatrix = m_ui.ActiveSceneCamera->GetWorldToViewMatrix();
		}
		else
		{
			viewMatrix = GetActiveCamera().GetWorldToViewMatrix();
		}

		bool topologyChanged = false;

		{
			if (!planarView)
			{
				m_View = planarView = std::make_shared<PlanarView>();
				m_ViewPrevious = std::make_shared<PlanarView>();
				topologyChanged = true;
			}

			float4x4 projection = perspProjD3DStyleReverse(verticalFov, renderTargetSize.x / renderTargetSize.y, zNear);

			planarView->SetViewport(nvrhi::Viewport(renderTargetSize.x, renderTargetSize.y));

			planarView->SetMatrices(viewMatrix, projection);
			planarView->UpdateCache();

			m_ThirdPersonCamera.SetView(*planarView);

			if (topologyChanged)
			{
				*std::static_pointer_cast<PlanarView>(m_ViewPrevious) = *std::static_pointer_cast<PlanarView>(m_View);
			}
		}

		return topologyChanged;
	}

	virtual bool MouseScrollUpdate(double xoffset, double yoffset) override
	{
		GetActiveCamera().MouseScrollUpdate(xoffset, yoffset);

		return true;
	}

	virtual bool KeyboardUpdate(int key, int scancode, int action, int mods) override
	{
		if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		{
			m_ui.ShowUI = !m_ui.ShowUI;
			return true;
		}

		if (key == GLFW_KEY_GRAVE_ACCENT && action == GLFW_PRESS)
		{
			m_ui.ShowConsole = !m_ui.ShowConsole;
			return true;
		}

		if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
		{
			m_ui.EnableAnimations = !m_ui.EnableAnimations;
			return true;
		}

		if (key == GLFW_KEY_T && action == GLFW_PRESS)
		{
			CopyActiveCameraToFirstPerson();
			if (m_ui.ActiveSceneCamera)
			{
				m_ui.UseThirdPersonCamera = false;
				m_ui.ActiveSceneCamera = nullptr;
			}
			else
			{
				m_ui.UseThirdPersonCamera = !m_ui.UseThirdPersonCamera;
			}
			return true;
		}

		if (!m_ui.ActiveSceneCamera)
			GetActiveCamera().KeyboardUpdate(key, scancode, action, mods);
		return true;
	}

	uint2                               m_PickPosition = 0u;

	virtual bool MousePosUpdate(double xpos, double ypos) override
	{
		GetActiveCamera().MousePosUpdate(xpos, ypos);

		m_PickPosition = uint2(static_cast<uint>(xpos), static_cast<uint>(ypos));

		return true;
	}

	bool m_Pick = false;
	virtual bool MouseButtonUpdate(int button, int action, int mods) override
	{
		GetActiveCamera().MouseButtonUpdate(button, action, mods);

		if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_2)
			m_Pick = true;

		return true;
	}

	void Render(nvrhi::IFramebuffer* framebuffer) override
	{
		const nvrhi::FramebufferInfo& fbinfo = framebuffer->getFramebufferInfo();

		if (!m_Pipeline)
		{
			nvrhi::GraphicsPipelineDesc psoDesc;
			psoDesc.VS = m_VertexShader;
			psoDesc.PS = m_PixelShader;
			psoDesc.inputLayout = m_InputLayout;
			psoDesc.bindingLayouts = { m_BindingLayout };
			psoDesc.primType = nvrhi::PrimitiveType::TriangleList;
			psoDesc.renderState.depthStencilState.depthTestEnable = false;

			m_Pipeline = GetDevice()->createGraphicsPipeline(psoDesc, framebuffer);
		}

		m_CommandList->open();

		nvrhi::utils::ClearColorAttachment(m_CommandList, framebuffer, 0, nvrhi::Color(0.f));

		math::float4x4 projMatrix = math::perspProjD3DStyle(math::radians(60.f), float(fbinfo.width) / float(fbinfo.height), 0.1f, 10.f);

		SetupView(framebuffer);

		auto camera = GetActiveCamera();

		math::float4x4 viewProjMatrix = math::affineToHomogeneous(camera.GetWorldToViewMatrix()) * projMatrix;

		m_CommandList->writeBuffer(m_ConstantBuffer, &viewProjMatrix, sizeof(viewProjMatrix));

		nvrhi::GraphicsState state;
		state.bindings = { m_BindingSet };
		state.indexBuffer = { m_IndexBuffer, nvrhi::Format::R32_UINT, 0 };
		state.vertexBuffers = { { m_VertexBuffer, 0, 0 } };
		state.pipeline = m_Pipeline;
		state.framebuffer = framebuffer;
		state.viewport.addViewportAndScissorRect(fbinfo.getViewport());

		m_CommandList->setGraphicsState(state);

		nvrhi::DrawArguments args;
		args.vertexCount = g_Indices.size();
		m_CommandList->drawIndexed(args);

		m_CommandList->close();
		GetDevice()->executeCommandList(m_CommandList);
	}

	std::string GetCurrentSceneName() const
	{
		return m_CurrentSceneName;
	}

	std::vector<std::string> const& GetAvailableScenes() const
	{
		return m_SceneFilesAvailable;
	}

	std::vector<std::string>& GetAvailableScenes()
	{
		return m_SceneFilesAvailable;
	}

	using FEMType = StaticFEM2DAppP1;

	StaticFEM2DAppP1 fem_app;

	void BeginLoadingScene(const std::shared_ptr<vfs::RootFileSystem>& root_fs, const std::string& string)
	{
		OpenMesh::IO::read_mesh(fem_app.GetMesh(), string);
		auto mesh = fem_app.GetMesh();
		auto vertices = mesh.vertices();

		fem_app.evaluate();

		for (auto vertex : vertices)
		{
			auto point = mesh.point(vertex);

			auto faces = vertex.faces();

			float val = 0;
			if (faces.begin() != faces.end())
			{
				int id = faces.begin()->idx();
				val = fem_app.Value(id, Eigen::Vector3d(point[0], point[1], 0));
			}

			Vertex v = { {point[0],val,point[1] } };
			g_Vertices.push_back(v);
		}

		auto faces = mesh.faces();

		for (auto face : faces)
		{
			auto face_v = face.vertices();

			for (auto v : face_v)
			{
				g_Indices.push_back(v.idx());
			}
		}
	}

	void SetCurrentSceneName(const std::string& sceneName)
	{
		if (m_CurrentSceneName == sceneName)
			return;

		m_CurrentSceneName = sceneName;

		BeginLoadingScene(m_RootFs, m_CurrentSceneName);
	}

	std::shared_ptr<ShaderFactory> GetShaderFactory()
	{
		return m_ShaderFactory;
	}

private:
	UIData& m_ui;

	std::vector<std::string>            m_SceneFilesAvailable;

	FirstPersonCamera                   m_FirstPersonCamera;
	ThirdPersonCamera                   m_ThirdPersonCamera;
	std::string                         m_CurrentSceneName;

	std::shared_ptr<vfs::RootFileSystem> m_RootFs;
	std::shared_ptr<engine::ShaderFactory> m_ShaderFactory;
};

class UIRenderer : public app::ImGui_Renderer
{
private:
	std::shared_ptr<FEM2DVisualizer> m_app;

	std::unique_ptr<ImGui_Console> m_console;
	std::shared_ptr<engine::Light> m_SelectedLight;

	UIData& m_ui;
	nvrhi::CommandListHandle m_CommandList;

public:
	UIRenderer(DeviceManager* deviceManager, std::shared_ptr<FEM2DVisualizer> app, UIData& ui)
		: ImGui_Renderer(deviceManager)
		, m_app(app)
		, m_ui(ui)
	{
		m_CommandList = GetDevice()->createCommandList();

		ImGui_Console::Options opts;
		auto interpreter = std::make_shared<console::Interpreter>();
		// m_console = std::make_unique<ImGui_Console>(interpreter,opts);

		ImGui::GetIO().IniFilename = nullptr;
	}

protected:
	virtual void buildUI(void) override
	{
		if (!m_ui.ShowUI)
			return;

		const auto& io = ImGui::GetIO();

		int width, height;
		GetDeviceManager()->GetWindowDimensions(width, height);

		if (m_ui.ShowConsole && m_console)
		{
			m_console->Render(&m_ui.ShowConsole);
		}

		ImGui::SetNextWindowPos(ImVec2(10.f, 10.f), 0);
		ImGui::Begin("Settings", 0, ImGuiWindowFlags_AlwaysAutoResize);
		ImGui::Text("Renderer: %s", GetDeviceManager()->GetRendererString());
		double frameTime = GetDeviceManager()->GetAverageFrameTimeSeconds();
		if (frameTime > 0.0)
			ImGui::Text("%.3f ms/frame (%.1f FPS)", frameTime * 1e3, 1.0 / frameTime);

		const std::string currentScene = m_app->GetCurrentSceneName();
		if (ImGui::BeginCombo("Scene", currentScene.c_str()))
		{
			const std::vector<std::string>& scenes = m_app->GetAvailableScenes();
			for (const std::string& scene : scenes)
			{
				bool is_selected = scene == currentScene;
				if (ImGui::Selectable(scene.c_str(), is_selected))
					m_app->SetCurrentSceneName(scene);
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}

		if (ImGui::Button("Reload Shaders"))
			m_ui.ShaderReoladRequested = true;

		ImGui::Checkbox("VSync", &m_ui.EnableVsync);
		ImGui::Checkbox("Deferred Shading", &m_ui.UseDeferredShading);

		ImGui::Checkbox("Stereo", &m_ui.Stereo);
		ImGui::Checkbox("Animations", &m_ui.EnableAnimations);

		if (ImGui::BeginCombo("Camera (T)", m_ui.ActiveSceneCamera ? m_ui.ActiveSceneCamera->GetName().c_str()
			: m_ui.UseThirdPersonCamera ? "Third-Person" : "First-Person"))
		{
			if (ImGui::Selectable("First-Person", !m_ui.ActiveSceneCamera && !m_ui.UseThirdPersonCamera))
			{
				m_ui.ActiveSceneCamera.reset();
				m_ui.UseThirdPersonCamera = false;
			}
			if (ImGui::Selectable("Third-Person", !m_ui.ActiveSceneCamera && m_ui.UseThirdPersonCamera))
			{
				m_ui.ActiveSceneCamera.reset();
				m_ui.UseThirdPersonCamera = true;
				m_app->CopyActiveCameraToFirstPerson();
			}

			ImGui::EndCombo();
		}

		ImGui::SliderFloat("Ambient Intensity", &m_ui.AmbientIntensity, 0.f, 1.f);

		ImGui::Checkbox("Enable Light Probe", &m_ui.EnableLightProbe);
		if (m_ui.EnableLightProbe && ImGui::CollapsingHeader("Light Probe"))
		{
			ImGui::DragFloat("Diffuse Scale", &m_ui.LightProbeDiffuseScale, 0.01f, 0.0f, 10.0f);
			ImGui::DragFloat("Specular Scale", &m_ui.LightProbeSpecularScale, 0.01f, 0.0f, 10.0f);
		}

		ImGui::Checkbox("Enable Procedural Sky", &m_ui.EnableProceduralSky);
		if (m_ui.EnableProceduralSky && ImGui::CollapsingHeader("Sky Parameters"))
		{
			ImGui::SliderFloat("Brightness", &m_ui.SkyParams.brightness, 0.f, 1.f);
			ImGui::SliderFloat("Glow Size", &m_ui.SkyParams.glowSize, 0.f, 90.f);
			ImGui::SliderFloat("Glow Sharpness", &m_ui.SkyParams.glowSharpness, 1.f, 10.f);
			ImGui::SliderFloat("Glow Intensity", &m_ui.SkyParams.glowIntensity, 0.f, 1.f);
			ImGui::SliderFloat("Horizon Size", &m_ui.SkyParams.horizonSize, 0.f, 90.f);
		}
		ImGui::Checkbox("Enable SSAO", &m_ui.EnableSsao);
		ImGui::Checkbox("Enable Bloom", &m_ui.EnableBloom);
		ImGui::DragFloat("Bloom Sigma", &m_ui.BloomSigma, 0.01f, 0.1f, 100.f);
		ImGui::DragFloat("Bloom Alpha", &m_ui.BloomAlpha, 0.01f, 0.01f, 1.0f);
		ImGui::Checkbox("Enable Shadows", &m_ui.EnableShadows);
		ImGui::Checkbox("Enable Translucency", &m_ui.EnableTranslucency);

		ImGui::Separator();
		ImGui::Checkbox("Material Events", &m_ui.EnableMaterialEvents);
		ImGui::Separator();

		if (ImGui::Button("Screenshot"))
		{
			std::string fileName;
			if (FileDialog(false, "BMP files\0*.bmp\0All files\0*.*\0\0", fileName))
			{
				m_ui.ScreenshotFileName = fileName;
			}
		}

		ImGui::End();

		if (!m_ui.UseDeferredShading)
			m_ui.EnableSsao = false;
	}
};

int main(int argc, const char** argv)
{
	nvrhi::GraphicsAPI api = app::GetGraphicsAPIFromCommandLine(argc, argv);
	app::DeviceManager* deviceManager = app::DeviceManager::Create(api);

	app::DeviceCreationParameters deviceParams;
#ifdef _DEBUG
	deviceParams.enableDebugRuntime = true;
	deviceParams.enableNvrhiValidationLayer = true;
#endif

	if (!deviceManager->CreateWindowDeviceAndSwapChain(deviceParams, g_WindowTitle))
	{
		log::fatal("Cannot initialize a graphics device with the requested parameters");
		return 1;
	}

	{
		UIData uiData;

		std::shared_ptr<FEM2DVisualizer> example = std::make_shared<FEM2DVisualizer>(deviceManager, uiData);
		if (example->Init())
		{
			std::shared_ptr<UIRenderer> gui = std::make_shared<UIRenderer>(deviceManager, example, uiData);

			gui->Init(example->GetShaderFactory());

			deviceManager->AddRenderPassToBack(example.get());
			deviceManager->AddRenderPassToBack(gui.get());

			deviceManager->RunMessageLoop();
			deviceManager->RemoveRenderPass(example.get());
		}
	}

	deviceManager->Shutdown();

	delete deviceManager;

	return 0;
}