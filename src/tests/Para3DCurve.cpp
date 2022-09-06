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

#include "donut/engine/BindingCache.h"
#include "Geometry/ParameterDesign/Bezier.hpp"
#include "Geometry/ParameterDesign/ParameterCurve.h"

using namespace donut;
using namespace math;
using namespace app;
using namespace vfs;
using namespace engine;

static const char* g_WindowTitle = "Donut Example: Vertex Buffer";

struct Vertex
{
	float3 position;
	float3 normal;
};

static const Vertex Vertices[] = {
	{ {-0.5f,  0.5f, -0.5f}, {0.0f, 0.0f,-1.0f} }, // front face
	{ { 0.5f, -0.5f, -0.5f}, {0.0f, 0.0f,-1.0f} },
	{ {-0.5f, -0.5f, -0.5f}, {0.0f, 0.0f,-1.0f} },
	{ { 0.5f,  0.5f, -0.5f}, {0.0f, 0.0f,-1.0f} },

	{ { 0.5f, -0.5f, -0.5f}, {1.0f,0.0f,0.0f} }, // right side face
	{ { 0.5f,  0.5f,  0.5f}, {1.0f,0.0f,0.0f} },
	{ { 0.5f, -0.5f,  0.5f}, {1.0f,0.0f,0.0f} },
	{ { 0.5f,  0.5f, -0.5f}, {1.0f,0.0f,0.0f} },

	{ {-0.5f,  0.5f,  0.5f}, {-1.0,0.0f, 0.0f} }, // left side face
	{ {-0.5f, -0.5f, -0.5f}, {-1.0,0.0f, 0.0f} },
	{ {-0.5f, -0.5f,  0.5f}, {-1.0,0.0f, 0.0f} },
	{ {-0.5f,  0.5f, -0.5f}, {-1.0,0.0f, 0.0f} },

	{ { 0.5f,  0.5f,  0.5f}, {0.0f,0.0f,1.0f} }, // back face
	{ {-0.5f, -0.5f,  0.5f}, {0.0f,0.0f,1.0f} },
	{ { 0.5f, -0.5f,  0.5f}, {0.0f,0.0f,1.0f} },
	{ {-0.5f,  0.5f,  0.5f}, {0.0f,0.0f,1.0f} },

	{ {-0.5f,  0.5f, -0.5f}, {0.0f,-1.0f,0.0f} }, // top face
	{ { 0.5f,  0.5f,  0.5f}, {0.0f,-1.0f,0.0f} },
	{ { 0.5f,  0.5f, -0.5f}, {0.0f,-1.0f,0.0f} },
	{ {-0.5f,  0.5f,  0.5f}, {0.0f,-1.0f,0.0f} },

	{ { 0.5f, -0.5f,  0.5f},  {0.0f,1.0f,0.0f} }, // bottom face
	{ {-0.5f, -0.5f, -0.5f}, {0.0f,1.0f,0.0f} },
	{ { 0.5f, -0.5f, -0.5f}, {0.0f,1.0f,0.0f} },
	{ {-0.5f, -0.5f,  0.5f}, {0.0f,1.0f,0.0f} },
};

static const uint32_t Indices[] = {
	 0,  1,  2,   0,  3,  1, // front face
	 4,  5,  6,   4,  7,  5, // left face
	 8,  9, 10,   8, 11,  9, // right face
	12, 13, 14,  12, 15, 13, // back face
	16, 17, 18,  16, 19, 17, // top face
	20, 21, 22,  20, 23, 21, // bottom face
};

struct UIData
{
	bool ShowUI = true;
	bool ShowConsole = false;
	bool UseDeferredShading = true;
	bool Stereo = false;
	bool EnableSsao = true;

	render::SkyParameters SkyParams;

	bool EnableVsync = true;
	bool ShaderReoladRequested = false;
	bool EnableProceduralSky = true;
	bool EnableBloom = true;
	float BloomSigma = 32.f;
	float BloomAlpha = 0.05f;
	bool EnableTranslucency = true;
	bool EnableMaterialEvents = false;
	bool EnableShadows = true;
	float AmbientIntensity = 1.0f;
	bool EnableLightProbe = true;
	float LightProbeDiffuseScale = 1.f;
	float LightProbeSpecularScale = 1.f;
	float CsmExponent = 4.f;
	bool DisplayShadowMap = false;
	bool UseThirdPersonCamera = true;
	bool EnableAnimations = false;
	std::shared_ptr<Material> SelectedMaterial;
	std::shared_ptr<SceneGraphNode> SelectedNode;
	std::string ScreenshotFileName;
	std::shared_ptr<SceneCamera> ActiveSceneCamera;
};

auto rhs_func = [](Eigen::Vector3d vector)
{
	Float x = vector.x();
	Float y = vector.y();

	return -10 * (-2 * cos(1 - y) * cos(y) * sin(1 - x) * sin(x) - 2 * cos(1 - x) * cos(x) * sin(1 - y) * sin(y) -
		4 * sin(1 - x) * sin(x) * sin(1 - y) * sin(y));
};

struct CameraPara
{
	float4x4 projview;
	float3 cam_pos;
};

class FEM2DVisualizer : public IRenderPass
{
	std::vector<Vertex> g_Vertices;

	std::vector<uint32_t> g_Indices;

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

	FEM2DVisualizer(DeviceManager* deviceManager, UIData& ui) : IRenderPass(deviceManager),
		m_ui(ui)
	{
	}

	std::shared_ptr<NativeFileSystem> nativeFS;

	void LoadSceneToGPU()
	{
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
	}

	std::shared_ptr<CommonRenderPasses> m_CommonPasses;
	std::unique_ptr<BindingCache> m_BindingCache;

	bool Init()
	{
		nativeFS = std::make_shared<NativeFileSystem>();

		std::filesystem::path frameworkShaderPath = GetDirectoryWithExecutable() / "shaders/framework" /
			GetShaderTypeName(GetDevice()->getGraphicsAPI());
		std::filesystem::path appShaderPath = GetDirectoryWithExecutable() / "shaders/Numerical_test_Para3DCurve" /
			GetShaderTypeName(GetDevice()->getGraphicsAPI());

		m_RootFs = std::make_shared<RootFileSystem>();
		m_RootFs->mount("/shaders/donut", frameworkShaderPath);
		m_RootFs->mount("/shaders/app", appShaderPath);

		BeginLoadingScene(m_RootFs, "");

		m_ShaderFactory = std::make_shared<ShaderFactory>(GetDevice(), m_RootFs, "/shaders");
		m_CommonPasses = std::make_shared<CommonRenderPasses>(GetDevice(), m_ShaderFactory);
		m_BindingCache = std::make_unique<BindingCache>(GetDevice());

		m_VertexShader = m_ShaderFactory->CreateShader("app/shaders.hlsl", "main_vs", nullptr,
			nvrhi::ShaderType::Vertex);
		m_PixelShader = m_ShaderFactory->CreateShader("app/shaders.hlsl", "main_ps", nullptr, nvrhi::ShaderType::Pixel);

		if (!m_VertexShader || !m_PixelShader)
		{
			return false;
		}

		m_ConstantBuffer = GetDevice()->createBuffer(nvrhi::utils::CreateVolatileConstantBufferDesc(
			sizeof(CameraPara), "ConstantBuffer",
			c_MaxRenderPassConstantBufferVersions));

		nvrhi::VertexAttributeDesc attributes[] = {
			nvrhi::VertexAttributeDesc()
			.setName("POSITION")
			.setFormat(nvrhi::Format::RGB32_FLOAT)
			.setOffset(offsetof(Vertex, position))
			.setElementStride(sizeof(Vertex)),
			nvrhi::VertexAttributeDesc()
			.setName("NORMAL")
			.setFormat(nvrhi::Format::RGB32_FLOAT)
			.setOffset(offsetof(Vertex, normal))
			.setElementStride(sizeof(Vertex)),
		};

		m_InputLayout = GetDevice()->createInputLayout(attributes, static_cast<uint32_t>(std::size(attributes)),
			m_VertexShader);

		CommonRenderPasses commonPasses(GetDevice(), m_ShaderFactory);
		TextureCache textureCache(GetDevice(), nativeFS, nullptr);

		LoadSceneToGPU();

		nvrhi::BindingSetDesc bindingSetDesc;
		bindingSetDesc.bindings = {
			nvrhi::BindingSetItem::ConstantBuffer(0, m_ConstantBuffer),
		};

		if (!nvrhi::utils::CreateBindingSetAndLayout(GetDevice(), nvrhi::ShaderType::All, 0, bindingSetDesc,
			m_BindingLayout, m_BindingSet))
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
			affine3 viewToWorld = m_ui.ActiveSceneCamera->GetViewToWorldMatrix();
			float3 cameraPos = viewToWorld.m_translation;
			m_FirstPersonCamera.LookAt(cameraPos, cameraPos + viewToWorld.m_linear.row2, viewToWorld.m_linear.row1);
		}
		else if (m_ui.UseThirdPersonCamera)
		{
			m_FirstPersonCamera.LookAt(m_ThirdPersonCamera.GetPosition(),
				m_ThirdPersonCamera.GetPosition() + m_ThirdPersonCamera.GetDir(),
				m_ThirdPersonCamera.GetUp());
		}
	}

	BaseCamera& GetActiveCamera() const
	{
		return m_ui.UseThirdPersonCamera ? (BaseCamera&)m_ThirdPersonCamera : (BaseCamera&)m_FirstPersonCamera;
	}

	std::shared_ptr<IView> m_View;
	std::shared_ptr<IView> m_ViewPrevious;
	float m_CameraVerticalFov = 60.f;

	bool SetupView(nvrhi::IFramebuffer* framebuffer)
	{
		const nvrhi::FramebufferInfo& fbinfo = framebuffer->getFramebufferInfo();

		auto renderTargetSize = float2(fbinfo.width, fbinfo.height);

		std::shared_ptr<StereoPlanarView> stereoView = std::dynamic_pointer_cast<StereoPlanarView, IView>(m_View);
		std::shared_ptr<PlanarView> planarView = std::dynamic_pointer_cast<PlanarView, IView>(m_View);

		affine3 viewMatrix;
		float verticalFov = radians(m_CameraVerticalFov);
		float zNear = 0.01f;
		float zFar = 100.f;
		if (m_ui.ActiveSceneCamera)
		{
			auto perspectiveCamera = std::dynamic_pointer_cast<PerspectiveCamera>(m_ui.ActiveSceneCamera);
			if (perspectiveCamera)
			{
				zNear = perspectiveCamera->zNear;
				if (perspectiveCamera->zFar.has_value())
				{
					zFar = perspectiveCamera->zFar.value();
				}
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

			float4x4 projection = perspProjD3DStyle(verticalFov, renderTargetSize.x / renderTargetSize.y, zNear, zFar);

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

	bool MouseScrollUpdate(double xoffset, double yoffset) override
	{
		GetActiveCamera().MouseScrollUpdate(xoffset, yoffset);

		return true;
	}

	bool KeyboardUpdate(int key, int scancode, int action, int mods) override
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

	uint2 m_PickPosition = 0u;

	bool MousePosUpdate(double xpos, double ypos) override
	{
		GetActiveCamera().MousePosUpdate(xpos, ypos);

		m_PickPosition = uint2(static_cast<uint>(xpos), static_cast<uint>(ypos));

		return true;
	}

	bool m_Pick = false;

	bool MouseButtonUpdate(int button, int action, int mods) override
	{
		GetActiveCamera().MouseButtonUpdate(button, action, mods);

		if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_2)
			m_Pick = true;

		return true;
	}

	nvrhi::TextureHandle m_DepthBuffer;
	nvrhi::TextureHandle m_ColorBuffer;
	nvrhi::FramebufferHandle m_Framebuffer;

	void Render(nvrhi::IFramebuffer* framebuffer) override
	{
		const nvrhi::FramebufferInfo& fbinfo = framebuffer->getFramebufferInfo();

		if (!m_Pipeline)
		{
			nvrhi::TextureDesc textureDesc;
			textureDesc.format = nvrhi::Format::SRGBA8_UNORM;
			textureDesc.isRenderTarget = true;
			textureDesc.initialState = nvrhi::ResourceStates::RenderTarget;
			textureDesc.keepInitialState = true;
			textureDesc.clearValue = nvrhi::Color(0.f);
			textureDesc.useClearValue = true;
			textureDesc.debugName = "ColorBuffer";
			textureDesc.width = fbinfo.width;
			textureDesc.height = fbinfo.height;
			textureDesc.dimension = nvrhi::TextureDimension::Texture2D;
			m_ColorBuffer = GetDevice()->createTexture(textureDesc);

			textureDesc.format = nvrhi::Format::D24S8;
			textureDesc.debugName = "DepthBuffer";
			textureDesc.initialState = nvrhi::ResourceStates::DepthWrite;
			m_DepthBuffer = GetDevice()->createTexture(textureDesc);

			nvrhi::FramebufferDesc framebufferDesc;
			framebufferDesc.addColorAttachment(m_ColorBuffer, nvrhi::AllSubresources);
			framebufferDesc.setDepthAttachment(m_DepthBuffer);
			m_Framebuffer = GetDevice()->createFramebuffer(framebufferDesc);

			nvrhi::GraphicsPipelineDesc psoDesc;
			psoDesc.VS = m_VertexShader;
			psoDesc.PS = m_PixelShader;
			psoDesc.inputLayout = m_InputLayout;
			psoDesc.bindingLayouts = { m_BindingLayout };
			psoDesc.primType = nvrhi::PrimitiveType::TriangleList;
			psoDesc.renderState.depthStencilState.depthTestEnable = true;
			psoDesc.renderState.depthStencilState.depthFunc = nvrhi::ComparisonFunc::GreaterOrEqual;
			//psoDesc.renderState.rasterState.frontCounterClockwise = true;
			psoDesc.renderState.rasterState.setCullNone();
			m_Pipeline = GetDevice()->createGraphicsPipeline(psoDesc, m_Framebuffer);
		}

		m_CommandList->open();

		m_CommandList->clearTextureFloat(m_ColorBuffer, nvrhi::AllSubresources, nvrhi::Color(0.f));
		m_CommandList->clearDepthStencilTexture(m_DepthBuffer, nvrhi::AllSubresources, true, 0.f, true, 0);

		//nvrhi::utils::ClearColorAttachment(m_CommandList, framebuffer, 0, nvrhi::Color(1.0f));
		//nvrhi::utils::ClearDepthStencilAttachment(m_CommandList, framebuffer, 0, 0.0);

		float4x4 projMatrix = perspProjD3DStyleReverse(radians(60.f),
			static_cast<float>(fbinfo.width) / static_cast<float>(fbinfo.
				height), 0.1f);

		SetupView(framebuffer);

		auto camera = GetActiveCamera();

		CameraPara para;

		para.projview = affineToHomogeneous(camera.GetWorldToViewMatrix()) * projMatrix;
		para.cam_pos = camera.GetPosition();

		m_CommandList->writeBuffer(m_ConstantBuffer, &para, sizeof(para));

		nvrhi::GraphicsState state;
		state.bindings = { m_BindingSet };
		state.indexBuffer = { m_IndexBuffer, nvrhi::Format::R32_UINT, 0 };
		state.vertexBuffers = { {m_VertexBuffer, 0, 0} };
		state.pipeline = m_Pipeline;
		state.framebuffer = m_Framebuffer;
		state.viewport.addViewportAndScissorRect(fbinfo.getViewport());

		m_CommandList->setGraphicsState(state);

		nvrhi::DrawArguments args;
		args.vertexCount = g_Indices.size();
		m_CommandList->drawIndexed(args);
		m_CommonPasses->BlitTexture(m_CommandList, framebuffer, m_ColorBuffer, m_BindingCache.get());

		m_CommandList->close();
		GetDevice()->executeCommandList(m_CommandList);
	}

	std::string GetCurrentSceneName() const
	{
		return m_CurrentSceneName;
	}

	const std::vector<std::string>& GetAvailableScenes() const
	{
		return m_SceneFilesAvailable;
	}

	std::vector<std::string>& GetAvailableScenes()
	{
		return m_SceneFilesAvailable;
	}

	using FEMType = StaticFEM2DAppP1;

	void BeginLoadingScene(const std::shared_ptr<RootFileSystem>& root_fs, const std::string& string)
	{
		g_Indices.clear();
		g_Vertices.clear();
		std::vector<Eigen::Vector3d> points{ Eigen::Vector3d(2,0,1),Eigen::Vector3d(2,1,1),Eigen::Vector3d(0,2,2) };
		BezierCurveND<3> bezier(points);
		bezier.evaluate();

		for (int i = 0; i < sizeof(Vertices) / sizeof(Vertex); ++i)
		{
			g_Vertices.push_back(Vertices[i]);
		}

		for (int i = 0; i < sizeof(Indices) / sizeof(int); ++i)
		{
			g_Indices.push_back(Indices[i]);
		}
	}

	std::shared_ptr<ShaderFactory> GetShaderFactory()
	{
		return m_ShaderFactory;
	}

private:
	UIData& m_ui;

	std::vector<std::string> m_SceneFilesAvailable;

	FirstPersonCamera m_FirstPersonCamera;
	ThirdPersonCamera m_ThirdPersonCamera;
	std::string m_CurrentSceneName;

	std::shared_ptr<RootFileSystem> m_RootFs;
	std::shared_ptr<ShaderFactory> m_ShaderFactory;
};

class UIRenderer : public ImGui_Renderer
{
private:
	std::shared_ptr<FEM2DVisualizer> m_app;

	std::unique_ptr<ImGui_Console> m_console;
	std::shared_ptr<Light> m_SelectedLight;

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
	void buildUI(void) override
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
		ImGui::Begin("Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
		ImGui::Text("Renderer: %s", GetDeviceManager()->GetRendererString());
		double frameTime = GetDeviceManager()->GetAverageFrameTimeSeconds();
		if (frameTime > 0.0)
			ImGui::Text("%.3f ms/frame (%.1f FPS)", frameTime * 1e3, 1.0 / frameTime);

		const std::string currentScene = m_app->GetCurrentSceneName();

		if (ImGui::Button("Reload Shaders"))
			m_ui.ShaderReoladRequested = true;

		ImGui::Checkbox("VSync", &m_ui.EnableVsync);
		ImGui::Checkbox("Deferred Shading", &m_ui.UseDeferredShading);

		ImGui::Checkbox("Stereo", &m_ui.Stereo);
		ImGui::Checkbox("Animations", &m_ui.EnableAnimations);

		if (ImGui::BeginCombo("Camera (T)", m_ui.ActiveSceneCamera
			? m_ui.ActiveSceneCamera->GetName().c_str()
			: m_ui.UseThirdPersonCamera
			? "Third-Person"
			: "First-Person"))
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
	nvrhi::GraphicsAPI api = GetGraphicsAPIFromCommandLine(argc, argv);
	DeviceManager* deviceManager = DeviceManager::Create(api);

	DeviceCreationParameters deviceParams;
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

		auto example = std::make_shared<FEM2DVisualizer>(deviceManager, uiData);
		if (example->Init())
		{
			auto gui = std::make_shared<UIRenderer>(deviceManager, example, uiData);

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