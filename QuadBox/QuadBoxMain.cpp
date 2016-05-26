#include "pch.h"
#include "QuadBoxMain.h"
#include "Common\DirectXHelper.h"
#include "Content\camera.h"

using namespace QuadBox;
using namespace Windows::Foundation;
using namespace Windows::System::Threading;
using namespace Concurrency;
using namespace Windows::UI::Core;
using namespace Windows::System;
using namespace Windows::Gaming::Input;
using namespace Windows::Devices::Input;

CCameraObject g_Camera;

// Loads and initializes application assets when the application is loaded.
QuadBoxMain::QuadBoxMain(const std::shared_ptr<DX::DeviceResources>& deviceResources) :
	m_deviceResources(deviceResources)
{
	// Register to be notified if the Device is lost or recreated
	m_deviceResources->RegisterDeviceNotify(this);

	// TODO: Replace this with your app's content initialization.
	m_sceneRenderer = std::unique_ptr<Sample3DSceneRenderer>(new Sample3DSceneRenderer(m_deviceResources));

	m_fpsTextRenderer = std::unique_ptr<SampleFpsTextRenderer>(new SampleFpsTextRenderer(m_deviceResources));

	// TODO: Change the timer settings if you want something other than the default variable timestep mode.
	// e.g. for 60 FPS fixed timestep update logic, call:
	/*
	m_timer.SetFixedTimeStep(true);
	m_timer.SetTargetElapsedSeconds(1.0 / 60);
	*/

	CreateWindowSizeDependentResources();

	CoreWindow^ curWindow = CoreWindow::GetForCurrentThread();
	curWindow->SetPointerCapture();
}

void QuadBoxMain::OnMouseMove(MouseDevice^ sender, MouseEventArgs^ args)
{
	
}

QuadBoxMain::~QuadBoxMain()
{
	// Deregister device notification
	m_deviceResources->RegisterDeviceNotify(nullptr);
}

// Updates application state when the window size changes (e.g. device orientation change)
void QuadBoxMain::CreateWindowSizeDependentResources() 
{
	// TODO: Replace this with the size-dependent initialization of your app's content.
	m_sceneRenderer->CreateWindowSizeDependentResources();

	Size outputSize = m_deviceResources->GetOutputSize();
	float aspectRatio = outputSize.Width / outputSize.Height;
	float fovAngleY = 70.0f * 3.14159F / 180.0f;
	if (aspectRatio < 1.0f)
	{
		fovAngleY *= 2.0f;
	}

	vec3f camPos(0.0f, 0.7f, 1.5f);
	vec3f lookPos(0.0f, -0.1f, 0.0f);
	vec3f up(0.0f, 1.0f, 0.0f);
	g_Camera.SetCamera(camPos, lookPos, up);
	g_Camera.SetPerspective(fovAngleY, aspectRatio, 0.1f, 100.0f);

	mat4 V, P;
	g_Camera.GetViewMatrix(V);
	g_Camera.GetProjectionMatrix(P);

	m_sceneRenderer->SetViewMatrix(V);
	m_sceneRenderer->SetProjectionMatrix(P);
}

// Updates the application state once per frame.
void QuadBoxMain::Update() 
{
	ProcessInput();

	// Update scene objects.
	m_timer.Tick([&]()
	{
		// TODO: Replace this with your app's content update functions.
		m_sceneRenderer->Update(m_timer);
		m_fpsTextRenderer->Update(m_timer);
	});
}

// Process all input from the user before updating game state
void QuadBoxMain::ProcessInput()
{
	const float DEAD_ZONE = 0.45f;
	auto count = Gamepad::Gamepads->Size;
	
	if (count > 0)
	{
		auto curGamepad = Gamepad::Gamepads->First()->Current;
		auto reading = curGamepad->GetCurrentReading();

		if (reading.LeftThumbstickX < -DEAD_ZONE)
			g_Camera.MoveLeft();
		else if (reading.LeftThumbstickX > DEAD_ZONE)
			g_Camera.MoveRight();

		if (reading.LeftThumbstickY > DEAD_ZONE)
			g_Camera.MoveForward();
		else if (reading.LeftThumbstickY < -DEAD_ZONE)
			g_Camera.MoveBack();
	}

	CoreWindow^ curWindow = CoreWindow::GetForCurrentThread();
	if(curWindow->GetAsyncKeyState(VirtualKey::W) == CoreVirtualKeyStates::Down)
	{
		g_Camera.MoveForward();
	}

	if (curWindow->GetAsyncKeyState(VirtualKey::S) == CoreVirtualKeyStates::Down)
	{
		g_Camera.MoveBack();
	}

	if (curWindow->GetAsyncKeyState(VirtualKey::A) == CoreVirtualKeyStates::Down)
	{
		g_Camera.MoveLeft();
	}

	if (curWindow->GetAsyncKeyState(VirtualKey::D) == CoreVirtualKeyStates::Down)
	{
		g_Camera.MoveRight();
	}
}

// Renders the current frame according to the current application state.
// Returns true if the frame was rendered and is ready to be displayed.
bool QuadBoxMain::Render() 
{
	// Don't try to render anything before the first Update.
	if (m_timer.GetFrameCount() == 0)
	{
		return false;
	}

	auto context = m_deviceResources->GetD3DDeviceContext();

	// Reset the viewport to target the whole screen.
	auto viewport = m_deviceResources->GetScreenViewport();
	context->RSSetViewports(1, &viewport);

	// Reset render targets to the screen.
	ID3D11RenderTargetView *const targets[1] = { m_deviceResources->GetBackBufferRenderTargetView() };
	context->OMSetRenderTargets(1, targets, m_deviceResources->GetDepthStencilView());

	// Clear the back buffer and depth stencil view.
	context->ClearRenderTargetView(m_deviceResources->GetBackBufferRenderTargetView(), DirectX::Colors::CornflowerBlue);
	context->ClearDepthStencilView(m_deviceResources->GetDepthStencilView(), D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0f, 0);

	// Update the camera matrices 
	mat4 V, P;
	g_Camera.GetViewMatrix(V);
	g_Camera.GetProjectionMatrix(P);

	m_sceneRenderer->SetViewMatrix(V);
	m_sceneRenderer->SetProjectionMatrix(P);

	// Render the scene objects.
	// TODO: Replace this with your app's content rendering functions.
	m_sceneRenderer->Render();
	m_fpsTextRenderer->Render();

	return true;
}

// Notifies renderers that device resources need to be released.
void QuadBoxMain::OnDeviceLost()
{
	m_sceneRenderer->ReleaseDeviceDependentResources();
	m_fpsTextRenderer->ReleaseDeviceDependentResources();
}

// Notifies renderers that device resources may now be recreated.
void QuadBoxMain::OnDeviceRestored()
{
	m_sceneRenderer->CreateDeviceDependentResources();
	m_fpsTextRenderer->CreateDeviceDependentResources();
	CreateWindowSizeDependentResources();
}
