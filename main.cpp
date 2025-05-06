#include <QApplication>
#include <QColorSpace>
#include <QCommandLineParser>
#include <QSurfaceFormat>

#include <QComboBox>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QPushButton>
#include <QTabWidget>
#include <QCheckBox>

#include "GLMainWindow.hpp"
#include "ExampleRenderer.hpp"


#ifdef _WIN32
// always use (proper) hardware acceleration if available, since Intel's iGPUs have extremely buggy OpenGL support
extern "C"
{
	__declspec(dllexport) DWORD NvOptimusEnablement = 1;
	__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
}
#endif

void addExampleTab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto tab = new QWidget{tabWidget};
	tabWidget->addTab(tab, "Beispiel");
	QWidget::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index != current)
			return;

		mainWindow->setRendererFactory(
			[] (QObject * parent) {
				// ExampleRenderer contains our OpenGL code
				return new ExampleRenderer{parent};
			}
		);
	});
}



void addSimulationControls(GLMainWindow * mainWindow)
{
	auto dock = new QDockWidget("Simulation controls", mainWindow);
	dock->setAllowedAreas(Qt::DockWidgetArea::RightDockWidgetArea);

	auto tabWidget = new QTabWidget(dock);
	tabWidget->setTabPosition(QTabWidget::TabPosition::West);
	dock->setWidget(tabWidget);

	addExampleTab(mainWindow, tabWidget);
	// * *
	// TODO: add a tabs for each programming task
	// * *

	mainWindow->addDockWidget(Qt::DockWidgetArea::RightDockWidgetArea, dock);
}

int main(int argc, char ** argv)
{
	// set up OpenGL surface format
	auto surfaceFormat = QSurfaceFormat::defaultFormat();
	surfaceFormat.setVersion(3, 3);
	surfaceFormat.setProfile(QSurfaceFormat::CoreProfile);
	surfaceFormat.setOption(QSurfaceFormat::DebugContext);
	surfaceFormat.setColorSpace(QColorSpace::NamedColorSpace::SRgb);
	surfaceFormat.setSamples(4);
	QSurfaceFormat::setDefaultFormat(surfaceFormat);

	// set up a Qt application
	using App = QApplication;
	App::setAttribute(Qt::AA_UseDesktopOpenGL);
	App::setApplicationName("SimulationFramework");
	App::setApplicationDisplayName(App::translate("main", "Simulation Framework"));
	App::setApplicationVersion("1.0");
	App app(argc, argv);

	// configure command line parser
	QCommandLineParser parser;
	parser.setApplicationDescription(App::translate("main", "Simulation framework for the TU Darmstadt lecture on physically based simulation and animation."));
	parser.addHelpOption();
	parser.addVersionOption();

	// provide a flag for OpenGL debugging (outputs can be viewed with debugger!)
	QCommandLineOption debugGLOption({ "g", "debug-gl" }, App::translate("main", "Enable OpenGL debug logging"));
	parser.addOption(debugGLOption);

	// parse command line
	parser.process(app);

	// create main window. modify GLMainWindow.ui to add widgets etc.
	GLMainWindow widget;

	// enable OpenGL error logging (look at debugger output!) if flag is passed
	if(parser.isSet(debugGLOption))
	{
		widget.setOpenGLLoggingSynchronous(true);
		widget.setOpenGLLoggingEnabled(true);
	}

	// set up which renderer to use. factory to create renderer when OpenGL context exists
	widget.setRendererFactory(
		[] (QObject * parent) {
			// ExampleRenderer contains our OpenGL code
			return new ExampleRenderer{parent};
		}
	);

	addSimulationControls(&widget);

	// show the main window
	widget.show();

	// run the event loop (do not write your own!)
	return app.exec();
}
