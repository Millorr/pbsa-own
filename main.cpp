#include <QApplication>
#include <QColorSpace>
#include <QCommandLineParser>
#include <QSurfaceFormat>
#include <QDebug>

#include <QComboBox>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QPushButton>
#include <QTabWidget>
#include <QCheckBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>

#include "GLMainWindow.hpp"
#include "ExampleRenderer.hpp"
#include "Prog1_1Simulation.hpp"
#include "Prog1_2Simulation.hpp"
#include "Prog1_3Simulation.hpp"
#include "Prog2_1Simulation.hpp"
#include "Prog2_2Simulation.hpp"
#include "Prog2_3Simulation.hpp"

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

void addEx1_1Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto tab = new QWidget{tabWidget};
	tabWidget->addTab(tab, "Exercise 1.1");

	// Layouts
	auto mainLayout = new QVBoxLayout{tab};
	auto controlLayout = new QVBoxLayout{};
	auto controlGroup = new QGroupBox{"Simulation Controls"};
	controlGroup->setLayout(controlLayout);

	// r0 (initial position)
	auto r0Layout = new QHBoxLayout{};
	auto r0Label = new QLabel{"r₀ (initial position)"};
	auto r0SpinX = new QDoubleSpinBox{};
	auto r0SpinY = new QDoubleSpinBox{};
	auto r0SpinZ = new QDoubleSpinBox{};
	for(auto* s : {r0SpinX, r0SpinY, r0SpinZ}) {
		s->setRange(-1000.0, 1000.0);
		s->setDecimals(4);
	}
	r0SpinX->setValue(0.0);
	r0SpinY->setValue(15.0);
	r0SpinZ->setValue(0.0);
	r0Layout->addWidget(r0Label);
	r0Layout->addWidget(r0SpinX);
	r0Layout->addWidget(r0SpinY);
	r0Layout->addWidget(r0SpinZ);
	controlLayout->addLayout(r0Layout);

	// v0 (initial velocity)
	auto v0Layout = new QHBoxLayout{};
	auto v0Label = new QLabel{"v₀ (initial velocity)"};
	auto v0SpinX = new QDoubleSpinBox{};
	auto v0SpinY = new QDoubleSpinBox{};
	auto v0SpinZ = new QDoubleSpinBox{};
	for(auto* s : {v0SpinX, v0SpinY, v0SpinZ}) {
		s->setRange(-1000.0, 1000.0);
		s->setDecimals(4);
	}
	v0SpinX->setValue(20.0);
	v0SpinY->setValue(0.0);
	v0SpinZ->setValue(0.0);
	v0Layout->addWidget(v0Label);
	v0Layout->addWidget(v0SpinX);
	v0Layout->addWidget(v0SpinY);
	v0Layout->addWidget(v0SpinZ);
	controlLayout->addLayout(v0Layout);

	// g (Gravitationskonstante * M)
	auto gLayout = new QHBoxLayout{};
	auto gLabel = new QLabel{"g = G * M"};
	auto gSpin = new QDoubleSpinBox{};
	gSpin->setRange(0.0001, 1000000.0);
	gSpin->setDecimals(6);
	gSpin->setValue(20.0);
	gLayout->addWidget(gLabel);
	gLayout->addWidget(gSpin);
	controlLayout->addLayout(gLayout);

	// delta t (Zeitschritt)
	auto dtLayout = new QHBoxLayout{};
	auto dtLabel = new QLabel{"Δt (Zeitschritt)"};
	auto dtSpin = new QDoubleSpinBox{};
	dtSpin->setRange(0.0001, 100.0);
	dtSpin->setDecimals(5);
	dtSpin->setValue(0.01);
	dtLayout->addWidget(dtLabel);
	dtLayout->addWidget(dtSpin);
	controlLayout->addLayout(dtLayout);

	// Integrationsschema Auswahl
	auto methodLayout = new QHBoxLayout{};
	auto methodLabel = new QLabel{"Integrationsverfahren"};
	auto methodCombo = new QComboBox{};
	methodCombo->addItem("Expliziter Euler");
	methodCombo->addItem("Impliziter Euler");
	methodCombo->addItem("Verlet");
	methodLayout->addWidget(methodLabel);
	methodLayout->addWidget(methodCombo);
	controlLayout->addLayout(methodLayout);

	// Reset-Button
	auto resetButton = new QPushButton{"Reset"};
	controlLayout->addWidget(resetButton);

	auto glWidget = new QWidget{}; // Platzhalter
	// -> Das tatsächliche Widget wird über das RendererFactory erstellt und dann hier gesetzt

	mainLayout->addWidget(controlGroup);
	mainLayout->addWidget(glWidget, 1);
	QWidget::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index != current)
			return;

		mainWindow->setRendererFactory(
			[=] (QObject * parent) {
				auto sim = new Prog1_1Simulation{};
				sim->setParent(parent);

				QObject::connect(resetButton, &QPushButton::clicked, [=]() {
					sim->reset(
						static_cast<Prog1_1Simulation::Integration>(methodCombo->currentIndex()),
						Eigen::Vector3d(r0SpinX->value(), r0SpinY->value(), r0SpinZ->value()),
						Eigen::Vector3d(v0SpinX->value(), v0SpinY->value(), v0SpinZ->value()),
						gSpin->value(),
						dtSpin->value()
					);
				});
				return sim;
			}
		);
	});
}

void addEx1_2Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	int current = tabWidget->count();
	QWidget* tab = new QWidget{tabWidget};
	tabWidget->addTab(tab, "Exercise 1.2");

	// Layouts
	auto mainLayout = new QVBoxLayout{tab};
	auto controlLayout = new QVBoxLayout{};
	auto controlGroup = new QGroupBox{"Simulation Controls"};
	controlGroup->setLayout(controlLayout);

	// Number of bodies (minimum 3)
	auto nLayout = new QHBoxLayout{};
	auto nLabel = new QLabel{"Number of bodies (n ≥ 3)"};
	auto nSpin = new QSpinBox{};
	nSpin->setRange(3, 20);
	nSpin->setValue(3);
	nLayout->addWidget(nLabel);
	nLayout->addWidget(nSpin);
	controlLayout->addLayout(nLayout);

	// Per-body controls: position, velocity, mass
	std::vector<QDoubleSpinBox*> rSpinX, rSpinY, rSpinZ, vSpinX, vSpinY, vSpinZ, mSpin;
	auto bodiesGroup = new QGroupBox{"Bodies"};
	auto bodiesLayout = new QVBoxLayout{};
	bodiesGroup->setLayout(bodiesLayout);
	std::vector<QGroupBox *> bodyGroups;

	for (int i = 0; i < 20; ++i) {
		// Gruppenbox für jeden Körper
		auto bodyGroup = new QGroupBox{QString("Body %1").arg(i + 1)};
		auto bodyVBox = new QVBoxLayout{};

		// Position
		auto rLayout = new QHBoxLayout{};
		auto rLabel = new QLabel{"r₀:"};
		auto rx = new QDoubleSpinBox{}; rx->setRange(-10000, 10000); rx->setDecimals(3); rx->setValue(i == 0 ? 0.0 : (i == 1 ? 0.8660254 : -0.8660254));
		auto ry = new QDoubleSpinBox{}; ry->setRange(-10000, 10000); ry->setDecimals(3); ry->setValue(i == 0 ? 1.0 : (i == 1 ? -0.5 : -0.5));
		auto rz = new QDoubleSpinBox{}; rz->setRange(-10000, 10000); rz->setDecimals(3); rz->setValue(i == 0 ? 0.0 : (i == 1 ? 0.0 : 0.0));
		rLayout->addWidget(rLabel);
		rLayout->addWidget(rx); rLayout->addWidget(ry); rLayout->addWidget(rz);
		bodyVBox->addLayout(rLayout);

		// Velocity
		auto vLayout = new QHBoxLayout{};
		auto vLabel = new QLabel{"v₀:"};
		auto vx = new QDoubleSpinBox{}; vx->setRange(-10000, 10000); vx->setDecimals(3); vx->setValue(i == 0 ? -1.224745 : (i == 1 ? 0.612372 : 0.612372));
		auto vy = new QDoubleSpinBox{}; vy->setRange(-10000, 10000); vy->setDecimals(3); vy->setValue(i == 0 ? 0.0 : (i == 1 ? 1.06066 : -1.06066));
		auto vz = new QDoubleSpinBox{}; vz->setRange(-10000, 10000); vz->setDecimals(3); vz->setValue(i == 0 ? 0.0 : (i == 1 ? 0.0 : 0.0));
		vLayout->addWidget(vLabel);
		vLayout->addWidget(vx); vLayout->addWidget(vy); vLayout->addWidget(vz);
		bodyVBox->addLayout(vLayout);

		// Mass
		auto mLayout = new QHBoxLayout{};
		auto mLabel = new QLabel{"m:"};
		auto mass = new QDoubleSpinBox{}; mass->setRange(0.0001, 1000000.0); mass->setDecimals(6); mass->setValue(1.0);
		mLayout->addWidget(mLabel);
		mLayout->addWidget(mass);
		bodyVBox->addLayout(mLayout);

		// hide body group if index exceeds nSpin value
		if(i >= nSpin->value())
		{
			bodyGroup->setHidden(true);
		}
		bodyGroups.push_back(bodyGroup);

		bodyGroup->setLayout(bodyVBox);
		bodiesLayout->addWidget(bodyGroup);

		rSpinX.push_back(rx); rSpinY.push_back(ry); rSpinZ.push_back(rz);
		vSpinX.push_back(vx); vSpinY.push_back(vy); vSpinZ.push_back(vz);
		mSpin.push_back(mass);
	}
	controlLayout->addWidget(bodiesGroup);

	// G (Gravitational constant)
	auto gLayout = new QHBoxLayout{};
	auto gLabel = new QLabel{"G (Gravitational constant)"};
	auto gSpin = new QDoubleSpinBox{};
	gSpin->setRange(0.0001, 1000000.0);
	gSpin->setDecimals(6);
	gSpin->setValue(1.0);
	gLayout->addWidget(gLabel);
	gLayout->addWidget(gSpin);
	controlLayout->addLayout(gLayout);

	// delta t (Zeitschritt)
	auto dtLayout = new QHBoxLayout{};
	auto dtLabel = new QLabel{"Δt (Zeitschritt)"};
	auto dtSpin = new QDoubleSpinBox{};
	dtSpin->setRange(0.0001, 100.0);
	dtSpin->setDecimals(5);
	dtSpin->setValue(0.01);
	dtLayout->addWidget(dtLabel);
	dtLayout->addWidget(dtSpin);
	controlLayout->addLayout(dtLayout);

	// Reset-Button
	auto resetButton = new QPushButton{"Reset"};
	controlLayout->addWidget(resetButton);

	auto glWidget = new QWidget{};
	mainLayout->addWidget(controlGroup);
	mainLayout->addWidget(glWidget, 1);

	QWidget::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index != current)
			return;

		mainWindow->setRendererFactory(
			[=] (QObject * parent) {
				auto sim = new Prog1_2Simulation{};
				sim->setParent(parent);

				QObject::connect(resetButton, &QPushButton::clicked, [=]() {
					std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> bodies;
					for (int i = 0; i < nSpin->value(); ++i) {
						Eigen::Vector3d r(rSpinX[i]->value(), rSpinY[i]->value(), rSpinZ[i]->value());
						Eigen::Vector3d v(vSpinX[i]->value(), vSpinY[i]->value(), vSpinZ[i]->value());
						double m = mSpin[i]->value();
						bodies.emplace_back(r, v, m);
					}
					sim->reset(
						bodies,
						gSpin->value(),
						dtSpin->value()
					);
				});
				return sim;
			}
		);
	});

	QObject::connect(nSpin, QOverload<int>::of(&QSpinBox::valueChanged), [=] (int n) mutable {
		for(int i = 0; i < 20; ++i)
		{
			bodyGroups[i]->setVisible(i < n);
		}
	});
}

void addEx1_3Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto tab = new QWidget{tabWidget};
	tabWidget->addTab(tab, "Exercise 1.3");

	// Layouts
	auto mainLayout = new QVBoxLayout{tab};
	auto controlLayout = new QVBoxLayout{};
	auto controlGroup = new QGroupBox{"Simulation Controls"};
	controlGroup->setLayout(controlLayout);

	// scale 
	auto scaleLayout = new QHBoxLayout{};
	auto scaleLabel = new QLabel{"Scale"};
	auto scaleSpin = new QDoubleSpinBox{};
	scaleSpin->setRange(0.0001, 100.0);
	scaleSpin->setDecimals(5);
	scaleSpin->setValue(5.00);
	scaleLayout->addWidget(scaleLabel);
	scaleLayout->addWidget(scaleSpin);
	controlLayout->addLayout(scaleLayout);

	// a (Tempteraturleitfähigkeit)
	auto aLayout = new QHBoxLayout{};
	auto aLabel = new QLabel{"a (Tempteraturleitfähigkeit)"};
	auto aSpin = new QDoubleSpinBox{};
	aSpin->setRange(0.0001, 1.0);
	aSpin->setDecimals(5);
	aSpin->setValue(0.1);
	aLayout->addWidget(aLabel);
	aLayout->addWidget(aSpin);
	controlLayout->addLayout(aLayout);

	// delta t (Zeitschritt)
	auto dtLayout = new QHBoxLayout{};
	auto dtLabel = new QLabel{"Δt (Zeitschritt)"};
	auto dtSpin = new QDoubleSpinBox{};
	dtSpin->setRange(0.0001, 100.0);
	dtSpin->setDecimals(5);
	dtSpin->setValue(0.1);
	dtLayout->addWidget(dtLabel);
	dtLayout->addWidget(dtSpin);
	controlLayout->addLayout(dtLayout);

	// Reset-Button
	auto resetButton = new QPushButton{"Reset"};
	controlLayout->addWidget(resetButton);

	auto glWidget = new QWidget{}; // Platzhalter
	// -> Das tatsächliche Widget wird über das RendererFactory erstellt und dann hier gesetzt

	mainLayout->addWidget(controlGroup);
	mainLayout->addWidget(glWidget, 1);
	QWidget::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index != current)
			return;

		mainWindow->setRendererFactory(
			[=] (QObject * parent) {
			auto sim = new Prog1_3Simulation{};
			sim->setParent(parent);

			QObject::connect(resetButton, &QPushButton::clicked, [=] () {
				sim->reset(
					dtSpin->value(),
					aSpin->value()
				);
			});

			QObject::connect(scaleSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), [=] (double scale) mutable {
				sim->setTemperatureScale(
					scale
				);
			});

			return sim;
		}
		);
	});

}

void addExercise1Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto exerciseTab = new QWidget{tabWidget};
	tabWidget->addTab(exerciseTab, "Exercise 1");

	// Unter-TabWidget für die einzelnen Exercises
	auto subTabWidget = new QTabWidget{exerciseTab};
	subTabWidget->setTabPosition(QTabWidget::TabPosition::North);

	// Layout für den Übung 1 Tab
	auto layout = new QVBoxLayout{exerciseTab};
	layout->addWidget(subTabWidget);

	// Die bisherigen Exercise-Tabs als Unter-Tabs hinzufügen
	addEx1_1Tab(mainWindow, subTabWidget);
	addEx1_2Tab(mainWindow, subTabWidget);
	addEx1_3Tab(mainWindow, subTabWidget);

	subTabWidget->setCurrentIndex(0);

	// Beim Wechsel auf diesen Haupt-Tab das aktuelle Sub-Tab triggern
	QObject::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index == current)
		{
			emit subTabWidget->currentChanged(subTabWidget->currentIndex());
		}
	});
}

QWidget * createSpringBlock(const QString & springType, QDoubleSpinBox *& ksSpin, QDoubleSpinBox *& kdSpin)
{
	auto block = new QWidget;
	auto layout = new QVBoxLayout{block};

	// Feder-Art oben
	auto typeLabel = new QLabel{springType};
	layout->addWidget(typeLabel);

	// ks und kd unten
	auto paramsLayout = new QHBoxLayout{};
	auto ksLabel = new QLabel{"ks:"};
	ksSpin = new QDoubleSpinBox{};
	ksSpin->setRange(-1000.0, 1000.0);
	ksSpin->setDecimals(4);
	ksSpin->setValue(10.0);

	auto kdLabel = new QLabel{"kd:"};
	kdSpin = new QDoubleSpinBox{};
	kdSpin->setRange(-1000.0, 1000.0);
	kdSpin->setDecimals(4);
	kdSpin->setValue(0.1);

	paramsLayout->addWidget(ksLabel);
	paramsLayout->addWidget(ksSpin);
	paramsLayout->addWidget(kdLabel);
	paramsLayout->addWidget(kdSpin);

	layout->addLayout(paramsLayout);

	return block;
}
QWidget * createSpringBlock(const QString & springType, QDoubleSpinBox *& ksSpin, QDoubleSpinBox *& kdSpin, double ks, double kd)
{
	auto block = new QWidget;
	auto layout = new QVBoxLayout{block};

	// Feder-Art oben
	auto typeLabel = new QLabel{springType};
	layout->addWidget(typeLabel);

	// ks und kd unten
	auto paramsLayout = new QHBoxLayout{};
	auto ksLabel = new QLabel{"ks:"};
	ksSpin = new QDoubleSpinBox{};
	ksSpin->setRange(0.0, 10000.0);
	ksSpin->setDecimals(4);
	ksSpin->setValue(ks);

	auto kdLabel = new QLabel{"kd:"};
	kdSpin = new QDoubleSpinBox{};
	kdSpin->setRange(0.0, 1000.0);
	kdSpin->setDecimals(4);
	kdSpin->setValue(kd);

	paramsLayout->addWidget(ksLabel);
	paramsLayout->addWidget(ksSpin);
	paramsLayout->addWidget(kdLabel);
	paramsLayout->addWidget(kdSpin);

	layout->addLayout(paramsLayout);

	return block;
}

void addEx2_1Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto tab = new QWidget{tabWidget};
	tabWidget->addTab(tab, "Exercise 2.1 a)");

	// Layouts
	auto mainLayout = new QVBoxLayout{tab};
	auto controlLayout = new QVBoxLayout{};
	auto controlGroup = new QGroupBox{"Simulation Controls"};
	controlGroup->setLayout(controlLayout);

	// delta t (Zeitschritt)
	auto dtLayout = new QHBoxLayout{};
	auto dtLabel = new QLabel{"Δt (Zeitschritt)"};
	auto dtSpin = new QDoubleSpinBox{};
	dtSpin->setRange(0.0001, 100.0);
	dtSpin->setDecimals(5);
	dtSpin->setValue(0.001);
	dtLayout->addWidget(dtLabel);
	dtLayout->addWidget(dtSpin);
	controlLayout->addLayout(dtLayout);


	auto boundaryConditionLayout = new QHBoxLayout{};
	auto boundaryConditionLabel = new QLabel{"Randbedingung"};
	auto boundaryConditionCombo = new QComboBox{};
	boundaryConditionCombo->addItem("Eine Ecke", QVariant::fromValue(0));
	boundaryConditionCombo->addItem("Zwei Ecken", QVariant::fromValue(1));
	boundaryConditionCombo->addItem("Diagonale Ecken", QVariant::fromValue(2));
	boundaryConditionCombo->addItem("Drei Ecken", QVariant::fromValue(3));
	boundaryConditionCombo->addItem("Alle Ecken", QVariant::fromValue(4));
	boundaryConditionLayout->addWidget(boundaryConditionLabel);
	boundaryConditionLayout->addWidget(boundaryConditionCombo);
	controlLayout->addLayout(boundaryConditionLayout);


	QDoubleSpinBox * structKs, * structKd;
	auto structSpringBlock = createSpringBlock("Struktur-Feder", structKs, structKd, 100.0, 0.5);
	controlLayout->addWidget(structSpringBlock);

	QDoubleSpinBox * shearKs, * shearKd;
	auto shearSpringBlock = createSpringBlock("Scher-Feder", shearKs, shearKd, 50.0, 0.25);
	controlLayout->addWidget(shearSpringBlock);

	QDoubleSpinBox * bendKs, * bendKd;
	auto bendSpringBlock = createSpringBlock("Biege-Feder", bendKs, bendKd, 10.0, 0.05);
	controlLayout->addWidget(bendSpringBlock);


	// Reset-Button
	auto resetButton = new QPushButton{"Reset"};
	controlLayout->addWidget(resetButton);

	auto glWidget = new QWidget{}; // Platzhalter
	// -> Das tatsächliche Widget wird über das RendererFactory erstellt und dann hier gesetzt

	mainLayout->addWidget(controlGroup);
	mainLayout->addWidget(glWidget, 1);
	QWidget::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index != current)
			return;

		mainWindow->setRendererFactory(
			[=] (QObject * parent) {
			auto sim = new Prog2_1Simulation{};
			sim->setParent(parent);

			QObject::connect(resetButton, &QPushButton::clicked, [=] () {
				sim->reset(
					dtSpin->value(),
					static_cast<Prog2_1Simulation::BoundaryCondition>(boundaryConditionCombo->currentData().value<int>()),
					structKs->value(),
					structKd->value(),
					shearKs->value(),
					shearKd->value(),
					bendKs->value(),
					bendKd->value()
				);
			});

			return sim;
		}
		);
	});

}

void addEx2_2Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto tab = new QWidget{tabWidget};
	tabWidget->addTab(tab, "Exercise 2.1 b)");

	// Layouts
	auto mainLayout = new QVBoxLayout{tab};
	auto controlLayout = new QVBoxLayout{};
	auto controlGroup = new QGroupBox{"Simulation Controls"};
	controlGroup->setLayout(controlLayout);

	// delta t (Zeitschritt)
	auto dtLayout = new QHBoxLayout{};
	auto dtLabel = new QLabel{"Δt (Zeitschritt)"};
	auto dtSpin = new QDoubleSpinBox{};
	dtSpin->setRange(0.0001, 100.0);
	dtSpin->setDecimals(5);
	dtSpin->setValue(0.001);
	dtLayout->addWidget(dtLabel);
	dtLayout->addWidget(dtSpin);
	controlLayout->addLayout(dtLayout);


	auto boundaryConditionLayout = new QHBoxLayout{};
	auto boundaryConditionLabel = new QLabel{"Randbedingung"};
	auto boundaryConditionCombo = new QComboBox{};
	boundaryConditionCombo->addItem("Eine Ecke", QVariant::fromValue(0));
	boundaryConditionCombo->addItem("Zwei Ecken", QVariant::fromValue(1));
	boundaryConditionCombo->addItem("Diagonale Ecken", QVariant::fromValue(2));
	boundaryConditionCombo->addItem("Drei Ecken", QVariant::fromValue(3));
	boundaryConditionCombo->addItem("Alle Ecken", QVariant::fromValue(4));
	boundaryConditionLayout->addWidget(boundaryConditionLabel);
	boundaryConditionLayout->addWidget(boundaryConditionCombo);
	controlLayout->addLayout(boundaryConditionLayout);


	QDoubleSpinBox * structKs, * structKd;
	auto structSpringBlock = createSpringBlock("Struktur-Feder", structKs, structKd, 100, 0.5);
	controlLayout->addWidget(structSpringBlock);

	QDoubleSpinBox * shearKs, * shearKd;
	auto shearSpringBlock = createSpringBlock("Scher-Feder", shearKs, shearKd, 50, 0.25);
	controlLayout->addWidget(shearSpringBlock);

	QDoubleSpinBox * bendKs, * bendKd;
	auto bendSpringBlock = createSpringBlock("Biege-Feder", bendKs, bendKd, 10, 0.05);
	controlLayout->addWidget(bendSpringBlock);


	// Reset-Button
	auto resetButton = new QPushButton{"Reset"};
	controlLayout->addWidget(resetButton);

	auto glWidget = new QWidget{}; // Platzhalter
	// -> Das tatsächliche Widget wird über das RendererFactory erstellt und dann hier gesetzt

	mainLayout->addWidget(controlGroup);
	mainLayout->addWidget(glWidget, 1);
	QWidget::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index != current)
			return;

		mainWindow->setRendererFactory(
			[=] (QObject * parent) {
			auto sim = new Prog2_2Simulation{};
			sim->setParent(parent);

			QObject::connect(resetButton, &QPushButton::clicked, [=] () {
				sim->reset(
					dtSpin->value(),
					static_cast<Prog2_2Simulation::BoundaryCondition>(boundaryConditionCombo->currentData().value<int>()),
					structKs->value(),
					structKd->value(),
					shearKs->value(),
					shearKd->value(),
					bendKs->value(),
					bendKd->value()
				);
			});

			return sim;
		}
		);
	});

}

void addEx2_3Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto tab = new QWidget{tabWidget};
	tabWidget->addTab(tab, "Exercise 2.1 c)");

	// Layouts
	auto mainLayout = new QVBoxLayout{tab};
	auto controlLayout = new QVBoxLayout{};
	auto controlGroup = new QGroupBox{"Simulation Controls"};
	controlGroup->setLayout(controlLayout);

	// delta t (Zeitschritt)
	auto dtLayout = new QHBoxLayout{};
	auto dtLabel = new QLabel{"Δt (Zeitschritt)"};
	auto dtSpin = new QDoubleSpinBox{};
	dtSpin->setRange(0.0001, 100.0);
	dtSpin->setDecimals(5);
	dtSpin->setValue(0.001);
	dtLayout->addWidget(dtLabel);
	dtLayout->addWidget(dtSpin);
	controlLayout->addLayout(dtLayout);

	// max Simulation Steps
	auto mssLayout = new QHBoxLayout{};
	auto mssLabel = new QLabel{"Max Simulation Steps"};
	auto mssSpin = new QSpinBox{};
	mssSpin->setRange(1, 1000);
	mssSpin->setSingleStep(1);
	mssSpin->setValue(10);
	mssLayout->addWidget(mssLabel);
	mssLayout->addWidget(mssSpin);
	controlLayout->addLayout(mssLayout);


	auto boundaryConditionLayout = new QHBoxLayout{};
	auto boundaryConditionLabel = new QLabel{"Randbedingung"};
	auto boundaryConditionCombo = new QComboBox{};
	boundaryConditionCombo->addItem("Eine Ecke", QVariant::fromValue(0));
	boundaryConditionCombo->addItem("Zwei Ecken", QVariant::fromValue(1));
	boundaryConditionCombo->addItem("Diagonale Ecken", QVariant::fromValue(2));
	boundaryConditionCombo->addItem("Drei Ecken", QVariant::fromValue(3));
	boundaryConditionCombo->addItem("Alle Ecken", QVariant::fromValue(4));
	boundaryConditionLayout->addWidget(boundaryConditionLabel);
	boundaryConditionLayout->addWidget(boundaryConditionCombo);
	controlLayout->addLayout(boundaryConditionLayout);


	QDoubleSpinBox * structKs, * structKd;
	auto structSpringBlock = createSpringBlock("Struktur-Feder", structKs, structKd, 100, 0.5);
	controlLayout->addWidget(structSpringBlock);

	QDoubleSpinBox * shearKs, * shearKd;
	auto shearSpringBlock = createSpringBlock("Scher-Feder", shearKs, shearKd, 50, 0.25);
	controlLayout->addWidget(shearSpringBlock);

	QDoubleSpinBox * bendKs, * bendKd;
	auto bendSpringBlock = createSpringBlock("Biege-Feder", bendKs, bendKd, 10, 0.05);
	controlLayout->addWidget(bendSpringBlock);


	// Reset-Button
	auto resetButton = new QPushButton{"Reset"};
	controlLayout->addWidget(resetButton);

	auto glWidget = new QWidget{}; // Platzhalter
	// -> Das tatsächliche Widget wird über das RendererFactory erstellt und dann hier gesetzt

	mainLayout->addWidget(controlGroup);
	mainLayout->addWidget(glWidget, 1);
	QWidget::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index != current)
			return;

		mainWindow->setRendererFactory(
			[=] (QObject * parent) {
			auto sim = new Prog2_3Simulation{};
			sim->setParent(parent);

			QObject::connect(resetButton, &QPushButton::clicked, [=] () {
				sim->reset(
					dtSpin->value(),
					static_cast<Prog2_3Simulation::BoundaryCondition>(boundaryConditionCombo->currentData().value<int>()),
					structKs->value(),
					structKd->value(),
					shearKs->value(),
					shearKd->value(),
					bendKs->value(),
					bendKd->value(),
					mssSpin->value()
				);
			});

			return sim;
		}
		);
	});

}

void addExercise2Tab(GLMainWindow * mainWindow, QTabWidget * tabWidget)
{
	auto current = tabWidget->count();
	auto exerciseTab = new QWidget{tabWidget};
	tabWidget->addTab(exerciseTab, "Exercise 2");

	// Unter-TabWidget für die einzelnen Exercises
	auto subTabWidget = new QTabWidget{exerciseTab};
	subTabWidget->setTabPosition(QTabWidget::TabPosition::North);

	// Layout für den Übung 1 Tab
	auto layout = new QVBoxLayout{exerciseTab};
	layout->addWidget(subTabWidget);

	// Die bisherigen Exercise-Tabs als Unter-Tabs hinzufügen
	addEx2_1Tab(mainWindow, subTabWidget);
	addEx2_2Tab(mainWindow, subTabWidget);
	addEx2_3Tab(mainWindow, subTabWidget);

	subTabWidget->setCurrentIndex(0);

	QObject::connect(tabWidget, &QTabWidget::currentChanged, [=] (int index) {
		if(index == current)
		{
			emit subTabWidget->currentChanged(subTabWidget->currentIndex());
		}
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
    addExercise1Tab(mainWindow, tabWidget); // Nur noch ein Tab für alle Übungen
	addExercise2Tab(mainWindow, tabWidget); // Übung 2 Tab
    tabWidget->setCurrentIndex(0);

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
