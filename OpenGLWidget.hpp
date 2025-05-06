#pragma once

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#include <glad/gl.h>

#include <QOpenGLWidget>

#include <functional>

class QOpenGLDebugLogger;
class OpenGLRenderer;

class OpenGLWidget : public QOpenGLWidget
{
	Q_OBJECT

public:
	OpenGLWidget(QWidget * parent = nullptr, Qt::WindowFlags f = Qt::WindowFlags());

	void setRendererFactory(std::function<OpenGLRenderer * (QObject * parent)> rendererFactory);

	bool event(QEvent * e) override;

public slots:
	void setLoggingEnabled(bool enable);
	void setLoggingSynchronous(bool synchronous);

signals:
	void loggingEnabledChanged(bool enable);
	void loggingSynchronousChanged(bool synchronous);

protected:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int w, int h) override;

private:
	QOpenGLDebugLogger * logger = nullptr;

	std::function<OpenGLRenderer * (QObject * parent)> rendererFactory;
	OpenGLRenderer * renderer = nullptr;

	bool loggingEnabled = false, loggingSynchronous = false;
};
