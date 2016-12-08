#include "Game.h"


bool flip = false;
int current = 1;

Game::Game() : window(VideoMode(800, 600), "OpenGL Cube"), primitives(1)
{
	index = glGenLists(primitives);
}

Game::~Game(){}

void Game::run()
{

	initialize();

	Event event;

	while (isRunning){

		cout << "Game running..." << endl;

		while (window.pollEvent(event))
		{
			if (event.type == Event::Closed)
			{
				isRunning = false;
			}
		}
		update();
		draw();
	}

}

void Game::initialize()
{
	isRunning = true;

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, window.getSize().x / window.getSize().y, 1.0, 500.0);
	glMatrixMode(GL_MODELVIEW);

	// glNewList(index, GL_COMPILE);
	// Creates a new Display List
	// Initalizes and Compiled to GPU
	// https://www.opengl.org/sdk/docs/man2/xhtml/glNewList.xml
	glNewList(index, GL_COMPILE);
	glBegin(GL_QUADS);
	{
		//Front Face
		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(1.0f, 1.0f, -5.0f);
		glVertex3f(-1.0f, 1.0f, -5.0f);
		glVertex3f(-1.0f, -1.0f, -5.0f);
		glVertex3f(1.0f, -1.0f, -5.0f);

		//Back Face
		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(1.0f, 1.0f, -15.0f);
		glVertex3f(-1.0f, 1.0f, -15.0f);
		glVertex3f(-1.0f, -1.0f, -15.0f);
		glVertex3f(1.0f, -1.0f, -15.0f);

		//Complete the faces of the Cube

		//left face
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(-1.0f, 1.0f, -5.0f);
		glVertex3f(-1.0f, 1.0f, -15.0f);
		glVertex3f(-1.0f, -1.0f, -15.0f);
		glVertex3f(-1.0f, -1.0f, -5.0f);
		
		//right face
		glColor3f(1.0, 1.0, 0.0);
		glVertex3f(1.0f, 1.0f, -5.0f);
		glVertex3f(1.0f, 1.0f, -15.0f);
		glVertex3f(1.0f, -1.0f, -15.0f);
		glVertex3f(1.0f, -1.0f, -5.0f);
		
		//top face
		glColor3f(0.0, 1.0, 1.0);
		glVertex3f(1.0f, 1.0f, -15.0f);
		glVertex3f(-1.0f, 1.0f, -15.0f);
		glVertex3f(-1.0f, 1.0f, -5.0f);
		glVertex3f(1.0f, 1.0f, -5.0f);
		
		//bottom face
		glColor3f(1.0, 0.0, 1.0);
		glVertex3f(1.0f, -1.0f, -15.0f);
		glVertex3f(-1.0f, -1.0f, -15.0f);
		glVertex3f(-1.0f, -1.0f, -5.0f);
		glVertex3f(1.0f, -1.0f, -5.0f);
	}
	glEnd();
	glEndList();
}

void Game::update()
{
	elapsed = clock.getElapsedTime();

	if (elapsed.asSeconds() >= 1.0f)
	{
		clock.restart();

		if (!flip)
		{
			flip = true;
			current++;
			if (current > primitives)
			{
				current = 1;
			}
		}
		else
			flip = false;
	}

	if (flip)
	{
		rotationAngle += 0.005f;

		if (rotationAngle > 360.0f)
		{
			rotationAngle -= 360.0f;
		}
	}
	
	cout << "Update up" << endl;
}

void Game::draw()
{
	cout << "Drawing" << endl;
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	cout << "Drawing Cube " << current << endl;
	glLoadIdentity();
	

	//glRotatef(rotationAngle, 0, 1, 0); // Rotates the camera on Y Axis

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glCallList(current);

	window.display();

}

void Game::unload()
{
	cout << "Cleaning up" << endl;
}