#include <stdio.h>
#include <math.h>
#include <fstream>
#include "camera.h"
#include "fps.h"
#include <GL/glut.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

#include "jelloMesh.h"
#include "world.h"

JelloMesh theJello;
Camera theCamera;
//World theWorld("worlds/ground.xml");
//World theWorld("worlds/cylinders.xml");
//World theWorld("worlds/sphere.xml");
World theWorld("worlds/cube.xml");
mmc::FpsTracker theFpsTracker;

// UI Helpers
int lastX = 0, lastY = 0;
int theMenu = 0;
int theButtonState = 0;
int theModifierState = 0;
int theFrameNum = 0;
bool isRunning = true;
bool isRecording = false;

void initCamera()
{
   double w = theJello.GetWidth()*2;   
   double h = theJello.GetHeight()*2;   
   double d = theJello.GetDepth()*2;   
   double angle = 0.5*theCamera.dfltVfov*M_PI/180.0;
   double dist;
   if (w > h) dist = w*0.5/tan(angle);  // aspect is 1, so i can do this
   else dist = h*0.5/tan(angle);
   theCamera.dfltEye.set(w*0.5, h, -(dist+d));
   theCamera.dfltLook.set(0.0, 0.0, 0.0);
   theCamera.reset();
}

void grabScreen()  
{
    unsigned int image;
    ilGenImages(1, &image);
    ilBindImage(image);

    ILenum error = ilGetError();
    assert(error == IL_NO_ERROR);

    ilTexImage(640, 480, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);

    error = ilGetError();
    assert(error == IL_NO_ERROR);

    unsigned char* data = ilGetData();

    error = ilGetError();
    assert(error == IL_NO_ERROR);

    for (int i=479; i>=0; i--) 
    {
	    glReadPixels(0,i,640,1,GL_RGB, GL_UNSIGNED_BYTE, 
		    data + (640 * 3 * i));
    }

    char anim_filename[2048];
    sprintf_s(anim_filename, 2048, "output/%04d.png", theFrameNum++); 

    ilSave(IL_PNG, anim_filename);

    error = ilGetError();
    assert(error == IL_NO_ERROR);

    ilDeleteImages(1, &image);

    error = ilGetError();
    assert(error == IL_NO_ERROR);
}

void onMouseMotionCb(int x, int y)
{
   int deltaX = lastX - x;
   int deltaY = lastY - y;
   bool moveLeftRight = abs(deltaX) > abs(deltaY);
   bool moveUpDown = !moveLeftRight;

   if (theButtonState == GLUT_LEFT_BUTTON)  // Rotate
   {
      if (moveLeftRight && deltaX > 0) theCamera.orbitLeft(deltaX);
      else if (moveLeftRight && deltaX < 0) theCamera.orbitRight(-deltaX);
      else if (moveUpDown && deltaY > 0) theCamera.orbitUp(deltaY);
      else if (moveUpDown && deltaY < 0) theCamera.orbitDown(-deltaY);
   }
   else if (theButtonState == GLUT_MIDDLE_BUTTON) // Zoom
   {
       if (theModifierState & GLUT_ACTIVE_ALT) // camera move
       {
            if (moveLeftRight && deltaX > 0) theCamera.moveLeft(deltaX);
            else if (moveLeftRight && deltaX < 0) theCamera.moveRight(-deltaX);
            else if (moveUpDown && deltaY > 0) theCamera.moveUp(deltaY);
            else if (moveUpDown && deltaY < 0) theCamera.moveDown(-deltaY);
       }
       else
       {
           if (moveUpDown && deltaY > 0) theCamera.moveForward(deltaY);
           else if (moveUpDown && deltaY < 0) theCamera.moveBack(-deltaY);
       }

   }    
 
   lastX = x;
   lastY = y;
   glutPostRedisplay();
}

void onMouseCb(int button, int state, int x, int y)
{
   theButtonState = button;
   theModifierState = glutGetModifiers();
   lastX = x;
   lastY = y;
   glutSetMenu(theMenu);
}


void onKeyboardCb(unsigned char key, int x, int y)
{
   unsigned int mask = 0x0;

   if (key == ' ') theCamera.reset();
   else if (key == 27) exit(0); // ESC Key
   else if (key == '8') theJello.SetIntegrationType(JelloMesh::EULER);
   else if (key == '9') theJello.SetIntegrationType(JelloMesh::MIDPOINT);
   else if (key == '0') theJello.SetIntegrationType(JelloMesh::RK4);
   else if (key == '>') isRunning = true;
   else if (key == '=') isRunning = false;
   else if (key == '<') theJello.Reset();
   else if (key == 'r') isRecording = !isRecording; if (isRecording) theFrameNum = 0;
   else if (key == '1') mask = theJello.MESH;
   else if (key == '2') mask = theJello.FORCES;
   else if (key == '3') mask = theJello.NORMALS;
   else if (key == '4') mask = theJello.STRUCTURAL;
   else if (key == '5') mask = theJello.SHEAR;
   else if (key == '6') mask = theJello.BEND;

   if (mask)
   {
       if (theJello.GetDrawFlags() & mask)
       {
           theJello.SetDrawFlags(theJello.GetDrawFlags() & ~mask);
       }
       else
       {
           theJello.SetDrawFlags(theJello.GetDrawFlags() | mask);
       }
   }

   glutPostRedisplay();
}

void onMenuCb(int value)
{
   switch (value)
   {
   case -1: exit(0);
   default: onKeyboardCb(value, 0, 0); break;
   }
}

void onKeyboardSpecialCb(int key, int x, int y)
{
}

void onTimerCb(int value)
{
   glutTimerFunc(100, onTimerCb, 0);
   if (isRunning) 
   {
       theJello.Update(0.01, theWorld); 
       if (isRecording) grabScreen();
   }

   glutPostRedisplay();
}

void onResizeCb(int width, int height)
{
   // Update viewport
   glViewport(0, 0, width, height);

   // Update camera projection's aspect ratio
   float vfov, aspect, zNear, zFar;
   theCamera.getProjection(&vfov, &aspect, &zNear, &zFar);
   theCamera.setProjection(vfov, ((GLfloat) width)/height, zNear, zFar);
}

void drawOverlay()
{
  // Draw Overlay
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glPushAttrib(GL_LIGHTING_BIT);
     glDisable(GL_LIGHTING);

     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluOrtho2D(0.0, 1.0, 0.0, 1.0);

     glMatrixMode(GL_MODELVIEW);
     glLoadIdentity();
     glRasterPos2f(0.01, 0.01);
     
     char* intstr;
     switch (theJello.GetIntegrationType())
     {
     case JelloMesh::EULER: intstr = "Euler"; break;
     case JelloMesh::MIDPOINT: intstr = "Midpoint"; break;
     case JelloMesh::RK4: intstr = "RK4"; break;
     }

     char info[1024];
     sprintf(info, "Framerate: %3.1f %s %s", 
         theFpsTracker.fpsAverage(),
         intstr, isRecording? "(Recording ON)" : "");
 
     for (unsigned int i = 0; i < strlen(info); i++)
     {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
     }
  glPopAttrib();
}

void drawAxes()
{
  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);

      glLineWidth(2.0); 
      glBegin(GL_LINES);
         glColor3f(1.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(1.0, 0.0, 0.0);

         glColor3f(0.0, 1.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 1.0, 0.0);

         glColor3f(0.0, 0.0, 1.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 1.0);
      glEnd();
  glPopAttrib();
}

void onDrawCb()
{
    // Keep track of time
    theFpsTracker.timestamp();

    // Draw Scene and overlay
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    theCamera.draw();
    //drawAxes();

    vec3 cpos = theCamera.getPosition();
    float pos[4] = {cpos[0], cpos[1]+2.0, cpos[2],0.0};
    glLightfv(GL_LIGHT0, GL_POSITION, pos);

    theWorld.Draw();
    theJello.Draw(cpos);
    drawOverlay();
    glutSwapBuffers();
}

// Load jello parameters from the config file
int loadJelloParameters(char* filename) throw (char*)
{
	ifstream config;
	config.open(filename);
	if (config.is_open())
	{
		string line;
		while (config.good())
		{
			getline(config, line);
			char* name = strtok((char*)line.c_str(), " ");
			char* value = strtok(NULL, " ");

			if (strcmp(name, "STRUCTURAL_KS") == 0)
			{
				JelloMesh::g_structuralKs = atof(value);
			}
			else if (strcmp(name, "STRUCTURAL_KD") == 0)
			{
				JelloMesh::g_structuralKd = atof(value);
			}
			else if (strcmp(name, "SHEAR_KS") == 0)
			{
				JelloMesh::g_shearKs = atof(value);
			}
			else if (strcmp(name, "SHEAR_KD") == 0)
			{
				JelloMesh::g_shearKd = atof(value);
			}
			else if (strcmp(name, "BEND_KS") == 0)
			{
				JelloMesh::g_bendKs = atof(value);
			}
			else if (strcmp(name, "BEND_KD") == 0)
			{
				JelloMesh::g_bendKd = atof(value);
			}
			else if (strcmp(name, "PENALTY_KS") == 0)
			{
				JelloMesh::g_penaltyKs = atof(value);
			}
			else if (strcmp(name, "PENALTY_KD") == 0)
			{
				JelloMesh::g_penaltyKd = atof(value);
			}
		}
		cout << "Loaded jello parameters from " << filename << endl;
		return 0;    // normal
	} else {
		return 1;    // fail to open file
	}
}

void init(void)
{
    initCamera();
    glClearColor(0.8, 0.8, 0.8, 1.0);

    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_NORMALIZE);
    glCullFace(GL_FRONT);


    float white[4] = {1.0,1.0,1.0,1.0};
    float black[4] = {0.0,0.0,0.0,1.0};
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
    glLightfv(GL_LIGHT0, GL_SPECULAR, white);
    glLightfv(GL_LIGHT0, GL_AMBIENT, black);

    GLfloat fogColor[4]= {0.8f, 0.8f, 0.8f, 1.0f};	
    glFogi(GL_FOG_MODE, GL_LINEAR);		// Fog Mode
    glFogfv(GL_FOG_COLOR, fogColor);			// Set Fog Color
    glFogf(GL_FOG_DENSITY, 0.35f);				// How Dense Will The Fog Be
    glHint(GL_FOG_HINT, GL_DONT_CARE);			// Fog Hint Value
    glFogf(GL_FOG_START, 10.0f);				// Fog Start Depth
    glFogf(GL_FOG_END, 40.0f);				// Fog End Depth
    glEnable(GL_FOG);					// Enables GL_FOG
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		cout << "No configuration file" << endl;
		system("pause");
		return 1;
	}

    ilInit();
    iluInit();
    ilEnable(IL_FILE_OVERWRITE);
    ilutRenderer(ILUT_OPENGL);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("JelloSim by AlineN");
    glutDisplayFunc(onDrawCb);
    glutKeyboardFunc(onKeyboardCb);
    glutSpecialFunc(onKeyboardSpecialCb);
    glutMouseFunc(onMouseCb);
    glutMotionFunc(onMouseMotionCb); 
    glutTimerFunc(100, onTimerCb, 0);
    glutReshapeFunc(onResizeCb);

    int intMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Euler\t'8'", '8');
    glutAddMenuEntry("Midpoint\t'9'", '9');
    glutAddMenuEntry("RK4\t'0'", '0');

    int displayMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Mesh\t'1'", '1');
    glutAddMenuEntry("Forces\t'2'", '2');
    glutAddMenuEntry("Collision Normals\t'3'", '3');
    glutAddMenuEntry("Structural Springs\t'4'", '4');
    glutAddMenuEntry("Shear Springs\t'5'", '5');
    glutAddMenuEntry("Bend Springs\t'6'", '6');

    theMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Start\t'>'", '>');
    glutAddMenuEntry("Pause\t'='", '=');
    glutAddMenuEntry("Reset\t'<'", '<');
    glutAddMenuEntry("Record\t'r'", 'r');
    glutAddSubMenu("Integration Type", intMenu);
    glutAddSubMenu("Draw Settings", displayMenu);
    glutAddMenuEntry("_________________", -1);
    glutAddMenuEntry("Exit", 27);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

	init();
	if (loadJelloParameters(argv[1]) != 0)
	{
		cout << "Failed to open configuration file" << endl;
		system("pause");
		return 1;
	}
	theJello.Reset();

    glutMainLoop();
    return 0;             
}

