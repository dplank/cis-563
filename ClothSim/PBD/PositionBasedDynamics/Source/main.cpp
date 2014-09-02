#include <iostream>
#include "openGL_headers.h"
#include "math_headers.h"
#include "cloth_sim.h"
#include "scene.h"
#include "stb_image_write.h"

int window_width = 1024;
int window_height = 768;

//----------State Control----------//
bool pause = false;
bool record = false;

//----------Time----------//
double now, lastTime;
float delta_t = 0.0f;
int frame_num = 0;

//----------Mouse Control----------//
int mouse_old_x, mouse_old_y;
unsigned char button_mask = 0x00;

//----------OpenGL Render Control----------//
GLuint m_uniform_location[2];
GLuint m_vert_handle, m_frag_handle, m_shaderprog_handle;

//----------Camera Control----------//
GLuint vbo_handle[4];
float eye_distance = 20.0f;
float head = 45.0f, pitch = 45.0f;
glm::vec3 cam_pos, up(0.0f, 1.0f, 0.0f), lookat(0.0f, 4.0f, 0.0f);

//----------Simulation----------//
Scene scene("../Scene/test_scene.xml");
// TODO: change here if you want to use a smaller iteration number.
ClothSim cloth_sim(10);

//----------functions----------//
// declare
bool initGL(void);
void initShader(const char* vert_path, const char* frag_path);
void cleanupShader();
void mouseClick(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void keypress(unsigned char key, int x, int y);
void display(void);
void aimCamera(void);
void drawFps(void);
void drawAxes(void);
void grabScreen(void);

void activate_shaderprog(GLuint shaderprog);
void deactivate_shaderprog(GLuint shaderprog);

// courtesy of Swiftless
char* textFileRead(const char* fileName);
void printLinkInfoLog(int prog);
void printShaderInfoLog(int shader);

// define
bool initGL(void)
{
    glewInit();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glViewport(0, 0, window_width, window_height);

    // projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // TODO (maybe) :: depending on your parameters, you may need to change
    // near and far view distances (1, 500), to better see the simulation.
    // If you do this, probably also change translate_z initial value at top.
    gluPerspective(60.0, (GLfloat)window_width / (GLfloat) window_height, 0.1, 100);
    return true;
}

void initShader(const char* vert_path, const char* frag_path)
{
    // generate vbos
    glGenBuffers(4, vbo_handle);

    // create shaders and shader program
    m_vert_handle = glCreateShader(GL_VERTEX_SHADER);
    m_frag_handle = glCreateShader(GL_FRAGMENT_SHADER);
    m_shaderprog_handle = glCreateProgram();

    // load shader source from file
    const char* vert_source = textFileRead(vert_path);
    const char* frag_source = textFileRead(frag_path);
    glShaderSource(m_vert_handle, 1, &vert_source, NULL);
    glShaderSource(m_frag_handle, 1, &frag_source, NULL);
    glCompileShader(m_vert_handle);
    glCompileShader(m_frag_handle);

    // compile shader source
    GLint compiled;
    glGetShaderiv(m_vert_handle, GL_COMPILE_STATUS, &compiled);
    if(!compiled)
        printShaderInfoLog(m_vert_handle);
    glGetShaderiv(m_frag_handle, GL_COMPILE_STATUS, &compiled);
    if(!compiled)
        printShaderInfoLog(m_frag_handle);

    // TODO: customize this part if you want to modify the glsl shader.
    // bind attribute locations for the shaders
    // 0 for position, 1 for color, 2 for normal.
    glBindAttribLocation(m_shaderprog_handle, 0, "v_position");
    glBindAttribLocation(m_shaderprog_handle, 1, "v_color");
    glBindAttribLocation(m_shaderprog_handle, 2, "v_normal");

    // attach shader to the shader program
    glAttachShader(m_shaderprog_handle, m_vert_handle);
    glAttachShader(m_shaderprog_handle, m_frag_handle);
    glLinkProgram(m_shaderprog_handle);
    GLint linked;
    glGetProgramiv(m_shaderprog_handle, GL_LINK_STATUS, &linked);
    if(!linked)
        printLinkInfoLog(m_shaderprog_handle);

    // TODO: customize this part if you want to modify the glsl shader.
    // query uniform locations from openGL.
    m_uniform_location[0] = glGetUniformLocation(m_shaderprog_handle, "u_modelviewMatrix");
    m_uniform_location[1] = glGetUniformLocation(m_shaderprog_handle, "u_projMatrix");

    // activate the shader program.
    glUseProgram(m_shaderprog_handle);
}

void cleanupShader()
{
    glDeleteBuffers(4, vbo_handle);
    glDetachShader(m_shaderprog_handle, m_vert_handle);
    glDetachShader(m_shaderprog_handle, m_frag_handle);
    glDeleteShader(m_vert_handle);
    glDeleteShader(m_frag_handle);
    glDeleteProgram(m_shaderprog_handle);
}

void mouseClick(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN) {
        button_mask |= 0x01 << button;
    } else if (state == GLUT_UP) {
        unsigned char mask_not = ~button_mask;
        mask_not |= 0x01 << button;
        button_mask = ~mask_not;
    }

    mouse_old_x = x;
    mouse_old_y = y;
    //glutPostRedisplay();
}

void mouseMotion(int x, int y)
{
    float dx, dy;
    dx = (float)(x - mouse_old_x);
    dy = (float)(y - mouse_old_y);

    if (button_mask & 0x01) 
    {// left button
        head += dy * 0.2f;
        pitch += dx * 0.2f;
    } 
    else if (button_mask & 0x02) 
    {// right button
        eye_distance -= dy * 0.01f;
    }
    else if (button_mask & 0x04)
    {// middle button
        glm::vec3 vdir(lookat - cam_pos);
        glm::vec3 u(glm::normalize(glm::cross(vdir, up)));
        glm::vec3 v(glm::normalize(glm::cross(u, vdir)));

        lookat += 0.01f * (dy * v - dx * u);
    }

    mouse_old_x = x;
    mouse_old_y = y;
}

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    //operations to reset camera's projection matrix.
    glm::mat4 projection = glm::perspective(60.0f, static_cast<float>(width) / static_cast<float>(height), 0.1f, 100.0f);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(&projection[0][0]);

    activate_shaderprog(m_shaderprog_handle);
    glUniformMatrix4fv(m_uniform_location[1], 1, false, &projection[0][0]);
    glutPostRedisplay();
}

void keypress(unsigned char key, int x, int y)
{
    switch(key)
    {
    case 'q':
    case 'Q':
        cleanupShader();
        exit(0);
        break;
    case 'r':
    case 'R':
        record = !record;
        break;
    case 'f':
    case 'F':
        cloth_sim.flip_draw_mode();
        break;
    case ' ':
        pause = !pause;
        break;
    }
    glutPostRedisplay();
}

void drawFps(void)
{
    glColor4f(1.0, 1.0, 1.0, 1.0);
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    char info[1024];
    glRasterPos2f(0.00f, 0.98f);
    sprintf(info, "Total Frame: %u", 
        frame_num);
    for (unsigned int i = 0; i < strlen(info); i++)
    {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
    }

    glRasterPos2f(0.00f, 0.96f);
    sprintf(info, "Frame Rate: %3.1f", 
        1.0f / delta_t);
    for (unsigned int i = 0; i < strlen(info); i++)
    {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
    }

    glRasterPos2f(0.00f, 0.94f);
    sprintf(info, "%s", 
        record ? "Recording" : "");
    for (unsigned int i = 0; i < strlen(info); i++)
    {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
    }

    glPopAttrib();
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(lastTime == 0)
    {
        lastTime = glutGet(GLUT_ELAPSED_TIME);
    }
    now = glutGet(GLUT_ELAPSED_TIME);

    delta_t = (now - lastTime) / 1000.0f;
    lastTime = now;

    aimCamera();

    if(!pause)
        cloth_sim.update(&scene, 0.006325f);

    activate_shaderprog(m_shaderprog_handle);
    cloth_sim.draw(vbo_handle);
    scene.draw(vbo_handle);
    deactivate_shaderprog(m_shaderprog_handle);

    drawAxes();
    drawFps();
    if(!pause && record)
        grabScreen();
    frame_num++;
    glutPostRedisplay();
    glutSwapBuffers();	
}

void drawAxes(void)
{
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    //draw axis.
    GLfloat modelview[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, modelview);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, &viewport[0]);
    GLint width = viewport[2] / 16;
    GLint height = viewport[3] / 16;
    glViewport(0, 0, width, height);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    //get the camera position and up vector from the modelview matrix.
    double campos[3] = {0.0 + 2.0f * modelview[2], 0.0 + 2.0f * modelview[6], 0.0 + 2.0f * modelview[10]};
    double up[3] = {modelview[1], modelview[5], modelview[9]};
    //set up the view matrix.
    gluLookAt(campos[0], campos[1], campos[2], 
        0.0, 0.0, 0.0,
        up[0], up[1], up[2]);

    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);

    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();

    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    glPopMatrix();
    glPopAttrib();
}

void aimCamera(void)
{
    float r_head = glm::radians(head), r_pitch = glm::radians(pitch);
    cam_pos.x = lookat.x + eye_distance * glm::cos(r_head) * glm::cos(r_pitch);
    cam_pos.y = lookat.y + eye_distance * glm::sin(r_head);
    cam_pos.z = lookat.z + eye_distance * glm::cos(r_head) * glm::sin(r_pitch);

    glMatrixMode(GL_MODELVIEW);
    up = glm::vec3(0.0f, (glm::cos(r_head) > 0.0f) ? 1.0f : -1.0f, 0.0f);
    glm::mat4 modelview = glm::lookAt(cam_pos, lookat, up);
    glLoadMatrixf(&modelview[0][0]);
    
    glMatrixMode(GL_PROJECTION);
    glm::mat4 projection = glm::perspective(60.0f, static_cast<float>(window_width) / static_cast<float>(window_height), 0.1f, 100.0f);
    glLoadMatrixf(&projection[0][0]);

    activate_shaderprog(m_shaderprog_handle);
    glUniformMatrix4fv(m_uniform_location[0], 1, false, &modelview[0][0]);
    glUniformMatrix4fv(m_uniform_location[1], 1, false, &projection[0][0]);
}

void grabScreen(void)
{
    unsigned char* bitmapData = new unsigned char[3 * window_width * window_height];

    for (int i=0; i < window_height; i++) 
    {
        glReadPixels(0, i, window_width, 1, GL_RGB, GL_UNSIGNED_BYTE, 
            bitmapData + (window_width * 3 * ((window_height - 1) - i)));
    }

    char anim_filename[2048];
    sprintf_s(anim_filename, 2048, "output/PBD_no_self_intersect_%04d.png", frame_num);

    stbi_write_png(anim_filename, window_width, window_height, 3, bitmapData, window_width * 3);

    delete [] bitmapData;
}

void activate_shaderprog(GLuint shaderprog)
{
    GLint current_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &current_prog);
    if(current_prog != (GLint)shaderprog)
        glUseProgram(shaderprog);
}

void deactivate_shaderprog(GLuint shaderprog)
{
    GLint current_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &current_prog);
    if(current_prog == (GLint)shaderprog)
        glUseProgram(0);
}

// global scope
int g_max_fps = 60;
int g_timestep = 1000 / g_max_fps;

void timeout(int);
void timeout(int value)
{
	glutTimerFunc(g_timestep, timeout, 0);
	glutPostRedisplay();
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("Position Based Dynamics");
    
    if(initGL() == false)
        return 0;

    initShader("./Shader/vert.glsl", "./Shader/frag.glsl");

    // TODO: change here if you want to modify the dimension.
    cloth_sim.initialize(30, 30, glm::vec3(-5.0f, 10.0f, -5.0f), glm::vec3(5.0f, 10.0f, 5.0f));

    glutDisplayFunc(display);
    glutKeyboardFunc(keypress);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMotion);
    glutReshapeFunc(reshape);

	//before glutMainLoop();
	glutTimerFunc(g_timestep, timeout, 0);

    glutMainLoop();
    return 0;
}

// helper function to read shader source and put it in a char array
// thanks to Swiftless
char* textFileRead(const char* fileName) 
{
    char* text;

    if (fileName != NULL) {
        FILE *file = fopen(fileName, "rt");

        if (file != NULL) {
            fseek(file, 0, SEEK_END);
            int count = ftell(file);
            rewind(file);

            if (count > 0) {
                text = (char*)malloc(sizeof(char) * (count + 1));
                count = fread(text, sizeof(char), count, file);
                text[count] = '\0';	//cap off the string with a terminal symbol, fixed by Cory
            }
            fclose(file);
        }
    }
    return text;
}

void printLinkInfoLog(int prog) 
{
    int infoLogLen = 0;
    int charsWritten = 0;
    GLchar *infoLog;

    glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &infoLogLen);

    // should additionally check for OpenGL errors here

    if (infoLogLen > 0)
    {
        infoLog = new GLchar[infoLogLen];
        // error check for fail to allocate memory omitted
        glGetProgramInfoLog(prog,infoLogLen, &charsWritten, infoLog);
        std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
        delete [] infoLog;
    }
}

void printShaderInfoLog(int shader)
{
    int infoLogLen = 0;
    int charsWritten = 0;
    GLchar *infoLog;

    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

    // should additionally check for OpenGL errors here

    if (infoLogLen > 0)
    {
        infoLog = new GLchar[infoLogLen];
        // error check for fail to allocate memory omitted
        glGetShaderInfoLog(shader,infoLogLen, &charsWritten, infoLog);
        std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
        delete [] infoLog;
    }

    // should additionally check for OpenGL errors here
}