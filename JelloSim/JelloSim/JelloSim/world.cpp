#include "World.h"
#include <tinyxml.h>
#include <map>
#include <GL/glut.h>

class XMLWorldVisitor : public TiXmlVisitor
{
public:
	XMLWorldVisitor(World* world) : m_world(world), m_curBody(NULL) { }

	bool IsBody(TiXmlNode const* node)
	{
		if(!node->ToElement()) return false;
		if(node->ValueStr() == "ground") return true;
		if(node->ValueStr() == "box") return true;
		if(node->ValueStr() == "sphere") return true;
		if(node->ValueStr() == "cylinder") return true;
		return false;
	}

	/// Visit an element.
	virtual bool VisitEnter( const TiXmlElement& element, const TiXmlAttribute* attribute)
	{
		if(element.ValueStr() == "world")
		{
			if(element.Parent() != element.GetDocument()) return false;
			return true;
		}
		else if(element.ValueStr() == "materials")
		{
			return true;
		}
		else if(element.ValueStr() == "material")
		{
            return true;
		}
		else if(element.ValueStr() == "bodies")
		{
			if(element.Parent()->ValueStr() != "world") return false;
			return true;
		}
		else if(element.ValueStr() == "ground")
		{
			if(element.Parent()->ValueStr() != "bodies") return false;
			assert(m_curBody == NULL);

            m_curBody = new World::Ground();
			return true;
		}
		else if(element.ValueStr() == "box")
		{
			if(element.Parent()->ValueStr() != "bodies") return false;
			assert(m_curBody == NULL);

			double hx = 0.5;
			double hy = 0.5;
			double hz = 0.5;
			
			element.Attribute("hx", &hx);
			element.Attribute("hy", &hy);
			element.Attribute("hz", &hz);

            World::Cube* cube = new World::Cube();
            cube->hx = hx;
            cube->hy = hy;
            cube->hz = hz;
            m_curBody = cube;
			return true;
		}
		else if(element.ValueStr() == "sphere")
		{
			if(element.Parent()->ValueStr() != "bodies") return false;
			assert(m_curBody == NULL);

			double r = 1;	
			double cx, cy, cz;

			element.Attribute("r", &r);
			element.Attribute("cx", &cx);
			element.Attribute("cy", &cy);
			element.Attribute("cz", &cz);

            World::Sphere* sphere = new World::Sphere();
            sphere->r = r;
			sphere->center = vec3(cx,cy,cz);
            m_curBody = sphere;
			return true;
		}
		else if(element.ValueStr() == "cylinder")
		{
			if(element.Parent()->ValueStr() != "bodies") return false;
			assert(m_curBody == NULL);

			double r = 1;
            double startx, starty, startz;
            double endx, endy, endz;
			element.Attribute("r", &r);

            element.Attribute("sx", &startx);
			element.Attribute("sy", &starty);
			element.Attribute("sz", &startz);

            element.Attribute("ex", &endx);
			element.Attribute("ey", &endy);
			element.Attribute("ez", &endz);

            World::Cylinder* cylinder = new World::Cylinder();
            cylinder->start = vec3(startx, starty, startz);
            cylinder->end = vec3(endx, endy, endz);
            cylinder->r = r;
            m_curBody = cylinder;
			return true;
		}
		else if(element.ValueStr() == "pos")
		{
			if(!IsBody(element.Parent())) return false;
			assert(m_curBody != NULL);
			double x = 0;
			double y = 0;
			double z = 0;

			element.Attribute("x", &x);
			element.Attribute("y", &y);
			element.Attribute("z", &z);

			m_curBody->pos = vec3(float(x), float(y), float(z));
			return true;
		}
		else if(element.ValueStr() == "vel")
		{
            return false;
		}
		else if(element.ValueStr() == "ori")
		{
            return false;
		}
		else if(element.ValueStr() == "avel")
		{
            return false;
		}
		else if(element.ValueStr() == "bodymaterial")
		{
            return false;
		}
		else
		{
			return false;
		}

		assert(false); // we should never get here
		return false;
	}
	/// Visit an element.
	virtual bool VisitExit( const TiXmlElement& element)
	{
		if(element.ValueStr() == "world")
		{
			return true;
		}
		else if(element.ValueStr() == "materials")
		{
			return true;
		}
		else if(element.ValueStr() == "material")
		{
			return true;
		}
		else if(element.ValueStr() == "bodies")
		{
			return true;
		}
		else if(element.ValueStr() == "ground")
		{
			m_world->m_shapes.push_back(m_curBody);
			m_curBody = NULL;
			return true;
		}
		else if(element.ValueStr() == "box")
		{
			m_world->m_shapes.push_back(m_curBody);
			m_curBody = NULL;
			return true;
		}
		else if(element.ValueStr() == "sphere")
		{
			m_world->m_shapes.push_back(m_curBody);
			m_curBody = NULL;
			return true;
		}
		else if(element.ValueStr() == "cylinder")
		{
			m_world->m_shapes.push_back(m_curBody);
			m_curBody = NULL;
			return true;
		}
		else if(element.ValueStr() == "pos")
		{
			return true;
		}
		else if(element.ValueStr() == "vel")
		{
			return true;
		}
		else if(element.ValueStr() == "ori")
		{
			return true;
		}
		else if(element.ValueStr() == "avel")
		{
			return true;
		}
		else if(element.ValueStr() == "bodymaterial")
		{
			return true;
		}
		else
		{
			return false;
		}
	}

private:
	World* m_world;
    World::Shape* m_curBody;
};

World::World(const std::string& filename) 
{
    LoadFromFile(filename);
}

World::~World() 
{
    for (unsigned int i = 0; i < m_shapes.size(); i++)
    {
        delete m_shapes[i];
    }
    m_shapes.clear();
}

void World::Draw()
{
	static GLUquadricObj *quadObj = gluNewQuadric();

    float white[4] = {1.0,1.0,1.0,1.0};
    float grey[4] = {0.8,0.8,0.8,1.0};
    float black[4] = {0.0,0.0,0.0,1.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, grey);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, black);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, black);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, black);

    for (unsigned int i = 0; i < m_shapes.size(); i++)
    {
        vec3 pos = m_shapes[i]->pos;
        if (m_shapes[i]->GetType() == SPHERE)
        {
            Sphere* c = (Sphere*) m_shapes[i];
            glPushMatrix();
			glTranslatef(pos[0]+c->center[0], pos[1]+c->center[1], pos[2]+c->center[2]);
	        glutSolidSphere(c->r, 20, 20);   
            glPopMatrix();
        }
        else if (m_shapes[i]->GetType() == CUBE)
        {
            Cube* c = (Cube*) m_shapes[i];
            glPushMatrix();
            glTranslatef(pos[0], pos[1], pos[2]);
            glScalef(c->hx*2, c->hy*2, c->hz*2);
	        glutSolidCube(1.0);   
            glPopMatrix();
        }
        else if (m_shapes[i]->GetType() == CYLINDER)
        {
            Cylinder* c = (Cylinder*) m_shapes[i];
            vec3 forward = c->end - c->start;
            double height = forward.Length();
            double radius = c->r;

            forward.Normalize();

            vec3 left = vec3(0,1,0)^forward;
            vec3 up;
            if (left.Length() < 0.0001)
            {
                up = forward^vec3(1,0,0);
                left = up^forward;
            }
            else
            {
                up = forward^left;
            }

            float m[16];
            m[0] = left[0]; m[4] = up[0]; m[8] = forward[0];  m[12] = 0; 
            m[1] = left[1]; m[5] = up[1]; m[9] = forward[1];  m[13] = 0; 
            m[2] = left[2]; m[6] = up[2]; m[10] = forward[2]; m[14] = 0; 
            m[3] = 0.0;  m[7] = 0.0;  m[11] = 0.0;  m[15] = 1.0;

	        glPushMatrix();
			glTranslated(c->start[0], c->start[1], c->start[2]);
            glMultMatrixf(m); 
	        gluQuadricDrawStyle(quadObj, GLU_FILL);
	        gluQuadricNormals(quadObj, GLU_SMOOTH);
	        gluCylinder(quadObj, radius, radius, height, 12, 12);            
			//endCaps
			glPushMatrix();
			gluDisk(quadObj, 0, radius, 12, 12);
			glTranslated(0, 0, height);
			gluDisk(quadObj, 0, radius, 12, 12);
			glPopMatrix();
			glPopMatrix();
        }
        else if (m_shapes[i]->GetType() == GROUND)
        {
            glBegin(GL_QUADS);
                glNormal3d(0,1,0);
                glVertex3f(100, 0.0, -100.0);
                glVertex3f(100, 0.0, 100.0);
                glVertex3f(-100, 0.0, 100);
                glVertex3f(-100, 0.0, -100);
            glEnd();
        }
    }
}

void World::LoadFromFile(const std::string& filename)
{
	TiXmlDocument doc(filename);
	bool success = doc.LoadFile();
	if(success)
	{
	    XMLWorldVisitor v(this);
	    doc.Accept(&v);
	}
}

/*
          glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
              glDisable(GL_LIGHTING);
              glBegin(GL_LINES);
                 glColor3f(0.25, 0.25, 0.25);
                 for (int i = -20; i < 20; i++)
                 {
                    glVertex3f(i, 0.0, -20.0);
                    glVertex3f(i, 0.0, 20.0);

                    glVertex3f(-20, 0.0, i);
                    glVertex3f(20, 0.0, i);
                 }
              glEnd();
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
          */