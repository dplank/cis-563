#ifndef WORLD_H
#define WORLD_H

#include <string>
#include <vector>
#include "vec.h"

class World
{
public:
    World(const std::string& filename);
    virtual ~World();

    void LoadFromFile(const std::string& filename);
    virtual void Draw();

public:
    enum ShapeType {SPHERE, GROUND, CUBE, CYLINDER};
    class Shape
    {
    public:
        vec3 pos;    
        ShapeType GetType() const { return type; }
    protected:
        Shape(ShapeType t) : type(t), pos(0,0,0) {}
        ShapeType type;
    };

    class Ground : public Shape
    {
    public:
        Ground() : Shape(GROUND) {}
    };

    class Sphere : public Shape
    {
    public:
        Sphere() : Shape(SPHERE), r(1.0) {}
        double r;
		vec3 center;
    };

    class Cube : public Shape
    {
    public:
        Cube() : Shape(CUBE), hx(0.5), hy(0.5), hz(0.5) {}
        double hx, hy, hz;
    };

    class Cylinder : public Shape
    {
    public:
        Cylinder() : Shape(CYLINDER), r(1.0) {}
        vec3 start, end;
        double r;
    };

    std::vector<Shape*> m_shapes;
};

#endif