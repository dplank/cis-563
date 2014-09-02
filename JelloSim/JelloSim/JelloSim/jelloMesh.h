#ifndef jelloMesh_H_
#define jelloMesh_H_

#include <vector>
#include <map>
#include "vec.h"
#include "World.h"

class JelloMesh
{
public:
    JelloMesh();
    virtual ~JelloMesh();

    // Update our JelloMesh (handles collisions and forces)
    virtual void Update(double dt, 
        const World& world,
        const vec3& externalForces = vec3(0.0, -9.8, 0));

    // Draw our JelloMesh
    virtual void Draw(const vec3& eyePos);

    // Reset to the initial state
    virtual void Reset();

    // Set/Get our JelloMesh resolution
    virtual void SetGridSize(int cols, int rows, int stacks);
    virtual int GetGridCols() const;
    virtual int GetGridRows() const;
    virtual int GetGridStacks() const;

    // Set/Get our JelloMesh size (in world coordinates)
    virtual void SetSize(float width, float height, float depth);
    virtual float GetWidth() const;
    virtual float GetHeight() const;
    virtual float GetDepth() const;

    // Set/Get our numerical integration type
    enum IntegrationType { EULER, MIDPOINT, RK4 };
    virtual void SetIntegrationType(IntegrationType type);
    virtual IntegrationType GetIntegrationType() const;

    // Set/Get our draw flags
    virtual void SetDrawFlags(unsigned int flags);
    virtual unsigned int GetDrawFlags() const;

    // Spring types
    enum SpringType { STRUCTURAL = 0x1, SHEAR = 0x2, BEND = 0x4 }; 

    int GetIndex(int i, int j, int k) const;
    void GetCell(int idx, int& i, int &j, int &k) const;


protected:

    class Particle;
    class Spring;
    friend class FaceMesh;
    friend class TestJelloMesh;

    typedef std::vector<std::vector<std::vector<Particle>>> ParticleGrid;
    Particle& GetParticle(ParticleGrid& grid, int i, int j, int k);
    Particle& GetParticle(ParticleGrid& grid, int idx);
    const Particle& GetParticle(const ParticleGrid& grid, int i, int j, int k) const;
    const Particle& GetParticle(const ParticleGrid& grid, int idx) const;

    bool isInterior(const Spring& s) const;
    bool isInterior(int idx) const;
    bool isInterior(int i, int j, int k) const;

    virtual void InitJelloMesh();
    virtual void AddStructuralSpring(Particle& p1, Particle& p2);
    virtual void AddBendSpring(Particle& p1, Particle& p2);
    virtual void AddShearSpring(Particle& p1, Particle& p2);

    class Intersection;
	virtual void CheckForCollisions(ParticleGrid& grid, const World& world);
	virtual void ResolveCollisions(ParticleGrid& grid);
	virtual void ResolveContacts(ParticleGrid& grid);
    virtual bool FloorIntersection(Particle& p, Intersection& intersection);
    virtual bool CylinderIntersection(Particle& p, World::Cylinder* cylinder, Intersection& intersection);
	virtual bool SphereIntersection(Particle& p, World::Sphere* sphere, Intersection& intersection);
	virtual bool CubeIntersection(Particle& p, World::Cube* cube, Intersection& intersection);

    virtual void ComputeForces(ParticleGrid& grid);
	virtual void EulerIntegrate(double dt);
	virtual void MidPointIntegrate(double dt);
	virtual void RK4Integrate(double dt);

    enum Face {XLEFT, XRIGHT, YTOP, YBOTTOM, ZFRONT, ZBACK};
    class FaceMesh
    {
    public:
        std::map<int,std::vector<int>> m_neighbors; 
        std::vector<std::vector<int>> m_strips;
        double distToEye;

        FaceMesh(const JelloMesh& m, Face f);
        void Draw(const JelloMesh& m);
        void DrawNormals(const JelloMesh& m);

        void CalcDistToEye(const JelloMesh& m, const vec3& eyePos);
        static bool compare(const FaceMesh& one, const FaceMesh& other);
    };

    void DrawMesh(const vec3& eyePos);
    void DrawSprings(double a);
    void DrawCollisionNormals();
    void DrawForces();

protected:

    int m_cols, m_rows, m_stacks;
    float m_width, m_height, m_depth;
    unsigned int m_drawflags;
    vec3 m_externalForces;

    IntegrationType m_integrationType;
    std::vector<FaceMesh> m_mesh;
    ParticleGrid m_vparticles;

    std::vector<Spring> m_vsprings;
    std::vector<Intersection> m_vcontacts;
    std::vector<Intersection> m_vcollisions;

public:

    static double g_structuralKs;
    static double g_structuralKd;
    static double g_shearKs;
    static double g_shearKd;
    static double g_bendKs;
    static double g_bendKd;
    static double g_penaltyKs;
    static double g_penaltyKd;

    static const unsigned int MESH = 0x10;
    static const unsigned int NORMALS = 0x100;
    static const unsigned int FORCES = 0x1000;

protected:

    class Particle
    {
    public:
        Particle();
        Particle(const Particle& p);
        Particle& operator=(const Particle& p);
        Particle(int idx, const vec3& pos, const vec3& vel = vec3(0,0,0), double m = 1);

        int index;
        vec3 position;
        vec3 velocity;
        vec3 force;
        double mass;

        static Particle EMPTY;
    };

    class Spring
    {
    public:
        Spring();
        Spring(const Spring& p);
        Spring& operator=(const Spring& p);
        Spring(SpringType t, int p1, int p2, 
            double Ks, double Kd, double restLen);

        SpringType m_type;
        int m_p1;
        int m_p2;
        double m_Ks;
        double m_Kd;
        double m_restLen;
    };

    enum IntersectionType { CONTACT, COLLISION };
    class Intersection
    {
    public:
        Intersection();
        Intersection(const Intersection& p);
        Intersection& operator=(const Intersection& p);
        Intersection(IntersectionType type, int p, const vec3& normal, double d = 0);

        int m_p;
        vec3 m_normal;
        double m_distance;
        IntersectionType m_type;
    };
};

#endif