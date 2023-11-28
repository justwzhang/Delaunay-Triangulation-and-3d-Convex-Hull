#include <time.h>
#include "PowerDiagram.h"

void PowerDiagram::CPowerDiagram::init(int num_pts)
{
    std::vector<CPoint*> pts;

    // 1. generate random points
    // https://mathworld.wolfram.com/DiskPointPicking.html
    srand((unsigned) time(NULL));
    double n[2];
    for (int i = 0; i < num_pts; ++i)
    {
        n[0] = rand() / double(RAND_MAX);               // [0, 1]
        n[1] = 2.0 * M_PI * rand() / double(RAND_MAX);  // [0, 2*pi]

        double x = std::sqrt(n[0]);
        double y = x;
        x *= std::cos(n[1]);
        y *= std::sin(n[1]);

        CPoint *p = new CPoint(x, y, 0);
        pts.push_back(p);
    }

    // 2. lift the points onto paraboloid z = (x^2 + y^2)/2

    //insert your code here
    for (CPoint* point: pts) {
        double x = (*point) * CPoint(1, 0, 0);
        double y = (*point) * CPoint(0, 1, 0);
        double z = (pow(x, 2)+ pow(y, 2)) / 2;
        CPoint *zPoint = new CPoint(0, 0, z);
        *point += *zPoint;
    }

    // 3. feed into CPowerDiagram
    init(pts);
}

void PowerDiagram::CPowerDiagram::init(const std::vector<CPoint*>& points) 
{
    m_pts.clear();
    
    for (const auto p : points)
    {
        m_pts.push_back(p);
    }
}

void PowerDiagram::CPowerDiagram::calc_delaunay()
{
    // 1. calculate the convex hull
    CConvexHull ch;
    ch.init(m_pts);
    ch.construct();

    // 2. remove the faces with upward normal vector
    using M = CConvexHullMesh;
    std::vector<M::CFace*> removed_faces;
    CPoint up(0, 0, 1);
    for (M::FaceIterator fiter(&ch.hull()); !fiter.end(); ++fiter)
    {
        M::CFace* pF = *fiter;
        //insert your code here
        if (pF->normal() * up > 0) removed_faces.push_back(pF);
    }
    
    ch.hull().remove_faces(removed_faces);

    // 3. copy the result
    m_mesh.copy(&ch.hull());
}

void PowerDiagram::CPowerDiagram::calc_voronoi() 
{ 
    for (CMesh::FaceIterator fiter(&m_mesh); !fiter.end(); ++fiter)
    {
        CMesh::CFace* pF = *fiter;
        CMesh::CDart* pD = NULL;
        pD = m_mesh.face_dart(pF);
        CPoint& a = m_mesh.dart_target(pD)->point();
        pD = m_mesh.dart_next(pD);
        CPoint& b = m_mesh.dart_target(pD)->point();
        pD = m_mesh.dart_next(pD);
        CPoint& c = m_mesh.dart_target(pD)->point();

        double x = 0, y = 0, z = 0;
        
        //insert your code here

        //version 3
        // These two links helped a lot for this part
        // https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/
        // https://www.geeksforgeeks.org/program-for-point-of-intersection-of-two-lines/
        //create the line info
        //ax + by = c where CPoint(a, b, c)
        CPoint abLineinfo = CPoint(b[1] - a[1], a[0] - b[0], (b[1] - a[1]) * a[0] + (a[0] - b[0]) * a[1]);

        CPoint bcLineinfo = CPoint(c[1] - b[1], b[0] - c[0], (c[1] - b[1]) * b[0] + (b[0] - c[0]) * b[1]);


        //calculate perp bisectors
        double newA, newB, newC, newD, newE, newF;
        //index of a = 0, b = 1, c = 2
        double mABx = (a[0] + b[0]) / 2;
        double mABy = (a[1] + b[1]) / 2;
        newC = -abLineinfo[1] * mABx + abLineinfo[0] * mABy;
        newA = -abLineinfo[1];
        newB = abLineinfo[0];
        abLineinfo = CPoint(newA, newB, newC);

        double mBCx = (c[0] + b[0]) / 2;
        double mBCy = (c[1] + b[1]) / 2;
        newF = -bcLineinfo[1] * mBCx + bcLineinfo[0] * mBCy;
        newD = -bcLineinfo[1];
        newE = bcLineinfo[0];
        bcLineinfo = CPoint(newD, newE, newF);

        //determinate to find point of intersection
        double a1 = abLineinfo[0];
        double b1 = abLineinfo[1];
        double c1 = abLineinfo[2];
        double a2 = bcLineinfo[0];
        double b2 = bcLineinfo[1];
        double c2 = bcLineinfo[2];

        double determinate = a1 * b2 - a2 * b1;
        x = (b2 * c1 - b1 * c2) / determinate;
        y = (a1 * c2 - a2 * c1) / determinate;

        z = (pow(x, 2) + pow(y, 2)) / 2;
        pF->dual_point() = CPoint(x,y,z);
        
    }

}
