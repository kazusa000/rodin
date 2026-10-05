/* Distributed under the Boost Software License, Version 1.0. */
#ifndef KELVIN_PERIODIC_LEVEL_SET_CUT_H
#define KELVIN_PERIODIC_LEVEL_SET_CUT_H
#include <map>
#include <set>
#include <Eigen/LU>
#include "PeriodicCuts.h"

namespace KelvinBall
{
  /** Cut the actual P1 zero set before MMG optimization. Polygon fans use
   * face barycentres, not orientation-dependent quad diagonals. The paired
   * input triangles and identified nodal values therefore induce the same
   * rotated subdivision. No interface point is projected or repaired.
   */
  template<class Field>
  MMG::Mesh periodicLevelSetCut(const Field& phi)
  {
    const auto& source = phi.getFiniteElementSpace().getMesh();
    using Polygon = std::vector<Index>;
    using Triangle = std::array<Index, 3>;
    using Tet = std::array<Index, 4>;
    std::vector<Math::SpatialPoint> points;
    for (Index i=0; i<source.getVertexCount(); ++i)
      points.push_back(source.getVertexCoordinates(i));
    std::map<std::array<Index,2>, Index> intersections;
    std::map<Polygon, Index> centres;
    std::map<Triangle, Attribute> boundary;
    std::map<Triangle, std::set<Attribute>> materials;
    std::vector<std::pair<Tet, Attribute>> cells;
    const auto sorted = [](Triangle key) { std::sort(key.begin(), key.end()); return key; };
    std::map<Triangle, Attribute> inputLabels;
    for (auto f=source.getPolytope(2); f; ++f)
      if (f->getAttribute() && f->getAttribute()!=Gamma)
        inputLabels[sorted({f->getVertices()[0],f->getVertices()[1],f->getVertices()[2]})]=*f->getAttribute();
    const auto value = [&](Index i) { return phi.getData()(phi.getFiniteElementSpace().getDOFs(0,i)(0)); };
    const auto intersect = [&](Index i, Index j) -> Index {
      if (value(i)==0) return i;
      if (value(j)==0) return j;
      if (i>j) std::swap(i,j);
      const std::array<Index,2> key{i,j};
      const auto found=intersections.find(key);
      if (found!=intersections.end()) return found->second;
      const Real t=value(i)/(value(i)-value(j));
      if (!(t>0 && t<1)) throw std::runtime_error("Invalid P1 cut intersection.");
      const Index index=points.size();
      points.push_back((1-t)*points[i]+t*points[j]);
      intersections.emplace(key,index);
      return index;
    };
    const auto triangulate = [&](const Polygon& polygon) {
      std::vector<Triangle> triangles;
      if (polygon.size()<3) return triangles;
      if (polygon.size()==3) { triangles.push_back({polygon[0],polygon[1],polygon[2]}); return triangles; }
      Polygon key=polygon; std::sort(key.begin(),key.end());
      auto found=centres.find(key);
      Index centre;
      if (found==centres.end()) {
        Math::SpatialPoint x=Math::SpatialPoint::Zero(3);
        for (Index i:key) x+=points[i];
        x/=Real(key.size()); centre=points.size(); points.push_back(x); centres.emplace(key,centre);
      } else centre=found->second;
      for (size_t i=0;i<polygon.size();++i)
        triangles.push_back({polygon[i],polygon[(i+1)%polygon.size()],centre});
      return triangles;
    };
    const std::array<std::array<size_t,3>,4> faces{{{1,2,3},{0,3,2},{0,1,3},{0,2,1}}};
    const auto addCell = [&](Tet tet, Attribute material) {
      Eigen::Matrix3d a;
      for (size_t i=0;i<3;++i) for (size_t j=0;j<3;++j)
        a(j,i)=points[tet[i+1]](j)-points[tet[0]](j);
      const Real determinant=a.determinant();
      if (determinant==0 || !std::isfinite(determinant)) throw std::runtime_error("Degenerate P1 cut cell.");
      if (determinant<0) std::swap(tet[0],tet[1]);
      cells.emplace_back(tet,material);
      for (const auto& f:faces) materials[sorted({tet[f[0]],tet[f[1]],tet[f[2]]})].insert(material);
    };
    for (auto cell=source.getCell();cell;++cell) {
      const auto& vertices=cell->getVertices();
      Tet tet{vertices[0],vertices[1],vertices[2],vertices[3]};
      size_t negative=0,positive=0;
      for (Index i:tet) { negative+=value(i)<0; positive+=value(i)>0; }
      if (!negative && !positive) throw std::runtime_error("P1 level set vanishes on a tetrahedron.");
      if (!negative || !positive) {
        addCell(tet,negative?Obstacle:Fluid);
        for (const auto& f:faces) {
          const Triangle key=sorted({tet[f[0]],tet[f[1]],tet[f[2]]});
          if (inputLabels.count(key)) boundary[key]=inputLabels.at(key);
        }
        continue;
      }
      std::set<Index> interface;
      for (Index i:tet) if (value(i)==0) interface.insert(i);
      for (size_t i=0;i<4;++i) for (size_t j=i+1;j<4;++j)
        if ((value(tet[i])<0 && value(tet[j])>0)||(value(tet[i])>0 && value(tet[j])<0))
          interface.insert(intersect(tet[i],tet[j]));
      Polygon contour(interface.begin(),interface.end());
      if (contour.size()<3) throw std::runtime_error("Incomplete P1 cut contour.");
      Math::SpatialPoint centre=Math::SpatialPoint::Zero(3);
      for (Index i:contour) centre+=points[i]; centre/=Real(contour.size());
      const Math::SpatialPoint u=(points[contour[0]]-centre).normalized();
      const Math::SpatialPoint n=(points[contour[1]]-points[contour[0]]).cross(points[contour[2]]-points[contour[0]]).normalized();
      const Math::SpatialPoint v=n.cross(u);
      std::sort(contour.begin(),contour.end(),[&](Index a,Index b) {
        const Math::SpatialPoint x=points[a]-centre,y=points[b]-centre;
        return std::atan2(x.dot(v),x.dot(u))<std::atan2(y.dot(v),y.dot(u));
      });
      for (const int side:{-1,1}) {
        std::vector<Triangle> shell=triangulate(contour);
        for (const auto& f:faces) {
          Polygon polygon;
          for (size_t k=0;k<3;++k) {
            const Index a=tet[f[k]], b=tet[f[(k+1)%3]];
            const bool insideA=side*value(a)>=0,insideB=side*value(b)>=0;
            if (insideA) polygon.push_back(a);
            if (insideA!=insideB) polygon.push_back(intersect(a,b));
          }
          Polygon unique;
          for (Index i:polygon) if (unique.empty()||unique.back()!=i) unique.push_back(i);
          if (unique.size()>1 && unique.front()==unique.back()) unique.pop_back();
          const auto triangles=triangulate(unique);
          shell.insert(shell.end(),triangles.begin(),triangles.end());
          const Triangle old=sorted({tet[f[0]],tet[f[1]],tet[f[2]]});
          if (inputLabels.count(old)) for (const auto& triangle:triangles) boundary[sorted(triangle)]=inputLabels.at(old);
        }
        std::set<Index> body;
        for (const auto& triangle:shell) for (Index i:triangle) body.insert(i);
        Math::SpatialPoint x=Math::SpatialPoint::Zero(3);
        for (Index i:body) x+=points[i]; x/=Real(body.size());
        const Index interior=points.size(); points.push_back(x);
        for (const auto& triangle:shell) addCell({triangle[0],triangle[1],triangle[2],interior},side<0?Obstacle:Fluid);
      }
    }
    Mesh::Builder builder; builder.initialize(3).nodes(points.size());
    for (const auto& x:points) builder.vertex(x);
    for (const auto& [tet,material]:cells) {
      IndexArray vertices(4);
      for (size_t i=0;i<4;++i) vertices(i)=tet[i];
      Index index; builder.polytope(Polytope::Type::Tetrahedron,std::move(vertices),index);
      builder.attribute({3,index},material);
    }
    Mesh mesh=builder.finalize();
    for (auto f=mesh.getPolytope(2);f;++f) {
      const auto& v=f->getVertices(); const Triangle key=sorted({v[0],v[1],v[2]});
      if (boundary.count(key)) mesh.setAttribute({2,f->getIndex()},boundary.at(key));
      else if (materials.at(key).size()==2) mesh.setAttribute({2,f->getIndex()},Gamma);
    }
    PeriodicCuts(mesh).checkClosedInterface();
    return MMG::Mesh(mesh);
  }
}
#endif
