#include <iostream>
#include <map>

#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <stdio.h>
//#include "libigl-main/include/igl/writeOBJ.h"

int sgn(double x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, const Vector& b) {
	return Vector(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
bool operator<(const Vector& a, const Vector& b) {
	return (a[0] < b[0]) && (a[1] < b[1]) && (a[2] < b[2]);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector pow(const Vector& a, double b) {
	return Vector(std::min(std::pow(a[0], b),255.), std::min(std::pow(a[1], b), 255.), std::min(std::pow(a[2], b),255.));
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
double min(const Vector& a){
	return std::min(std::min(a[0], a[1]), a[2]);
}

double max(const Vector& a){
	return std::max(std::max(a[0], a[1]), a[2]);
}

double max_idx(const Vector& a){
	int idx = 0;
	for (size_t i = 0; i < 3; i++){
		if (a[i] > a[idx]){
			idx = i;
		}
	}
	return idx;
}

class Edge {
public:
    int a, b;
    Edge(int u = -1, int v = -1): a(u), b(v) {}
};

bool operator<(const Edge& e1, const Edge& e2) {
    if (std::min(e1.a, e1.b) < std::min(e2.a, e2.b))
        return true; 
    
    return std::max(e1.a, e1.b) < std::max(e1.a, e2.b);
}

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};
 
 
class TriangleMesh {
public:
  ~TriangleMesh() {}
    TriangleMesh() {};
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        // navigate to folder 92-mask 
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
 
    }

    void writeOBJ(const char* filename) {
        FILE* f;
        f = fopen(filename, "w");
        fprintf(f, "mtllib default.mtl\n");
        fprintf(f, "usemtl material_0\n");
        std::cout << vertexcolors.size() << std::endl;
        for (int i = 0; i < (int)vertices.size(); i++) {
            fprintf(f, "v %f %f %f %f %f %f\n", vertices[i][0], vertices[i][1], vertices[i][2], vertexcolors[i][0], vertexcolors[i][1], vertexcolors[i][2]);
        }
        std::cout << "no seg 2" << std::endl;
        for (int i = 0; i < (int)normals.size(); i++) {
            fprintf(f, "vn %f %f %f\n", normals[i][0], normals[i][1], normals[i][2]);
        }
        std::cout << "no seg 3" << std::endl;
        for (int i = 0; i < (int)uvs.size(); i++) {
            fprintf(f, "vt %f %f\n", uvs[i][0], uvs[i][1]);
        }
        std::cout << "no seg 4" << std::endl;
        for (int i = 0; i < (int)indices.size(); i++) {
            fprintf(f, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", indices[i].vtxi + 1, indices[i].uvi + 1, indices[i].ni + 1, indices[i].vtxj + 1, indices[i].uvj + 1, indices[i].nj + 1, indices[i].vtxk + 1, indices[i].uvk + 1, indices[i].nk + 1);
        }
        std::cout << "no seg 5" << std::endl;
        fclose(f);
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    
};



void modifyMesh(TriangleMesh& mesh) {
    std::map<Edge, std::vector<size_t> > map_edge_tri;
    std::map<int, std::vector<int> > map_vertex_tri;

    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        int a = mesh.indices[i].vtxi;
        int b = mesh.indices[i].vtxj;
        int c = mesh.indices[i].vtxk;

        map_edge_tri[Edge(a, b)].push_back(i);
        map_edge_tri[Edge(b, c)].push_back(i);
        map_edge_tri[Edge(c, a)].push_back(i);

        map_vertex_tri[a].push_back(i);
        map_vertex_tri[b].push_back(i);
        map_vertex_tri[c].push_back(i);
    }

    std::vector<Edge> boundary_edges;
    std::vector<bool> is_boundary(mesh.vertices.size(), false);

    for (const auto& entry : map_edge_tri) {
        if (entry.second.size() == 1) {
            boundary_edges.push_back(entry.first);
            is_boundary[entry.first.a] = true;
            is_boundary[entry.first.b] = true;
        }
    }

    std::vector<Edge> ordered_boundary_edges(boundary_edges.size());
    ordered_boundary_edges[0] = boundary_edges[0];

    for (size_t i = 1; i < boundary_edges.size(); ++i) {
        for (size_t j = 0; j < boundary_edges.size(); ++j) {
            if (boundary_edges[j].a == ordered_boundary_edges[i - 1].b) {
                ordered_boundary_edges[i] = boundary_edges[j];
                break;
            }
        }
    }

    for (size_t i = 0; i < ordered_boundary_edges.size(); ++i) {
        double theta = i / static_cast<double>(ordered_boundary_edges.size()) * 2 * M_PI;
        Vector circle_vtx(0.5 * std::cos(theta), 0.5 + 0.5 * std::sin(theta), 0.);
        mesh.vertices[ordered_boundary_edges[i].a] = circle_vtx;
    }

    for (int iter = 0; iter < 10; ++iter) {
        std::vector<Vector> updated_vertices = mesh.vertices;

        for (size_t i = 0; i < mesh.vertices.size(); ++i) {
            if (is_boundary[i])
                continue;

            Vector sum_neighbors(0., 0., 0);
            size_t total_neighbors = 0;

            for (size_t j = 0; j < map_vertex_tri[i].size(); ++j) {
                size_t tri_idx = map_vertex_tri[i][j];
                const auto& indices = mesh.indices[tri_idx];

                if (indices.vtxi != i)
                    sum_neighbors = sum_neighbors + mesh.vertices[indices.vtxi];
                if (indices.vtxj != i)
                    sum_neighbors = sum_neighbors + mesh.vertices[indices.vtxj];
                if (indices.vtxk != i)
                    sum_neighbors = sum_neighbors + mesh.vertices[indices.vtxk];

                total_neighbors += 2;
            }

            updated_vertices[i] = sum_neighbors / static_cast<double>(total_neighbors);
        }

        mesh.vertices = updated_vertices;
    }
}

int main() {
    TriangleMesh mesh;
    std::cout << "Reading OBJ file..." << std::endl;
    mesh.readOBJ("Mask.obj");
    std::cout<< mesh.vertexcolors.size() << std::endl;
   
    std::cout << "Done reading OBJ file." << std::endl;
    modifyMesh(mesh);
    std::cout << "Done modifying." << std::endl;
  
    mesh.writeOBJ("Mask_result.obj");

    

    return 0;
}
