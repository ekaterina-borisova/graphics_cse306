#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <cfloat>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <thread>

#include <climits>
#include <thread>
#include <numeric>
#include <iterator>
#include <random>
#include <unistd.h>
#include <initializer_list>
#include <list>

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform (0 ,1);

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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

class Ray {
public:
	Ray(const Vector& O, const Vector& u) : O(O), u(u) {};
	Vector O;
	Vector u;
};
double sqr(double x){
	return x*x;
};


class Geometry{
public:
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const = 0;
	Vector rho;
	bool is_reflective;
	bool is_refractive;
	bool hollow;
};


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};
 

class BoundingBox {
	public:
		Vector Bmin;
		Vector Bmax;
		BoundingBox(){}
		BoundingBox(const Vector& Bmin, const Vector& Bmax): Bmin(Bmin), Bmax(Bmax) {};
		bool intersect(const Ray& r, double& inter_distance) const{
			Vector X = (Bmin - r.O) / r.u;
			Vector Y = (Bmax - r.O) / r.u;

			Vector t0(std::min(X[0], Y[0]), std::min(X[1], Y[1]), std::min(X[2], Y[2]));
			Vector t1(std::max(X[0], Y[0]), std::max(X[1], Y[1]), std::max(X[2], Y[2]));
			
			if(min(t1) > max(t0) && max(t0) > 0){
				inter_distance = max(t0);
				return true;
			}
			return false;
		}
 };

 class Node{
	public:
		int begin;
		int end;
		BoundingBox bbox;
		Node* left;
		Node* right;
};

 
class TriangleMesh : public Geometry{
public:
  	~TriangleMesh() {}
	std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
	BoundingBox bbox;
	Node* root;
	
    TriangleMesh(const Vector& rho, bool is_reflective = false, bool is_refractive = false, bool hollow = false){
		this->rho = rho;
		this->is_reflective = is_reflective;
		this->is_refractive = is_refractive;
		this->hollow = hollow;
	}

	BoundingBox set_bbox(int begin, int end)
	{
		double x_min = DBL_MAX;
		double y_min = DBL_MAX;
		double z_min = DBL_MAX;
		double x_max = -DBL_MAX;
		double y_max = -DBL_MAX;
		double z_max = -DBL_MAX;
		for (int i = begin; i < end; i++)
        {
			Vector x(vertices[indices[i].vtxi][0], vertices[indices[i].vtxj][0], vertices[indices[i].vtxk][0]);
			Vector y(vertices[indices[i].vtxi][1], vertices[indices[i].vtxj][1], vertices[indices[i].vtxk][1]);
			Vector z(vertices[indices[i].vtxi][2], vertices[indices[i].vtxj][2], vertices[indices[i].vtxk][2]);
			x_min = std::min(x_min, min(x));
			y_min = std::min(y_min, min(y));
			z_min = std::min(z_min, min(z));
			x_max = std::max(x_max, max(x));
			y_max = std::max(y_max, max(y));
			z_max = std::max(z_max, max(z));
        }
		Vector Bmin(x_min, y_min, z_min);
		Vector Bmax(x_max, y_max, z_max);
		BoundingBox bbox(Bmin, Bmax);
		return bbox;
	}

	void bvh(Node *node, int begin, int end)
    {

        BoundingBox bbox = set_bbox(begin, end);
        node->bbox = bbox;
        node->begin = begin;
        node->end = end;


        Vector diag = node->bbox.Bmax - node->bbox.Bmin;
        Vector middle_diag = node->bbox.Bmin + diag * 0.5;


        int longest_axis = max_idx(diag);
        int pivot_index = begin;
        for (int i = begin; i < end; i++)
        {
            Vector barycenter = (vertices[indices[i].vtxi] + vertices[indices[i].vtxj] + vertices[indices[i].vtxk]) / 3;

            if (barycenter[longest_axis] < middle_diag[longest_axis])
            {
                std::swap(indices[i], indices[pivot_index]);
                pivot_index++;
            }
        }

        if (pivot_index <= begin|| pivot_index >= (end - 1) || (end - begin) < 1 ) return ;

        node->left = new Node();
        node->right = new Node();
        bvh(node->left, begin, pivot_index);
        bvh(node->right, pivot_index, end);
    }

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const{
		double inter_distance = 0.0;
		if(!root->bbox.intersect(r, inter_distance)) return false;

		t = DBL_MAX;		
		bool intersection = false;
		std::list<Node*> nodes_to_visit;
		nodes_to_visit.push_front(root);
		double best_inter_distance = std::numeric_limits<double >::max();
		while(!nodes_to_visit.empty()){
			Node* curNode = nodes_to_visit.back();
			nodes_to_visit.pop_back();
			// if there is one child , then it is not a leaf , so test the bounding box
			if( curNode->left ) {
				if (curNode->left->bbox.intersect(r, inter_distance)) {
					if(inter_distance < best_inter_distance) { 
						nodes_to_visit.push_back(curNode->left);
					} 
			}
				if (curNode->right->bbox.intersect(r, inter_distance)) { 
					if(inter_distance < best_inter_distance) {
						nodes_to_visit.push_back(curNode->right) ; 
					}
				} 
			} else {
				for(int i = curNode->begin; i < curNode->end; i++){
					TriangleIndices triangle = indices[i];
					Vector A = vertices[triangle.vtxi];
					Vector B = vertices[triangle.vtxj];
					Vector C = vertices[triangle.vtxk];

					Vector e1 = B - A;
					Vector e2 = C - A;

					Vector N_cur = cross(e1, e2);

					double t_cur = dot(A - r.O, N_cur)/ dot(r.u, N_cur);

					if (t_cur > 0 && t_cur < t){
						double beta = dot(e2, cross(A - r.O, r.u))/ dot(r.u, N_cur);
						double gamma = -dot(e1, cross(A - r.O, r.u))/ dot(r.u, N_cur);
						double alpha = 1 - beta - gamma;
						if (beta > 0 && beta < 1 && gamma > 0 && gamma < 1 && alpha > 0 && alpha < 1){
							N = N_cur;
							N.normalize();
							t = t_cur;
							P = A + beta * e1 +  gamma * e2;
							intersection = true;
						}
					}
				}
			}
		}
		return intersection;
	}
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
		std::cout<< 1;
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
		// Scale and translate
		double scale = 0.6;
		Vector translation(0,-10,0);
		for (size_t i = 0; i < vertices.size(); i++)
		{
				vertices[i] = vertices[i] * scale + translation;
		}
		// Allocate memory for root
    	root = new Node();

		// Run the bvh
		bvh(root, 0, indices.size());
        fclose(f);
 
    }
    
};


Vector random_cos(Vector& N){
	double r1 = uniform(engine);
	double r2 = uniform(engine);
	double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
	double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
	double z = sqrt(r2);
	Vector T1;
	if (N[0] == min(N)){
		T1 = Vector(0, N[2], -N[1]);
	}else if(N[1] == min(N)){
		T1 = Vector(N[2], 0, -N[0]);
	}else if(N[2] == min(N)){
		T1 = Vector(N[1],-N[0], 0);
	}
	T1.normalize();
	Vector T2 = cross(N, T1);
	Vector V = x * T1 + y * T2 + z * N;
	return V;
}


class Sphere : public Geometry{
public:
	Vector C;
	double R;
	Sphere(Vector C, double R, Vector rho, bool is_reflective = false, bool is_refractive = false, bool hollow = false){
		this->C = C;
		this->R = R;
		this->rho = rho;
		this->is_reflective = is_reflective;
		this->is_refractive = is_refractive;
		this->hollow = hollow;
	}

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const{
		double delta = sqr(dot(r.u, r.O - C)) - (r.O - C, r.O - C).norm2() + sqr(R);
		if (delta < 0){return false;}
		else{
			double t1 = dot(r.u, C - r.O) - sqrt(delta);
			double t2 = dot(r.u, C - r.O) + sqrt(delta);
			if (t2 < 0){return false;}
			else if(t1 >= 0){t = t1;}
			else{t = t2;}
			P = r.O + t * r.u;
			N = P - C;
			N.normalize();
			if(hollow){
				N = -N;
			}
			return true;
		}
	}
};


class Scene {
public:
	Scene(const Vector& S, const double& I): S(S), I(I) {};
	Vector S;
	double I;
	bool intersect(const Ray& r, Vector& P, Vector& N, double& t, size_t& sphere_idx){
		t = DBL_MAX;
		bool res = false;
		for(size_t i = 0; i < objects.size(); i++){
			double t_tmp;
			Vector P_tmp, N_tmp;
			if(objects[i]->intersect(r, P_tmp, N_tmp, t_tmp)){
				if(t_tmp < t){
					P = P_tmp;
					N = N_tmp;
					t = t_tmp;
					sphere_idx = i;
					res = true;
				}

			}
		}
		return res;

	}

	Vector get_color(const Ray& r, int ray_depth){
		if (ray_depth < 0){
			return Vector(0., 0., 0.);
		}

		double t;
		Vector P, N;
		size_t sphere_idx;
		double eps = 0.0001;
		Vector L(0,0,0);
		if(intersect(r, P, N, t, sphere_idx)){

			if (objects[sphere_idx]->is_reflective){
				Vector w_r = r.u - 2 * dot(r.u, N) * N;
				Ray reflection_ray(P + eps * N, w_r);
				return get_color(reflection_ray, ray_depth - 1);
			}
			if (objects[sphere_idx]->is_refractive){
				double n1 = 1.0;
				double n2 = 1.5;
				Vector Ntmp = N;
				if(dot(r.u, N) > 0){
					Ntmp = -Ntmp;
					std::swap(n1,n2);
				}
				Vector w_tT = (n1 / n2) * (r.u - dot(r.u, Ntmp) * Ntmp);
				double tmp = 1 - pow((n1 / n2), 2) * (1 - pow(dot(r.u, N), 2));
				double k0 = pow((n1 - n2),2) / pow((n1 + n2),2);
				double R = k0 + (1 - k0) * pow((1 - std::abs(dot(N, r.u))), 5);
				double u = uniform(engine);
				if (tmp < 0 or u < R){
					Vector w_r = r.u - 2 * dot(r.u, N) * N;
					Ray reflection_ray(P + eps * N, w_r);
					return get_color(reflection_ray, ray_depth - 1);
				}
				Vector w_tN = -Ntmp * sqrt(tmp);
				Vector w_t = w_tT + w_tN;
				Ray refraction_ray(P + eps * (-Ntmp), w_t);
				return get_color(refraction_ray, ray_depth - 1);

			}
			double d = (S - P).norm();
			Vector w_i = S - P;
			w_i.normalize();
			Ray light(P + eps * N, w_i);
			bool visible = true;
			size_t inter_sphere_idx;
			double t_inter;
			Vector P_inter, N_inter;
			if (intersect(light, P_inter, N_inter, t_inter, inter_sphere_idx)){
				if(t_inter <= d){
					visible  = false;
				}
			}
			L = I / (4 * M_PI * sqr(d)) * int(visible) * objects[sphere_idx]->rho / M_PI * (std::max(0., dot(N, w_i)));
			Vector V = random_cos(N);
			Ray random_ray(P + eps * N, V);
			L = L + objects[sphere_idx]->rho * get_color(random_ray, ray_depth - 1);
		}
		return L;

	}
	void add_object(Geometry* geo) {
		objects.push_back(geo);
	}
	std::vector<Geometry* > objects;
};


void boxMuller(double stdev , double &x, double &y) { 
	double r1 = uniform ( engine ) ;
	double r2 = uniform ( engine ) ;
	x = sqrt(-2 * log(r1))*cos(2 * M_PI*r2)*stdev; 
	y = sqrt(-2 * log(r1))*sin(2 * M_PI*r2)*stdev;
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " [scenario number, where 1 - spheres, 2 - cat]\n";
        return 1;
    }
	int W = 512;
	int H = 512;
	Sphere* ceiling = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
	Sphere* floor = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
	Sphere* front = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
	Sphere* back = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 0.5));
	Sphere* left = new Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0));
	Sphere* right = new Sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1));
	Vector camera_center(0, 0, 55);
	const double alpha = 60.*M_PI/180.;
	const double I = 2E10;
	const int K = 64;
	Vector S(-10, 20, 40);
	Scene scene(S, I);
	Sphere* sphere1;
	Sphere* sphere2;
	Sphere* sphere3;
	Sphere* sphere3b;
	TriangleMesh* cat;
	switch (std::atoi(argv[1])) {
        case 1:
			// Spheres image
            sphere1 = new Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1), false, true);
			sphere2 = new Sphere(Vector(-20, 0, 0), 10, Vector(1, 1, 1), true);
			sphere3 = new Sphere(Vector(20, 0, 0), 10, Vector(1, 1, 1), false, true);
			sphere3b = new Sphere(Vector(20, 0, 0), 9.5, Vector(1, 1, 1), false, true, true);
			scene.add_object(sphere1);
			scene.add_object(sphere2);
			scene.add_object(sphere3);
			scene.add_object(sphere3b);
            break;
        case 2:
			// Cat image
            cat = new TriangleMesh(Vector(1, 1, 1));
			cat->readOBJ("models/cat.obj");
			scene.add_object(cat);
            break;
        default:
            std::cerr << "Invalid scenario number\n";
            return 1;
    }

	
	
	//scene.add_object(sphere1);
	// scene.add_object(sphere2);
	// scene.add_object(sphere3);
	// scene.add_object(sphere3b);

	//scene.add_object(sphere1);
	// scene.add_object(sphere2);
	// scene.add_object(sphere3);
	// scene.add_object(sphere3b);
	scene.add_object(ceiling);
	scene.add_object(floor);
	scene.add_object(front);
	scene.add_object(back);
	scene.add_object(left);
	scene.add_object(right);

	std::vector<unsigned char> image(W * H * 3, 0);
	
	#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			double x,y;
			Vector color(0,0,0);
			for (size_t k = 0; k < K; ++k){
				boxMuller(1,x,y);
				Vector direction(x + j + 0.5 - W / 2, H - (y + i) - 1 + 0.5 - H / 2, -W / (2 * tan(alpha / 2)));
				direction.normalize();
				Ray ray = Ray(camera_center, direction);
				Vector tmp_color = scene.get_color(ray, 5);
				color = color + tmp_color / K;
			}

			double gamma = 2.2;
			color = pow(color, 1/gamma);
			image[(i * W + j) * 3 + 0] = color[0];
			image[(i * W + j) * 3 + 1] = color[1];
			image[(i * W + j) * 3 + 2] = color[2];
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}
