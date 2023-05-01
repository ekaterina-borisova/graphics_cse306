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
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
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
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    
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
class Ray {
public:
	Ray(const Vector& O, const Vector& u) : O(O), u(u) {};
	Vector O;
	Vector u;
};
double sqr(double x){
	return x*x;
};

class Sphere {
public:
	Sphere(const Vector& C, double R, const Vector& rho, bool is_reflective = false, bool is_refractive = false, bool hollow = false):
	C(C), R(R), rho(rho), is_reflective(is_reflective), is_refractive(is_refractive), hollow(hollow){};
	Vector C;
	double R;
	Vector rho;
	bool is_reflective;
	bool is_refractive;
	bool hollow;

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const {
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
		for(int i = 0; i < objects.size(); i++){
			double t_tmp;
			Vector P_tmp, N_tmp;
			if(objects[i].intersect(r, P_tmp, N_tmp, t_tmp)){
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

			if (objects[sphere_idx].is_reflective){
				Vector w_r = r.u - 2 * dot(r.u, N) * N;
				Ray reflection_ray(P + eps * N, w_r);
				return get_color(reflection_ray, ray_depth - 1);
			}
			if (objects[sphere_idx].is_refractive){
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
			L = I / (4 * M_PI * sqr(d)) * int(visible) * objects[sphere_idx].rho / M_PI * (std::max(0., dot(N, w_i)));

			// Vector V = random_cos(N);
			// Ray random_ray(P + eps * N, V);
			// L = L + objects[sphere_idx].rho * get_color(random_ray, ray_depth - 1);

		}
		return L;

	}
	void add_sphere(const Sphere& s) {
		objects.push_back(s);
	}
	std::vector<Sphere> objects;
};


void boxMuller(double stdev , double &x, double &y) { 
	double r1 = uniform ( engine ) ;
	double r2 = uniform ( engine ) ;
	x = sqrt(-2 * log(r1))*cos(2 * M_PI*r2)*stdev; 
	y = sqrt(-2 * log(r1))*sin(2 * M_PI*r2)*stdev;
}

int main() {
	int W = 512;
	int H = 512;
	bool frensel = true;

	Sphere sphere1(Vector(0, 0, 0), 10, Vector(1, 1, 1), false, true);
	Sphere sphere2(Vector(-20, 0, 0), 10, Vector(1, 1, 1), true);
	Sphere sphere3(Vector(20, 0, 0), 10, Vector(1, 1, 1), false, true);
	Sphere sphere3b(Vector(20, 0, 0), 9.5, Vector(1, 1, 1), false, true, true);
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
	Sphere floor(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
	Sphere front(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
	Sphere back(Vector(0, 0, 1000), 940, Vector(1, 0, 0.5));
	Sphere left(Vector(-1000, 0, 0), 940, Vector(1, 1, 0));
	Sphere right(Vector(1000, 0, 0), 940, Vector(0, 1, 1));
	Vector camera_center(0, 0, 55);
	const double alpha = 60.*M_PI/180.;
	const double I = 2E10;
	const int K = 1000;
	Vector S(-10, 20, 40);
	Scene scene(S, I);
	scene.add_sphere(sphere1);
	scene.add_sphere(sphere2);
	scene.add_sphere(sphere3);
	scene.add_sphere(sphere3b);
	scene.add_sphere(ceiling);
	scene.add_sphere(floor);
	scene.add_sphere(front);
	scene.add_sphere(back);
	scene.add_sphere(left);
	scene.add_sphere(right);

	std::vector<unsigned char> image(W * H * 3, 0);

	#pragma omp parallel for or schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			double x,y;
			Vector color(0,0,0);
			Vector direction(x + j + 0.5 - W / 2, H - (y + i) - 1 + 0.5 - H / 2, -W / (2 * tan(alpha / 2)));
			direction.normalize();
			for (size_t k = 0; k < K; ++k){
				// boxMuller(1,x,y);
				// Vector direction(x + j + 0.5 - W / 2, H - (y + i) - 1 + 0.5 - H / 2, -W / (2 * tan(alpha / 2)));
				// direction.normalize();
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
