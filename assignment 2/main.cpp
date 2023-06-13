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
#include <sstream>
#include <climits>
#include <thread>
#include <numeric>
#include <iterator>
#include <random>
#include <unistd.h>
#include <initializer_list>
#include <list>

#include "liblbfgs/lbfgs.c"

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform (0 ,1);

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


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


// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name

double triangle_area(Vector x, Vector y, Vector z){
    Vector XY = y - x;
    Vector XZ = z - x;
    return std::abs(XY[0] * XZ[1] - XY[1] * XZ[0]) / 2;
        
}


class Polygon {  
public:
    std::vector<Vector> vertices;
    double cell_area(){
        double A = 0.;
        if (vertices.size() < 3) return 0;
        for (int i = 0; i < vertices.size(); i++){
            int j = ((i < vertices.size() - 1)?(i + 1):0);
            Vector x = vertices[i];
            Vector y = vertices[j];
            A += x[0] * y[1] - x[1] * y[0];
        }
        A /= 2;
        return std::abs(A);
    }


    double points_dist_integral(Vector& P) 
    {
		double res = 0.;
        int N = vertices.size();
        if (N < 3) return 0;
		for (int i = 1; i < N-1; i++) 
        {
			Vector C[3] = {vertices[0], vertices[i], vertices[i+1]};
            double T_area = triangle_area(C[0], C[1], C[2]);
			for (int k = 0; k < 3; k++) {
				for (int l = k; l < 3; l++) {
					res += T_area / 6 * dot(C[k]-P, C[l]-P);
				}

			}
		}
		
		return res;
	}

    Vector centroid() 
    {
		Vector C(0, 0, 0);
		double A = cell_area();
        int N = vertices.size();
		if (A == 0){
			if (N > 0) return vertices[0];
			else return Vector(0, 0, 0);
        }
		for (int i = 0; i < N - 1; i++) {
			C = C + 1/(6*A) * (vertices[i] + vertices[i+1]) * (vertices[i][0] * vertices[i + 1][1] - vertices[i + 1][0] * vertices[i][1]);
		}
        C = C + 1/(6*A) * (vertices[N - 1] + vertices[0]) * (vertices[N - 1][0] * vertices[0][1] - vertices[0][0] * vertices[N - 1][1]);
		return -C;
	}
};  
 
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+"); 
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}
 
 
// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
    FILE* f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    //if (i < N) {   // the N first particles may represent fluid, displayed in blue
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    //}
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }

                }
                
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}


/*
CLIPPING
*/


Polygon clip_edge(Polygon& subjectPolygon, const Vector& u, const Vector& v){
    Polygon outPolygon;
    int k = subjectPolygon.vertices.size();
    Vector N(v[1] - u[1], u[0] - v[0]);
    for (int i = 0; i < k; i++){
        Vector& curVertex = subjectPolygon.vertices[i];
        Vector& prevVertex =subjectPolygon.vertices[(i > 0)?(i - 1):(k-1)];
        double t = dot(u - prevVertex, N) / dot(curVertex - prevVertex, N);
        Vector P = prevVertex + t * (curVertex - prevVertex);

        if(dot(u - curVertex, N) <= 0){
            if (!(dot(u - prevVertex, N) <= 0)){
                outPolygon.vertices.push_back(P);
            }
            outPolygon.vertices.push_back(curVertex);

        }else if(dot(u - prevVertex, N) <= 0){
            outPolygon.vertices.push_back(P);

        }
    }
    return outPolygon;
}

Polygon voronoi_clipping(Polygon& subjectPolygon, int i, int j, const Vector* P, const double* w){
    Polygon outPolygon;
    int k = subjectPolygon.vertices.size();
    for (int a = 0; a < k; a++){
        Vector& curVertex = subjectPolygon.vertices[a];
        Vector& prevVertex = subjectPolygon.vertices[(a > 0)?(a - 1):(subjectPolygon.vertices.size()-1)];
        Vector M = (P[i] + P[j]) / 2;
        Vector diff = P[j] - P[i];
        M = M + (w[i] - w[j])/(2 * (-diff).norm2()) * (diff);
        double t = dot(M - prevVertex, diff) / dot(curVertex - prevVertex, diff);
        Vector inter = prevVertex + t * (curVertex - prevVertex);

        if((curVertex - P[i]).norm2() - w[i] <= (curVertex - P[j]).norm2() - w[j]){
            if (!((prevVertex - P[i]).norm2() - w[i] <= (prevVertex - P[j]).norm2() - w[j])){
                outPolygon.vertices.push_back(inter);
            }
            outPolygon.vertices.push_back(curVertex);

        }else if((prevVertex - P[i]).norm2() - w[i] <= (prevVertex - P[j]).norm2() - w[j]){
            outPolygon.vertices.push_back(inter);

        }
    }
    return outPolygon;
}


/*
DIAGRAM GENERATION
*/

std::vector<Polygon> build_diagram(const Vector* P, const double* w, int N){
    std::vector<Polygon> res(N);
	for (int i = 0; i < N; i++){
        res[i].vertices.push_back(Vector(0, 0, 0));
        res[i].vertices.push_back(Vector(0, 1, 0));
        res[i].vertices.push_back(Vector(1, 1, 0));
        res[i].vertices.push_back(Vector(1, 0, 0));
        for (int j = 0; j < N; j++){
            if (i == j) continue;
            res[i] = voronoi_clipping(res[i], i, j, P, w);
        }
	}
    return res;
}

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int N,
    const lbfgsfloatval_t step
    )
{
    int i;
    lbfgsfloatval_t fx = 0.0;

	Vector* P = static_cast<Vector*>(instance);
	std::vector<Polygon> diagram = build_diagram(&P[0], &x[0],N);
	double lambda = 1. / N;
    for (i = 0; i < N; i++ ) 
    {
		double poly_area = diagram[i].cell_area();
        g[i] = poly_area - lambda;
		fx += -diagram[i].points_dist_integral(P[i]) - x[i] * lambda + x[i] * poly_area;

	}
    return fx;
}

/*
FLUID DIAGRAM
*/


std::vector<Polygon> build_fluid(const Vector* P, const double* w, int N)
{
    std::vector<Polygon> diagram(N - 1);   
    int M = 200;
    #pragma omp parallel for 
	for (int i = 0; i < N - 1; i++){
        double r = sqrt(w[i] - w[N-1]);
        Polygon fluid_cells;
        fluid_cells.vertices.resize(M);
        for (int k = 0; k < M; k++){
            fluid_cells.vertices[k] = Vector(cos(k/ (double)M * 2 * M_PI), -sin(k/(double)M * 2 * M_PI)) * r + P[i];
        }
        diagram[i].vertices.push_back(Vector(0, 0, 0));
        diagram[i].vertices.push_back(Vector(0, 1, 0));
        diagram[i].vertices.push_back(Vector(1, 1, 0));
        diagram[i].vertices.push_back(Vector(1, 0, 0));
        for (int j = 0; j < N; j++){
            if (i == j) continue;
            diagram[i] = voronoi_clipping(diagram[i], i, j, P, w);
        }
        for (int l = 0; l < M - 1; l++){
            diagram[i] = clip_edge(diagram[i], fluid_cells.vertices[l], fluid_cells.vertices[l + 1]);
        }
        diagram[i] = clip_edge(diagram[i], fluid_cells.vertices[M - 1], fluid_cells.vertices[0]);
	}
    return diagram;
}

static lbfgsfloatval_t fluid_evaluate(

    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int N,
    const lbfgsfloatval_t step
    )
{
    lbfgsfloatval_t fx = 0.0;

	Vector* P = static_cast<Vector*>(instance);
	std::vector<Polygon> diagram = build_fluid(&P[0], &x[0], N);
    double total_fluid_area = 0;
    double fluid = 0.4;
    double air = 0.6;
	double lambda = fluid / (N - 1);
    for (int i = 0; i < N - 1; i++) 
    {
		double poly_area = diagram[i].cell_area();
        total_fluid_area += poly_area;
        g[i] = poly_area - lambda;
		fx += -diagram[i].points_dist_integral(P[i]) - x[i] * lambda + x[i] * poly_area;
	}
    g[N - 1] = 1 - total_fluid_area - air;
	fx +=  x[N - 1] * (1 - total_fluid_area) - x[N - 1] * air;
    return fx;
}

/*
RENDERING
*/


void gallouet_time_step(std::vector<Vector> &P, std::vector<Vector> &v, std::vector<double> &w, int t)
{
	double m = 200;
	double eps = 0.004;
	double dt = 0.002;

    int N = P.size();
    double fx;

	int ret = lbfgs(N + 1, &w[0], &fx, fluid_evaluate, NULL, &P[0], NULL);
    std::cout << "computed" << std::endl;
    std::cout<<ret<<std::endl;
    std::vector<Polygon> diagram = build_fluid(&P[0], &w[0], N + 1);
	save_frame(diagram, "frames/frame_", t);
	
    Vector gravity = Vector(0, -9.81);
    double eps2 = (std::pow(eps,2));
    for (int i = 0; i < N; i++)
	{
		Vector F_spring = 1 / eps2 * (diagram[i].centroid() - P[i]);
		Vector F = F_spring + m * gravity;
		v[i] = v[i] + dt / m * F;
		P[i] = P[i] + dt * v[i];
		P[i][0] = (P[i][0] < 0) ? -P[i][0] : P[i][0];
        P[i][1] = (P[i][1] < 0) ? -P[i][1] : P[i][1];
        P[i][0] = (P[i][0] >= 1) ? 2 - P[i][0] : P[i][0];
        P[i][1] = (P[i][1] >= 1) ? 2 - P[i][1] : P[i][1];
	}
}


////////////////////////
/*        MAIN        */
////////////////////////


int main()
{
	int N = 200;
	std::vector<Vector> P(N);
	std::vector<Vector> v(N);
	std::vector<double> w(N, 0);
	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < 2; j++)
		{
			P[i][j] = rand() / (double)RAND_MAX;
			v[i][j] = 0;
		}
		v[i][2] = 0;
		P[i][2] = 0;
		w[i] = 0.05;
	}
	for (int t = 0; t < 200; t++) {
        std::cout << t << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
		gallouet_time_step(P, v, w, t);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        double elapsed_secs = elapsed.count();

        std::cout << "Elapsed time: " << elapsed_secs << " seconds" << std::endl;
	}
	return 0;
}