#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <algorithm>

static std::default_random_engine engine(10);
static std::normal_distribution<double> normal(0, 1);

std::array<double, 3> randomDirectionGenerator() {
	double xDir = normal(engine);
	double yDir = normal(engine);
	double zDir = normal(engine);
	double norm = sqrt(xDir * xDir + yDir * yDir + zDir * zDir);
	xDir /= norm;
	yDir /= norm;
	zDir /= norm;
	return {xDir, yDir, zDir};
}


void modifiedSlicedMatching(const std::vector<double>& sourceImage, const std::vector<double>& sourceColor, int width, int height, std::vector<double>& resultImage) {
// Create a copy of the source image
	resultImage = sourceImage;

	for(int iteration = 0; iteration < 100; ++iteration) {
		auto [xDirection, yDirection, zDirection] =  randomDirectionGenerator();
		std::vector<std::pair<double, int> > sortedResults(width * height);
		std::vector<double> sortedSourceColor(width * height);

		for(int i = 0; i < width * height; ++i) {
			double dotProductResult = resultImage[i * 3] * xDirection + resultImage[i * 3 + 1] * yDirection + resultImage[i * 3 + 2] * zDirection;
			double dotProductSource = sourceColor[i * 3] * xDirection + sourceColor[i * 3 + 1] * yDirection + sourceColor[i * 3 + 2] * zDirection;
			sortedResults[i] = std::make_pair(dotProductResult, i);
			sortedSourceColor[i] = dotProductSource;
		}

		std::sort(sortedResults.begin(), sortedResults.end());
		std::sort(sortedSourceColor.begin(), sortedSourceColor.end());

		for(int i = 0; i < width * height; ++i) {
			double motionAmount = sortedSourceColor[i] - sortedResults[i].first;
			int index = sortedResults[i].second;
			resultImage[index * 3 + 0] += motionAmount * xDirection;
			resultImage[index * 3 + 1] += motionAmount * yDirection;
			resultImage[index * 3 + 2] += motionAmount * zDirection;
		}	
	}
}

unsigned char clamp(double x, double min, double max) {
	if(x < min) return min;
	if(x > max) return max;
	return x;
}

int main() {
	int W, H, C;
	
	// stbi_set_flip_vertically_on_load(false);
	unsigned char *image_source_ptr = stbi_load("img1.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);
	unsigned char *color_source_ptr = stbi_load("img2.jpg",
								 &W,
								 &H,
								 &C,
								 STBI_rgb);

	size_t total_size = W*H*3;
	std::vector<double> image_source(total_size);
	std::vector<double> color_source(total_size);
	for (int i = 0; i < total_size; ++i) {
		image_source[i] = (double)image_source_ptr[i];
		color_source[i] = (double)color_source_ptr[i];
	}

	std::vector<double> result(total_size);
	modifiedSlicedMatching(image_source, color_source, W, H, result);
	
	unsigned char *image_result = new unsigned char[W*H*3];
	for (int i = 0; i < H; ++i) {
		for (int j = 0; j < W; ++j) {
			image_result[(i*W + j) * 3] = clamp(result[(i*W + j) * 3], 0, 255);
			image_result[(i*W + j) * 3 + 1] = clamp(result[(i*W + j) * 3 + 1], 0, 255);
			image_result[(i*W + j) * 3 + 2] = clamp(result[(i*W + j) * 3 + 2], 0, 255);
		}
	}
	stbi_write_png("out.png", W, H, 3, image_result, 0);

	return 0;
}