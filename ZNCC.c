#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lodepng.h>
#include <limits.h>
#include <stdint.h>
#include <time.h>


unsigned error;
uint32_t width;
uint32_t height;
uint8_t* final_output;

uint32_t Width;
uint32_t Height;
unsigned int gray = 0;
int DISP = 0;
int win_x = 10;
int win_y = 17;
int win_size = 170;
uint8_t* LeftImage = NULL;
uint8_t* RightImage = NULL;
uint8_t* final_output;

uint8_t* image_input_vector2 = NULL;
uint8_t* resized_input_vector = NULL;
uint8_t* ResizedLeft;
uint8_t* ResizedRight;
uint8_t* resized_left_input_vector = NULL;
uint8_t* resized_right_input_vector = NULL;
uint8_t* left_mean_vector = NULL;
uint8_t* right_mean_vector = NULL;
int scale = 4;
unsigned intial_i, initial_j;
uint8_t* image_output_vector;
uint8_t* image_output_vector2 = NULL;

uint8_t* resized_Left;
uint8_t* resized_Right;
uint8_t result_image;
int32_t y;
int32_t x;
float zncc_value;
float den = 1;
int32_t current_sum;
int32_t disparity_value;
float upper_sum;
float lower_sum0;
float lower_sum1;
float left_value_diff_from_average;
float right_value_diff_from_average;
int32_t Right_average;
int32_t Left_average;
float best_zncc;
uint32_t best_disp;
uint32_t index = 0;
uint8_t ws = 7;
uint32_t WindowSum_L = 0;
uint32_t WindowSum_R = 0;
uint32_t left_mean;
uint32_t right_mean;
uint32_t lower_sum_0 = 0;
uint32_t lower_sum_1 = 0;
uint32_t upper_sum = 0;
uint32_t lower_sum_0_2 = 0;
uint32_t lower_sum_1_2 = 0;
uint32_t upper_sum_2 = 0;
uint8_t MAX_DISP = 32;
uint32_t ZNCC;
uint32_t ZNCC_2;
uint32_t MAX_ZNCC = 0;
uint32_t BEST_DISP = 0;
uint32_t MAX_ZNCC_2 = 0;
uint32_t BEST_DISP_2 = 0;
uint8_t* First_Dmap;
uint8_t* Second_Dmap;
int d, dd;

int32_t window_size = 25;

int32_t small_w, small_h;

uint8_t* map;
uint32_t threshold;
clock_t tt;

void GrayResize() {
    ResizedLeft = (uint8_t*)malloc((width * height) / 8);
    ResizedRight = (uint8_t*)malloc((width * height) / 8);

    small_w = width / 4;
    small_h = height / 4;
    for (uint32_t i = 0; i < small_h; i++) {
        for (uint32_t j = 0; j < small_w; j++) {

            ResizedLeft[i * small_w + j] =
                (LeftImage[(4 * i - 1 * (i > 0)) * (4 * width) + 4 * (4 * j - 1 * (j > 0))]
                    + LeftImage[(4 * i - 1 * (i > 0)) * (4 * width) + 4 * (4 * j - 1 * (j > 0)) + 1]
                    + LeftImage[(4 * i - 1 * (i > 0)) * (4 * width) + 4 * (4 * j - 1 * (j > 0)) + 2]) / 3;

            ResizedRight[i * small_w + j] =
                (RightImage[(4 * i - 1 * (i > 0)) * (4 * width) + 4 * (4 * j - 1 * (j > 0))]
                    + RightImage[(4 * i - 1 * (i > 0)) * (4 * width) + 4 * (4 * j - 1 * (j > 0)) + 1]
                    + RightImage[(4 * i - 1 * (i > 0)) * (4 * width) + 4 * (4 * j - 1 * (j > 0)) + 2]) / 3;
        }
    }

}
void zncc() {

    First_Dmap = (uint8_t*)malloc((Width * Height));
    Second_Dmap = (uint8_t*)malloc((Width * Height));



    for (int h = 0; h < (Height - ws); h++) {
        for (int w = 0; w < (Width); w++) {
            index = (Width * h) + w;

            BEST_DISP = 0;
            MAX_ZNCC = 0;
            BEST_DISP_2 = 0;
            MAX_ZNCC_2 = 0;

            for (d = 0; d < MAX_DISP; d++) {

                for (int y = 0; y < ws; y++) {
                    for (int x = 0;x < ws;x++) {
                        WindowSum_L = WindowSum_L + resized_left_input_vector[index + (Width * y) + x];
                        WindowSum_R = WindowSum_R + resized_right_input_vector[index + (Width * y) + x];
                    }
                }
                left_mean = WindowSum_L / (ws * ws);
                right_mean = WindowSum_R / (ws * ws);



                WindowSum_L = 0;
                WindowSum_R = 0;

                lower_sum_0 = 0;
                lower_sum_1 = 0;
                upper_sum = 0;
                lower_sum_0_2 = 0;
                lower_sum_1_2 = 0;
                upper_sum_2 = 0;
                for (int i = 0; i < ws; i++) {
                    for (int j = 0; j < ws; j++) {
                        lower_sum_0 += (left_mean - resized_left_input_vector[index + (Width * i) + j] * (left_mean - resized_left_input_vector[index + (Width * i) + j]));
                        lower_sum_1 += (right_mean - resized_right_input_vector[index + (Width * i) + j + d]) * (right_mean - resized_right_input_vector[index + (Width * i) + j + d]);
                        upper_sum += resized_left_input_vector[index + (Width * i) + j] * resized_right_input_vector[index + (Width * i) + j + d];

                        if (index < MAX_DISP) {
                            dd = index;
                            lower_sum_0_2 += (left_mean - resized_left_input_vector[index + (Width * i) + j - dd] * (left_mean - resized_left_input_vector[index + (Width * i) + j - dd]));
                            lower_sum_1_2 += (right_mean - resized_right_input_vector[index + (Width * i) + j]) * (right_mean - resized_right_input_vector[index + (Width * i) + j]);
                            upper_sum_2 += resized_left_input_vector[index + (Width * i) + j - dd] * resized_right_input_vector[index + (Width * i) + j];

                        }
                        else {
                            lower_sum_0_2 += (left_mean - resized_left_input_vector[index + (Width * i) + j - d] * (left_mean - resized_left_input_vector[index + (Width * i) + j - d]));
                            lower_sum_1_2 += (right_mean - resized_right_input_vector[index + (Width * i) + j]) * (right_mean - resized_right_input_vector[index + (Width * i) + j]);
                            upper_sum_2 += resized_left_input_vector[index + (Width * i) + j - d] * resized_right_input_vector[index + (Width * i) + j];

                        }
                    }
                }

                ZNCC = upper_sum / (sqrt(lower_sum_0) * sqrt(lower_sum_1));
                ZNCC_2 = upper_sum_2 / (sqrt(lower_sum_0_2) * sqrt(lower_sum_1_2));



                if (ZNCC > MAX_ZNCC) {
                    MAX_ZNCC = ZNCC;
                    BEST_DISP = d;
                }

                if (ZNCC_2 > MAX_ZNCC_2) {
                    MAX_ZNCC_2 = ZNCC_2;
                    BEST_DISP_2 = d;
                }


            }

            First_Dmap[index] = (BEST_DISP * 255) / MAX_DISP;
            Second_Dmap[index] = (BEST_DISP_2 * 255) / MAX_DISP;
        }
    }



}
void CrossCheck() {

    map = (uint8_t*)malloc((Width) * (Height));
    threshold = 8;

    for (index = 0;index < (Width * Height);index++) {
        if (abs((First_Dmap[index] - Second_Dmap[index])) < threshold)
            map[index] = 0;
        else
            map[index] = First_Dmap[index];
    }


}
void OcclusionFilling() {
    for (int h = 0; h < (Height - ws); h++) {
        for (int w = 0; w < (Width); w++) {
            index = (Width * h) + w;
            if (map[index] == 0) {
                map[index] = ((map[index + 1] + map[index + 2] + map[index + 3] + map[index + 4]) / 4);
            }
        }

    }
}

void occlusion_filling() {
    final_output = (uint8_t*)malloc((Width) * (Height));
    int index;  
    for (int h = 0; h < (Height - ws); h++) {
        for (int w = 0; w < (Width); w++) {
            index = (Width * h) + w;
            if (map[index] != 0)
                final_output[index] = map[index];
            else {
                if (map[index + 1] != 0)
                    final_output[index] = map[index + 1];


            }


        }
    }
}
int main(int argc, char* argv[])
{




    uint32_t Error;
    uint32_t w1, h1;
    uint32_t w2, h2;
    uint8_t* resizedL;
    const char* filename = "OriginalImageR.png";   
    const char* filename2 = "OriginalImageL.png";   
    const char* filename3;



    error = lodepng_decode32_file(&LeftImage, &width, &height, filename);
    if (error) printf("error loading image 1 %u: %s\n", error, lodepng_error_text(error));
    error = lodepng_decode32_file(&RightImage, &width, &height, filename2);
    if (error) printf("error loading image 2 %u: %s\n", error, lodepng_error_text(error));




    GrayResize();





    filename = "GrayscaledR.png";
    error = lodepng_encode_file(filename, ResizedLeft, small_w, small_h, LCT_GREY, 8);

    filename3 = "GrayscaledL.png";
    error = lodepng_encode_file(filename3, ResizedRight, small_w, small_h, LCT_GREY, 8);




    error = lodepng_decode_file(&resized_left_input_vector, &Width, &Height, filename, LCT_GREY, 8);
    if (error) printf("error encoding gray image 1 %u: %s\n", error, lodepng_error_text(error));
    error = lodepng_decode_file(&resized_right_input_vector, &Width, &Height, filename3, LCT_GREY, 8);
    if (error) printf("error loading gray image 2 %u: %s\n", error, lodepng_error_text(error));

    


    zncc();


    filename = "DmapedR.png";
    error = lodepng_encode_file(filename, First_Dmap, Width, Height, LCT_GREY, 8);

    filename = "DmapedL.png";
    error = lodepng_encode_file(filename, Second_Dmap, Width, Height, LCT_GREY, 8);


    
    free(LeftImage);
    free(RightImage);


    CrossCheck();




    filename = "Crosscheck.png";
    error = lodepng_encode_file(filename, map, Width, Height, LCT_GREY, 8);





    //OcclusionFilling();
    occlusion_filling();

  
    

    filename = "OcclusionFilled.png";
    error = lodepng_encode_file(filename, map, Width, Height, LCT_GREY, 8);


    free(First_Dmap);
    free(Second_Dmap);
    free(map);

    free(resized_left_input_vector);
    free(resized_right_input_vector);

    return 0;


}
