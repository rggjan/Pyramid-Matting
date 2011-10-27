#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <limits.h>

using namespace std;

const char* result_name = NULL;

void load_image (const char* filename, int* ret_width, int* ret_height, int* ret_num_colors, unsigned char** data_char, double** data_double) {
  FILE *fp = fopen (filename, "rb");
  
  char buffer[100];
  int num_colors = -1;
  int width = -1, height = -1;
  int max_value = -1;

  while(num_colors == -1 || width == -1 || height == -1 || max_value == -1) {
    if (fgets(buffer, 100, fp) == NULL) {
      cout << "Could not read file " << filename << endl;
      exit(1);
    }

    if (buffer[0] == '#')
      continue;

    if (num_colors == -1) {
      if (strncmp(buffer, "P6", 2) == 0) {
        num_colors = 3;
      } else if (strncmp(buffer, "P5", 2) == 0) {
        num_colors = 1;
      } else {
        cout << "Unknown format " << buffer << endl;
        exit(1);
      }
    } else if (width == -1) {
      if (sscanf(buffer, "%i %i", &width, &height) != 2) {
        cout << "Could not read width/height" << endl;
        exit(1);
      }
    } else {
      if (sscanf(buffer, "%i", &max_value) != 1) {
        cout << "Could not read max_value" << endl;
        exit(1);
      }
    }
  }

  unsigned char* data = new unsigned char[width*height*num_colors];
  
  fread(data, 1, width*height*num_colors, fp);
  fclose (fp);

  if (data_double != NULL) {
    *data_double = new double[width*height*num_colors];
    double max_value_d = max_value;

    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int c=0; c<num_colors; c++) {
          (*data_double)[(y*width+x)*num_colors+c] = 
            data[(y*width+x)*num_colors+c]/max_value_d;
        }
      }
    }
  }

  if (data_char == NULL)
    delete[] data_char;
  else
    *data_char = data;

  *ret_num_colors = num_colors;
  *ret_width = width;
  *ret_height = height;
}

void
save_image (const char* filename, int dimx, int dimy, int num_colors,
    double* data) {

  unsigned char* char_data = new unsigned char[dimx*dimy*num_colors];

  for (int y=0; y<dimy; y++) {
    for (int x=0; x<dimx; x++) {
      for (int c=0; c<num_colors; c++) {
        double d = data[(y*dimx+x)*num_colors+c]*255;
        if (d>255)
          d = 255;
        if (d<0)
          d = 0;

        if (isnan(d)) {
          cout << "NaN in image " << filename << " at position x/y: " << x << "/" << y << endl;
          exit(1);
        }
        char_data[(y*dimx+x)*num_colors+c] = d;
      }
    }
  }

  FILE *fp = fopen (filename, "wb");
  if (fp == NULL) {
    cout << "Could not open file " << endl;
    exit(0);
  }
  
  if (num_colors == 3)
    {
      fprintf (fp, "P6\n%d %d\n255\n", dimx, dimy);
    }
  else if (num_colors == 1)
    {
      fprintf (fp, "P5\n%d %d\n255\n", dimx, dimy);
    }
  else
    {
      printf ("Problem!\n");
      exit (1);
    }

  fwrite(char_data, 1, dimx*dimy*num_colors, fp);
  fclose (fp);

  delete[] char_data;
}

static void save_result(const char* name, double* data, int raise, int num_colors, int width, int height) {
  char buffer[100];
  snprintf(buffer, 100, "%s/%s_%i.ppm", result_name, name, raise);
  save_image(buffer, width, height, num_colors, data);
}

double projection (const double F[3],
                   const double B[3],
                   const double C[3],
                   double* alpha_pointer) {
  double BF[3];
  double BC[3];
  double dot_BF_BF = 0;
  double dot_BF_BC = 0;

  for (int c=0; c<3; c++) {
    BF[c] = F[c] - B[c];
    BC[c] = C[c] - B[c];
    dot_BF_BF += BF[c]*BF[c];
    dot_BF_BC += BF[c]*BC[c];
  }

  double alpha = dot_BF_BC/dot_BF_BF;

  if (alpha<0) 
    alpha = 0;

  if (alpha>1)
    alpha=1;

  double CC[3];
  for (int c=0; c<3; c++) {
    CC[c] = C[c] - (B[c] + alpha * BF[c]);
  }

  double quality = 0;
  for (int c=0; c<3; c++) {
    double diff = CC[c];
    quality += diff*diff;
  }

  *alpha_pointer = alpha;
  return quality;
}

int main(int argc, char* argv[]) {
  char* trimap_name = NULL;
  char* image_name = NULL;
  char* gt = NULL;

  bool debug = false;
  for (int i=1; i<argc; i++)
    if (strcmp(argv[i], "--debug") == 0)
      debug = true;
    else if (strcmp(argv[i], "--gt") == 0)
      gt = argv[++i];
    else if (image_name == NULL)
      image_name = argv[i];
    else if (trimap_name == NULL)
      trimap_name = argv[i];
    else if (result_name == NULL)
      result_name = argv[i];
    else {
      cout << "Usage: matting image.ppm trimap.pnm folder [--debug] [--gt ground_truth.pnm]\n";
      exit(0);
    }

  // ensure the "result" directory exists
  struct stat result_dir;
  if (stat(result_name, &result_dir) < 0) {
    int err = errno;
    if (err == ENOENT) {
      // create directory
      mkdir(result_name, 0777);
    } else {
      cerr << "Error while testing existence of 'results' directory: " << strerror(err) << '\n';
      return 1;
    }
  } else if (!S_ISDIR(result_dir.st_mode)) {
    cerr << "The file 'results' is not a directory.\n";
    return 1;
  }

  int original_width;
  int original_height;
  int num_colors;
  double* mask;

  load_image(trimap_name, &original_width, &original_height, &num_colors, NULL, &mask);
  if (num_colors != 1) {
    cerr << "num_colors != 1" << endl;
    exit(1);
  }

  int new_width, new_height;
  double* original;
  load_image(image_name, &new_width, &new_height, &num_colors, NULL, &original);
  if (num_colors != 3 || original_width != new_width || original_height != new_height) {
    cout << "original wrong" << endl;
  }

  int global_raise = original_width>original_height?
    original_width:original_height;
  global_raise = ceil(log2(global_raise));
  cout << "Width: " << original_width << endl;
  cout << "Height: " << original_height << endl;
  cout << "Raise: " << global_raise << endl;
  cout << "=============" << endl;

  int raise = global_raise;

  int width = original_width;
  int height = original_height;

  double* foreground = new double[width*height*3]();
  double* portion_foreground = new double[width*height]();
  double* background = new double[width*height*3]();
  double* portion_background = new double[width*height]();

  double* original_list[raise+1];
  double* mask_list[raise+1];
  double* color_list[raise+1][2];
  double* portion_list[raise+1][2];
  double* final_list[raise+1][2];
  double* alpha_list[raise+1];
  int width_list[raise+1];
  int height_list[raise+1];

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      int id=y*width+x;
      if (mask[id] == 1) {
        portion_foreground[id] = 1;
        for (int c=0; c<3; c++) {
          int idc = (width*y+x)*3+c;
          foreground[idc] = original[idc];
        }
      } else if (mask[id] == 0) {
        portion_background[id] = 1;
        for (int c=0; c<3; c++) {
          int idc = (width*y+x)*3+c;
          background[idc] = original[idc];
        }
      }
    }
  }
  
  original_list[raise] = original;
  color_list[raise][0] = foreground;
  portion_list[raise][0] = portion_foreground;
  color_list[raise][1] = background;
  portion_list[raise][1] = portion_background;
  mask_list[raise] = mask;
  width_list[raise] = width;
  height_list[raise] = height;

  if (!debug) {
    save_result("originals", original, raise, 3, width, height);
    save_result("foregrounds", color_list[raise][0], raise, 3, width, height);
    save_result("backgrounds", color_list[raise][1], raise, 3, width, height);
    save_result("masks", mask_list[raise], raise, 1, width, height);
  }

  while (raise > 0) {
    raise--;

    // Height half
    int old_width = width;
    int old_height = height;
    width = (width+1)/2;
    height = (height+1)/2;

    original_list[raise] = new double[width*height*3]();
    mask_list[raise] = new double[width*height]();
    for (int b=0; b<2; b++) {
      color_list[raise][b] = new double[width*height*3]();
      portion_list[raise][b] =  new double[width*height]();
    }
    height_list[raise] = height;
    width_list[raise] = width;

    double* mask = mask_list[raise];
    double* old_mask = mask_list[raise+1];
    double* original = original_list[raise];
    double* old_original = original_list[raise+1];
    double* color[2] = {color_list[raise][0], color_list[raise][1]};
    double* old_color[2] = {color_list[raise+1][0], color_list[raise+1][1]};
    double* portion[2] = {portion_list[raise][0], portion_list[raise][1]};
    double* old_portion[2] = {portion_list[raise+1][0],
                              portion_list[raise+1][1]};
    #pragma omp parallel for
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        int id = y*width+x;

        // Count number of available pixels
        int counter = 0;
        for (int n=0; n<4; n++) {
          int old_y = 2*y+(n&1);
          int old_x = 2*x+(n/2);
          if (old_x < old_width && old_y < old_height)
            counter++;
        }

        if (counter == 0) {
          cout << "This should not happen!" << endl;
          exit(1);
        }

        for (int n=0; n<4; n++) {
          int old_y = 2*y+(n&1);
          int old_x = 2*x+(n/2);
          if (old_x >= old_width || old_y >= old_height)
            continue;

          int idn = (old_y)*(old_width)+old_x;

          if (n == 0)
            mask[id] = old_mask[idn];
          else
            if (mask[id] != old_mask[idn])
              mask[id] = 0.5;

          for (int b=0; b<2; b++) {
            portion[b][id] += old_portion[b][idn]/counter;
          }
        }

        int id3 = (y*width+x)*3;
        for (int n=0; n<4; n++) {
          int old_y = 2*y+(n&1);
          int old_x = 2*x+(n/2);
          if (old_x >= old_width || old_y >= old_height)
            continue;

          int idn = old_y*old_width+old_x;
          int idn3 = idn*3;
          for (int c=0; c<3; c++) {
            original[id3+c] +=old_original[idn3+c]/counter;

            for (int b=0; b<2; b++) {
              double cur_portion = portion[b][id];
              if (cur_portion > 0)
                color[b][id3+c] +=
                  old_color[b][idn3+c] * old_portion[b][idn]
                  / cur_portion / counter;
            }
          }
        }
      }
    }

    if (!debug) {
      save_result("originals", original, raise, 3, width, height);
      save_result("foregrounds", color[0], raise, 3, width, height);
      save_result("backgrounds", color[1], raise, 3, width, height);
      save_result("masks", mask, raise, 1, width, height);
    }
  }
  
  // Upscaling
  raise = 0;
  width = 1;
  height = 1;

  alpha_list[0] = new double[1]();
  for (int b=0; b<2; b++) {
    final_list[0][b] = new double[3]();
    for (int c=0; c<3; c++) {
      final_list[0][b][c] = color_list[0][b][c];
    }
  }

  // Calculate alpha
  double merged_point[3];
  for (int c=0; c<3; c++) {
    double portion_foreground = portion_list[0][0][0];
    double portion_background = portion_list[0][1][0];
    merged_point[c] = (original_list[0][c]
      -portion_foreground*color_list[0][0][c]
      -portion_background*color_list[0][1][c])
      /(1-portion_foreground-portion_background);
  }

  projection(final_list[0][0], final_list[0][1],
      merged_point, &(alpha_list[0][0]));

  while (raise <= global_raise) {
    if (!debug || raise == global_raise) {
      double *tmp = new double[width*height*3]();
      for (int y=0; y<height; y++) {
        for (int x=0; x<width; x++) {
          int id = y*width+x;
          int id3 = id*3;
          double grey;
          if ((((int)(x/8))+((int)(y/8)))%2 == 0)
            grey = 153/255.;
          else
            grey = 102/255.;
          for (int c=0; c<3; c++) {
            double portion_part = portion_list[raise][0][id];
            double alpha_part = (1-portion_list[raise][0][id]-portion_list[raise][1][id])
              *alpha_list[raise][id];
            tmp[id3+c] = portion_part*color_list[raise][0][id3+c]
              +alpha_part*final_list[raise][0][id3+c]
              +(1-alpha_part-portion_part) * grey;
          }
        }
      }
      save_result("final", tmp, raise, 3, width, height);
      delete[] tmp;
      save_result("new_foregrounds", final_list[raise][0], raise, 3, width, height);
      save_result("new_backgrounds", final_list[raise][1], raise, 3, width, height);
      save_result("new_alphas", alpha_list[raise], raise, 1, width, height);
    }

    if (raise == global_raise)
      break;

    raise++;
    int old_width = width;
    int old_height = height;

    width = width_list[raise];
    height = height_list[raise];
    cout << "============= raise: " << raise << "=============" << endl;

    const double* original = original_list[raise];
    // const double* old_original = original_list[raise-1];

    const double* color[2] = {color_list[raise][0], color_list[raise][1]};
    // const double* old_color[2] = {color_list[raise-1][0], color_list[raise-1][1]};

    const double* portion[2] = {portion_list[raise][0], portion_list[raise][1]};
    // const double* old_portion[2] = {portion_list[raise-1][0],
    //   portion_list[raise-1][1]};

    const double* mask = mask_list[raise];
    const double* old_mask = mask_list[raise-1];

    alpha_list[raise] = new double[width*height]();
    for (int b=0; b<2; b++) {
      final_list[raise][b] = new double[width*height*3]();
    }

    double* alpha = alpha_list[raise];
    double* old_alpha = alpha_list[raise-1];
    
    double* final[2] = {final_list[raise][0], final_list[raise][1]};
    double* old_final[2] = {final_list[raise-1][0], final_list[raise-1][1]};

    #pragma omp parallel for
    for (int y=0; y<height; y+=2) {
      for (int x=0; x<width; x+=2) {
        int old_id = (y/2*old_width+x/2);
        int old_id3 = old_id*3;
        // int id = y*width+x;

        if (old_mask[old_id] == 1 || old_mask[old_id] == 0)
          continue;

        const double* best_color[4][2];
        double best_score[4] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX};

        double* test_color[2];

        int positions[] = {-16, -4, -1, 0, 1, 4, 16};
//        int positions[] = {-2, -1, 0, 1, 2};
        int num_pos = sizeof(positions)/sizeof(positions[0]);

        for (int bxdiff=0; bxdiff<num_pos; bxdiff++) {
          for (int bydiff=0; bydiff<num_pos; bydiff++) {
            // Check if inside old image
            if (!(x/2 + positions[bxdiff] >= 0 && y/2 + positions[bydiff] >= 0
                  && x/2 + positions[bxdiff] < old_width && y/2 + positions[bydiff] < old_height))
              continue;

            // Set test background
            int old_id_b = (y/2+positions[bydiff])*old_width+x/2+positions[bxdiff];
            int old_id3_b = old_id_b*3;
             
            if (old_mask[old_id_b] == 1 || old_mask[old_id_b] == 0)
              continue;

            test_color[1] = &(old_final[1][old_id3_b]);

            for (int fxdiff=0; fxdiff<num_pos; fxdiff++) {
              for (int fydiff=0; fydiff<num_pos; fydiff++) {
                // Check if inside old image
                if (!(x/2 + positions[fxdiff] >= 0 && y/2 + positions[fydiff] >= 0
                      && x/2 + positions[fxdiff] < old_width && y/2 + positions[fydiff] < old_height))
                  continue;

                // Set test foreground
                int old_id_f = (y/2+positions[fydiff])*old_width+x/2+positions[fxdiff];
                int old_id3_f = old_id_f*3;

                if (old_mask[old_id_f] == 1 || old_mask[old_id_f] == 0)
                  continue;

                test_color[0] = &(old_final[0][old_id3_f]);

                for (int n=0; n<4; n++) {
                  int ny = y+(n&1);
                  int nx = x+n/2;
                  if (nx >= width || ny >= height)
                    continue;

                  int idn = ny*width+nx;
                  if (mask[idn] == 1 || mask[idn] == 0)
                    continue;
                  int idn3 = idn*3;

                  double test_merged[2][3];
                  double proj_a;
                  double sum_portion = portion[0][idn] + portion[1][idn];

                  for (int c=0; c<3; c++) {
                    for (int b=0; b<2; b++){
                      double ratio = portion[b][idn] / (portion[b][idn]+(1-sum_portion)/2);
                      test_merged[b][c] = ratio*color[b][idn3+c]
                        + (1-ratio)*test_color[b][c];
                    }
                  }

                  double diff = positions[bxdiff]*positions[bxdiff]
                    +positions[bydiff]*positions[bydiff]
                    +positions[fxdiff]*positions[fxdiff]
                    +positions[fydiff]*positions[fydiff];
                  diff = sqrt(diff);
                  diff /= 5000;

                  double score = projection(test_merged[0], test_merged[1],
                      &(original[idn3]), &proj_a)+diff;

                  if (score < best_score[n]) {
                    best_score[n] = score;
                    for (int b=0; b<2; b++) {
                      best_color[n][b] = test_color[b];
                    }
                  }
                }
              }
            }
          }
        }
        
        for (int n=0; n<4; n++) {
          int ny = y+(n&1);
          int nx = x+n/2;
          if (nx >= width || ny >= height)
            continue;

          int idn = ny*width+nx;

          if (mask[idn] == 1 || mask[idn] == 0)
            continue;

          int idn3 = idn*3;
          double sum_portion = portion[0][idn] + portion[1][idn];

          for (int b=0; b<2; b++) {
            for (int c=0; c<3; c++) {
              double ratio = portion[b][idn] / (portion[b][idn]+(1-sum_portion)/2);

              final[b][idn3+c] = ratio*color[b][idn3+c]+
                (1-ratio)*
                (best_color[n][b][c]*raise+old_final[b][old_id3+c]*(global_raise-raise))/global_raise;
                //(best_color[n][b][c]+old_final[b][old_id3+c])/2;
                //(best_color[n][b][c]+old_final[b][old_id3+c])/2;
            }
          }

          projection(&(final[0][idn3]),
                     &(final[1][idn3]),
                     &(original[idn3]),
                     &(alpha[idn]));

          //alpha[idn] = (alpha[idn]*raise + old_alpha[old_id]*(global_raise-raise))/global_raise;
        }
      }
    }
  }

  if (gt != NULL) {
    // Compare to ground truth, if we have one
    char buffer[100];
    snprintf(buffer, 100, "%s/new_alphas_%i.ppm", result_name, raise);
    
    unsigned char* result;
    unsigned char* ground_truth;

    int new_width;
    int new_height;
    int num_colors;
    load_image(gt, &new_width, &new_height, &num_colors, &ground_truth, NULL);
    if (num_colors != 1 || new_width != original_width || new_height != original_height) {
      cerr << "ground truth wrong!" << endl;
      exit(1);
    }

    load_image(buffer, &width, &height, &num_colors, &result, NULL);

    int sum = 0;
    for (int x=0; x<height*width; x++) {
      if (mask[x] != 0 && mask[x] != 1) {
        int diff = ground_truth[x] - result[x];
        sum += abs(diff);
      }
    }
    cout << "SAD: " << ((double)sum)/(width*height) << endl;
  }

}
