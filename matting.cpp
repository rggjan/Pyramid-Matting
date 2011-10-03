#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <limits.h>

using namespace std;

double*
load_image (const char* filename, int* ret_width, int* ret_height, int* ret_num_colors) {

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
  double* data_double = new double[width*height*num_colors];
  double max_value_d = max_value;


  fread(data, 1, width*height*num_colors, fp);
  fclose (fp);

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      for (int c=0; c<num_colors; c++) {
        data_double[(y*width+x)*num_colors+c] = 
          data[(y*width+x)*num_colors+c]/max_value_d;
      }
    }
  }

  delete[] data;

  *ret_num_colors = num_colors;
  *ret_width = width;
  *ret_height = height;
  return data_double;
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

#define RESULTS "results/"

int main() {
  // ensure the "result" directory exists
  struct stat result_dir;
  if (stat(RESULTS, &result_dir) < 0) {
    int err = errno;
    if (err == ENOENT) {
      // create directory
      mkdir(RESULTS, 0777);
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

  double* mask = load_image("trimap.pnm", &original_width, &original_height, &num_colors);
  if (num_colors != 1) {
    cerr << "num_colors != 1" << endl;
    exit(1);
  }

  int new_width, new_height;

  double* original = load_image("test.ppm", &new_width, &new_height, &num_colors);
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

  static char buffer[100];
  snprintf(buffer, 100, RESULTS "originals_%i.ppm", raise);
  save_image(buffer, width, height, 3, original_list[raise]);

  snprintf(buffer, 100, RESULTS "foregrounds_%i.ppm", raise);
  save_image(buffer, width, height, 3, color_list[raise][0]);

  snprintf(buffer, 100, RESULTS "backgrounds_%i.ppm", raise);
  save_image(buffer, width, height, 3, color_list[raise][1]);

  snprintf(buffer, 100, RESULTS "masks_%i.ppm", raise);
  save_image(buffer, width, height, 1, mask_list[raise]);

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

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "originals_%i.ppm", raise);
    save_image(buffer, width, height, 3, original);
    
    snprintf(buffer, 100, RESULTS "foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, color[0]);
    
    snprintf(buffer, 100, RESULTS "backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, color[1]);
  
    snprintf(buffer, 100, RESULTS "masks_%i.ppm", raise);
    save_image(buffer, width, height, 1, mask);
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

  save_image(RESULTS "new_foregrounds_0.ppm", 1, 1, 3, final_list[0][0]);
  save_image(RESULTS "new_backgrounds_0.ppm", 1, 1, 3, final_list[0][1]);
  save_image(RESULTS "new_alphas_0.ppm", 1, 1, 1, alpha_list[0]);


  double final[3];
  for (int c=0; c<3; c++) {
    final[c] = portion_list[0][0][0]*color_list[0][0][c]
      +final_list[0][0][c]*alpha_list[0][0]
      *(1-portion_list[0][0][0]-portion_list[0][1][0]);
  }
  save_image(RESULTS "final_0.ppm", 1, 1, 3, final);

  while (raise < global_raise) {
    raise++;
    int old_width = width;
    int old_height = height;

    width = width_list[raise];
    height = height_list[raise];
    cout << "============= raise: " << raise << "=============" << endl;

    const double* original = original_list[raise];
    const double* old_original = original_list[raise-1];

    const double* color[2] = {color_list[raise][0], color_list[raise][1]};
    const double* old_color[2] = {color_list[raise-1][0], color_list[raise-1][1]};

    const double* portion[2] = {portion_list[raise][0], portion_list[raise][1]};
    const double* old_portion[2] = {portion_list[raise-1][0],
      portion_list[raise-1][1]};

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

    for (int y=0; y<height; y+=2) {
      for (int x=0; x<width; x+=2) {
        int old_id = (y/2*old_width+x/2);
        int old_id3 = old_id*3;
        int id = y*width+x;

        if (old_mask[old_id] == 1 || old_mask[old_id] == 0)
          continue;

        const double* best_color[4][2];
        double best_score[4] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX};

        double* test_color[2];

        const int f_radius = 2;
        const int b_radius = f_radius;
        for (int bxdiff=-b_radius; bxdiff<=b_radius; bxdiff++) {
          for (int bydiff=-b_radius; bydiff<=b_radius; bydiff++) {
            // Check if inside old image
            if (!(x/2 + bxdiff >= 0 && y/2 + bydiff >= 0
                  && x/2 + bxdiff < old_width && y/2 + bydiff < old_height))
              continue;

            // Set test background
            int old_id_b = (y/2+bydiff)*old_width+x/2+bxdiff;
            int old_id3_b = old_id_b*3;
             
            if (old_mask[old_id_b] == 1 || old_mask[old_id_b] == 0)
              continue;

            test_color[1] = &(old_final[1][old_id3_b]);

            for (int fxdiff=-f_radius; fxdiff<=f_radius; fxdiff++) {
              for (int fydiff=-f_radius; fydiff<=f_radius; fydiff++) {
                // Check if inside old image
                if (!(x/2 + fxdiff >= 0 && y/2 + fydiff >= 0
                      && x/2 + fxdiff < old_width && y/2 + fydiff < old_height))
                  continue;

                // Set test foreground
                int old_id_f = (y/2+fydiff)*old_width+x/2+fxdiff;
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

                  double score = projection(test_merged[0], test_merged[1],
                      &(original[idn3]), &proj_a);

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
                best_color[n][b][c];
            }
          }


          projection(&(final[0][idn3]),
                     &(final[1][idn3]),
                     &(original[idn3]),
                     &(alpha[idn]));
        }
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "final_%i.ppm", raise);
    double *tmp = new double[width*height*3]();
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        int id = y*width+x;
        int id3 = id*3;
        for (int c=0; c<3; c++) {
          tmp[id3+c] = portion[0][id]*color[0][id3+c]
            +final[0][id3+c]*alpha[id]
            *(1-portion[0][id]-portion[1][id]);
        }
      }
    }
    save_image(buffer, width, height, 3, tmp);
    delete[] tmp;

    snprintf(buffer, 100, RESULTS "new_foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, final[0]);
    
    snprintf(buffer, 100, RESULTS "new_backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, final[1]);
    
    snprintf(buffer, 100, RESULTS "new_alphas_%i.ppm", raise);
    save_image(buffer, width, height, 1, alpha);
    
  }
/*
  // cleanup (not strictly needed, but makes valgrind output cleaner)
  delete[] mask;

  for (int i = 0; i < raise; i++) {
    for (int j = 0; j < 2; j++) {
      delete[] originals[i][j];
      delete[] foregrounds[i][j];
      delete[] foreground_ps[i][j];
      delete[] backgrounds[i][j];
      delete[] background_ps[i][j];
      delete[] new_foregrounds[i][j];
      delete[] new_backgrounds[i][j];
      delete[] new_alphas[i][j];
    }
  }

  delete[] originals[9][0];
  delete[] foregrounds[9][0];
  delete[] foreground_ps[9][0];
  delete[] backgrounds[9][0];
  delete[] background_ps[9][0];
  delete[] new_foregrounds[9][0];
  delete[] new_backgrounds[9][0];
  delete[] new_alphas[9][0];*/
}
