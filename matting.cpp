#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <limits.h>

using namespace std;

#define ID(x, y) ((y*width)+(x))

double*
load_image (const char* filename, int dimx, int dimy, int num_colors) {
  unsigned char* data = new unsigned char[dimx*dimy*num_colors];
  double* data_double = new double[dimx*dimy*num_colors];

  FILE *fp = fopen (filename, "rb");
  fread(data, 1, dimx*dimy*num_colors, fp);
  fclose (fp);

  for (int y=0; y<dimy; y++) {
    for (int x=0; x<dimx; x++) {
      for (int c=0; c<num_colors; c++) {
        data_double[(y*dimx+x)*num_colors+c] = data[(y*dimx+x)*num_colors+c]/255.;
      }
    }
  }

  delete[] data;
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

double projection (double F[3], double B[3], double C[3], double* alpha_pointer) {
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

  int width = 512;
  int height = 512;
  int raise = 9;

  double* mask = load_image("trimap.pnm", width, height, 1);

  double* original_list[raise+1];
  double* color_list[raise+1][2];
  double* portion_list[raise+1][2];
  double* final_list[raise+1][2];
  double* alpha_list[raise+1];

  double* original = load_image("test.ppm", width, height, 3);
  double* foreground = new double[width*height*3]();
  double* portion_foreground = new double[width*height]();
  double* background = new double[width*height*3]();
  double* portion_background = new double[width*height]();

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
  
  original_list[9] = original;
  color_list[9][0] = foreground;
  portion_list[9][0] = portion_foreground;
  color_list[9][1] = background;
  portion_list[9][1] = portion_background;

  save_image(RESULTS "originals_9.ppm", width, height, 3, original_list[9]);
  save_image(RESULTS "foregrounds_9.ppm", width, height, 3, color_list[9][0]);
  save_image(RESULTS "backgrounds_9.ppm", width, height, 3, color_list[9][1]);

  while (raise > 0) {
    raise--;

    // Height half
    height = height/2;
    width = width/2;

    original_list[raise] = new double[width*height*3]();
    for (int b=0; b<2; b++) {
      color_list[raise][b] = new double[width*height*3]();
      portion_list[raise][b] =  new double[width*height]();
    }
  
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int xdiff = 0; xdiff <= 1; xdiff++) {
          for (int ydiff = 0; ydiff <= 1; ydiff++) {
            for (int b=0; b<2; b++) {
              portion_list[raise][b][y*width+x] +=
                portion_list[raise+1][b][(2*y+ydiff)*(2*width)+2*x+xdiff]/4;
            }
          }
        }

        for (int xdiff = 0; xdiff <= 1; xdiff++) {
          for (int ydiff = 0; ydiff <= 1; ydiff++) {
            for (int c=0; c<3; c++) {
              original_list[raise][(y*width+x)*3+c] +=
                original_list[raise+1][((2*y+ydiff)*(2*width)+2*x+xdiff)*3+c]/4;

              for (int b=0; b<2; b++) {
                double new_ps = portion_list[raise][b][y*width+x];
                if (new_ps > 0)
                  color_list[raise][b][(y*width+x)*3+c] +=
                    color_list[raise+1][b][((2*y+ydiff)*(2*width)+2*x+xdiff)*3+c]
                    *portion_list[raise+1][b][(2*y+ydiff)*(2*width)+2*x+xdiff]
                    /new_ps/4;
              }
            }
          }
        }
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "originals_%i.ppm", raise);
    save_image(buffer, width, height, 3, original_list[raise]);
    
    snprintf(buffer, 100, RESULTS "foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, color_list[raise][0]);
    
    snprintf(buffer, 100, RESULTS "backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, color_list[raise][1]);
  }
  

  // Upscaling
  raise = 0;
  width = 1;
  height = 1;

  alpha_list[0] = new double[1];
  for (int b=0; b<2; b++) {
    final_list[0][b] = new double[3];
    for (int c=0; c<3; c++) {
      final_list[0][b][c] = color_list[0][b][c];
    }
  }

  // Calculate alpha
  double merged_point[3];
  for (int c=0; c<3; c++) {
    double pf = portion_list[0][0][0];
    double pb = portion_list[0][1][0];
    merged_point[c] = (original_list[0][c]
      -pf*color_list[0][0][c]
      -pb*color_list[0][1][c])
      /(1-pf-pb);
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

  while (raise < 9) {
    cout << "============= raise: " << raise << "=============" << endl;
    width *= 2;
    height *= 2;
    raise++;

    for (int b=0; b<2; b++) {
      final_list[raise][b] = new double[width*height*3];
    }
    alpha_list[raise] = new double[width*height];

    for (int y=0; y<height; y+=2) {
      for (int x=0; x<width; x+=2) {
        double best_fbs[2][2][2][3];
        double best_score[2][2] = {{-1, -1}, {-1, -1}};
        
        double* fbt[2];
        for (int bxdiff=-3; bxdiff<=3; bxdiff++) {
          for (int bydiff=-3; bydiff<=3; bydiff++) {
            if (!(x/2 + bxdiff >= 0 && y/2 + bydiff >= 0
                  && x/2 + bxdiff < width/2 && y + bydiff < height/2))
              continue;

            fbt[0] = &(final_list[raise-1][1][((y/2+bydiff)*width/2+(x/2+bxdiff))*3]);
            for (int fxdiff=-3; fxdiff<=3; fxdiff++) {
              for (int fydiff=-3; fydiff<=3; fydiff++) {
                if (!(x/2 + fxdiff >= 0 && y + fydiff >= 0
                      && x/2 + fxdiff < width/2 && y + fydiff < height/2))
                  continue;

                fbt[1] = &(final_list[raise-1][0]
                    [((y/2+fydiff)*width/2+(x/2+fxdiff))*3]);

                for (int nx=0; nx<2; nx++) {
                  for (int ny=0; ny<2; ny++) {
                    double fbs_merged[2][3];
                    double proj_a;

                    for (int c=0; c<3; c++) {
                      for (int b=0; b<2; b++){
                        // TODO fix verhältnis?
                        double cur_ps = portion_list[raise][b][(y*ny)*width+x+nx];
                        fbs_merged[b][c] = color_list[raise][b][((y+ny)*width+x+nx)*3+c]
                          * cur_ps + (1 - cur_ps)*fbt[b][c];
                      }
                    }

                    double score = projection(fbs_merged[0], fbs_merged[1],
                        &(original_list[raise][((y+ny)*width+x+nx)*3]), &proj_a);

                    if (score < best_score[nx][ny]
                        || best_score[nx][ny] == -1) {
                      best_score[nx][ny] = score;
                      for (int c=0; c<3; c++) {
                        for (int b=0; b<2; b++) {
                          best_fbs[nx][ny][b][c] = fbt[b][c];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }

        for (int nx=0; nx<2; nx++) {
          for (int ny=0; ny<2; ny++) {
            for (int b=0; b<2; b++) {
              for (int c=0; c<3; c++) {
                double cur_ps = portion_list[raise][b][(y+ny)*width+x+nx];
                final_list[raise][b][((y+ny)*width+x+nx)*3+c] =
                  ((best_fbs[nx][ny][b][c]
                  + color_list[raise-1][b][(y/2*width/2+x/2)*3+c])/2)*(1-cur_ps)
                  + color_list[raise][b][((y+ny)*width+x+nx)]*cur_ps;
              }
              projection(&(final_list[raise][0][((y+ny)*width+x+nx)*3]),
                  &(final_list[raise][1][((y+ny)*width+x+nx)*3]),
                  &(original_list[raise][((y*ny)*width+x+nx)*3]),
                  &(alpha_list[raise][((y*ny)*width+x+nx)*3]));
            }
          }
        }
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "final_%i.ppm", raise);
    double *tmp = new double[width*height*3];
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int c=0; c<3; c++) {
          tmp[(y*width+x)*3+c] = portion_list[raise][0][y*width+1]
            *color_list[raise][0][(y*width+1)*3+c]
            +final_list[raise][0][(y*width+x)*3+c]*alpha_list[raise][(y*width+x)]
            *(1-portion_list[raise][0][y*width+x]-portion_list[raise][1][y*width+x]);
        }
      }
    }
    save_image(buffer, width, height, 3, tmp);
    delete[] tmp;

    snprintf(buffer, 100, RESULTS "new_foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, final_list[raise][0]);
    
    snprintf(buffer, 100, RESULTS "new_backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, final_list[raise][1]);
    
    snprintf(buffer, 100, RESULTS "new_alphas_%i.ppm", raise);
    save_image(buffer, width, height, 1, alpha_list[raise]);
    
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
