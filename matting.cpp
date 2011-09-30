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

// a1 and f1 known
double solve_equations(unsigned char* f0, unsigned char* f1, unsigned char* f2,
    unsigned char* b0, unsigned char* b1, unsigned char* b2,
    unsigned char* a0, unsigned char* a1, unsigned char* a2,
    unsigned char* c1, unsigned char* c2,
    double up0, double up1, double up2) {
#ifdef DEBUG_SOLVE_EQUATIONS
  cout << "\nSolve\n";
  cout << "a0: " << (int)a0[0] << "\n";
  cout << "a1: " << (int)a1[0] << "\n";
  cout << "up0: " << up0 << "\n";
  cout << "up1: " << up1 << "\n";
  cout << "up2: " << up2 << "\n";
  cout << "f0\tf1\tb0\tc1\tc2\n";
  for (int c=0; c<3; c++) {
    cout << (int)f0[c] << "\t" << (int)f1[c] << "\t" << (int)b0[c] << "\t" <<
      (int)c1[c] << "\t" << (int)c2[c] << "\n";
  }
#endif

  int nb1[3];
  int nf2[3];
  int nb2[3];
  a2[0] = (2*a0[0]*up0-a1[0]*up1)/up2;
  for (int c=0; c<3; c++) {
    // set the same fore/background as in the combined pixel if it is not used
    float alpha1 = a1[0]/255.;
    float alpha2 = a2[0]/255.;
    float alpha0 = a0[0]/255.;
    float background_alpha1 = 1-alpha1;
    float background_alpha2 = 1-alpha2;

    if (background_alpha1 != 0.0) {
      nb1[c] = (c1[c] - f1[c]*alpha1)/background_alpha1;
    } else {
      nb1[c] = b0[c];
    }

    if (alpha2 != 0.0) {
      nf2[c] = (2*f0[c]*up0*alpha0 - f1[c]*alpha1*up1)/(alpha2*up2);
    } else {
      nf2[c] = f0[c];
    }

    if (background_alpha2 != 0.0) {
      nb2[c] = (c2[c] - nf2[c]*alpha2)/background_alpha2;
    } else {
      nb2[c] = b0[c];
    }
  }

#ifdef DEBUG_SOLVE_EQUATIONS
  cout << "new_alpha2: " << (int)a2[0] << "\n";
  cout << "new_b1\tnew_b2\tnew_f2\n";
  for (int c=0; c<3; c++) {
    cout << (int)nb1[c] << "\t" << (int)nb2[c] << "\t" <<
      (int)nf2[c] << "\n";
  }
#endif

  double quality = 0;
  for (int c=0; c<3; c++) {
    double diff = nb1[c]-255;
    if (diff > 0) {
      quality -= diff*diff;
      nb1[c] = 255;
    } else if (nb1[c] < 0) {
      quality -= nb1[c]*nb1[c];
      nb1[c] = 0;
    }

    diff = nf2[c]-255;
    if (diff > 0) {
      quality -= diff*diff;
      nf2[c] = 255;
    } else if (nf2[c] < 0) {
      quality -= nf2[c]*nf2[c];
      nf2[c] = 0;
    }
    
    diff = nb2[c]-255;
    if (diff > 0) {
      quality -= diff*diff;
      nb2[c] = 255;
    } else if (nb2[c] < 0) {
      quality -= nb2[c]*nb2[c];
      nb2[c] = 0;
    }

    b1[c] = nb1[c];
    f2[c] = nf2[c];
    b2[c] = nb2[c];
  }

  //quality = 0;
  if (quality == 0) {
    for (int c=0; c<3; c++) {
      int diff = f1[c]-f0[c];
      quality += diff*diff;
      
      diff = f2[c]-f0[c];
      quality += diff*diff;
      
      diff = b1[c]-b0[c];
      quality += diff*diff;
      
      diff = b2[c]-b0[c];
      quality += diff*diff;
    }
      
    int diff = a1[0]-a0[0];
    quality += diff*diff;

    diff = a2[0]-a0[0];
    quality += diff*diff;
  }
  return quality;
}

// a1 and f1 unknown
void optimize(unsigned char* f0, unsigned char* f1, unsigned char* f2,
    unsigned char* b0, unsigned char* b1, unsigned char* b2,
    unsigned char* a0, unsigned char* a1, unsigned char* a2,
    unsigned char* c1, unsigned char* c2,
    double up0, double up1, double up2) {

  unsigned char best_f1[3] = {0};
  unsigned char best_a = 0;

  int best_result = INT_MIN;
  for (a1[0] = 0; a1[0] < 251; a1[0]+=4) {
    cout << (int)a1[0] << "/255" << endl;
    for (f1[0] = 0; f1[0] < 251; f1[0]+=4) {
      for (f1[1] = 0; f1[1] < 251; f1[1]+=4) {
        for (f1[2] = 0; f1[2] < 251; f1[2]+=4) {
          int result = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2, up0, up1, up2);
          if (best_result < 0) {
            if (result > best_result) {
              best_result = result;
              best_a = a1[0];
              for (int c=0; c<3; c++) {
                best_f1[c] = f1[c];
              }
#ifdef DEBUG_OPTIMIZE_LOOP
              cout << "New best result: " << best_result << "\n";
              cout << "a1: " << (int)a1[0] << "\n";
              cout << "a2: " << (int)a2[0] << "\n";
              cout << "f1\tf2\tb1\tb2\n";
              for (int c=0; c<3; c++) {
                cout << (int)f1[c] << "\t" << (int)f2[c] << "\t" << (int)b1[c] << "\t" <<
                  (int)b2[c] << "\n";
              }
#endif
            }
          } else {
            if (result > 0 && result < best_result) {
              best_result = result;
              best_a = a1[0];
              for (int c=0; c<3; c++) {
                best_f1[c] = f1[c];
              }
#ifdef DEBUG_OPTIMIZE_LOOP
              cout << "New best result: " << best_result << "\n";
              cout << "a1: " << (int)a1[0] << "\n";
              cout << "a2: " << (int)a2[0] << "\n";
              cout << "f1\tf2\tb1\tb2\n";
              for (int c=0; c<3; c++) {
                cout << (int)f1[c] << "\t" << (int)f2[c] << "\t" << (int)b1[c] << "\t" <<
                  (int)b2[c] << "\n";
              }
#endif
            }
          }
        }
      }
    }
  }

  /*f1[0] = f0[0];
  f1[1] = f0[1];
  f1[2] = f0[2];
  a1[0] = 128;*/

      /*
    f1[0] = rand()%256;
    f1[1] = rand()%256;
    f1[2] = rand()%256;
    a1[0] = rand()%256;

    for (int i=0; i<100; i++) {
      // alpha
      int qplus, qminus;
      if (a1[0] != 255)
        a1[0]++;
      qplus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
      if (a1[0] != 255)
        a1[0]--;

      if (a1[0] != 0)
        a1[0]--;
      qminus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
      if (a1[0] != 0)
        a1[0]++;

      cout << (int)f1[0] << "/" << (int)f1[1] << "/" << (int)f1[2] << "/" << (int)a1[0]+1 << ": " << qplus << endl;

      if (qplus > 0)
        if (qminus > 0)
          if (qplus < qminus)
            a1[0]++;
          else
            a1[0]--;
        else
          a1[0]++;
      else
        if (qminus > 0)
          a1[0]--;
        else
          if (qminus > qplus)
            a1[0]--;
          else
            a1[0]++;

      for (int c=0; c<3; c++) {
        int qplus, qminus;
        if (f1[c] != 255) {
          f1[c]++;
          qplus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
          f1[c]--;
        } else {
          qplus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
        }

        if (f1[c] != 0) {
          f1[c]--;
          qminus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
          f1[c]++;
        } else {
          qminus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
        }

        if (qplus > 0) {
          if (qminus > 0) {
            if (qplus < qminus) {
              if (f1[c] != 255)
                f1[c]++;
            } else {
              if (f1[c] != 0)
                f1[c]--;
            }
          } else {
            if (f1[c] != 255)
              f1[c]++;
          }
        } else {
          if (qminus > 0) {
            if (f1[c] != 0)
              f1[c]--;
          } else {
            if (qminus > qplus) {
              if (f1[c] != 0)
                f1[c]--;
            } else {
              if (f1[c] != 255)
                f1[c]++;
            }
          }
        }
      }
    }
  }*/
  
  a1[0] = best_a;
  for (int c=0; c<3; c++) {
    f1[c] = best_f1[c];
  }
  solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2, up0, up1, up2);

  //cout << "alpha 0/1/2: " << (int)a0[0] << "/" << (int)a1[0] << "/" << (int)a2[0] << "\n";
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
  double* originals[raise+1];

  double* new_alphas[raise+1];
  double* fbs[raise+1][2];
  double* ps[raise+1][2];
  double* new_colors[raise+1][2];

  originals[9] = load_image("test.ppm", width, height, 3);
  for (int b=0; b<2; b++) {
    fbs[9][b] = new double[width*height*3]();
    ps[9][b] = new double[width*height]();
  }

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      if (mask[width*y+x] == 1) {
        ps[9][0][width*y+x] = 1;
        for (int c=0; c<3; c++) {
          fbs[9][0][(width*y+x)*3+c] =
            originals[9][(width*y+x)*3+c];
        }
      } else if (mask[width*y+x] == 0) {
        ps[9][1][width*y+x] = 1;
        for (int c=0; c<3; c++) {
          fbs[9][1][(width*y+x)*3+c] =
            originals[9][(width*y+x)*3+c];
        }
      }
    }
  }

  save_image(RESULTS "originals_9.ppm", width, height, 3, originals[9]);
  save_image(RESULTS "foregrounds_9.ppm", width, height, 3, fbs[9][0]);
  save_image(RESULTS "backgrounds_9.ppm", width, height, 3, fbs[9][1]);

  while (raise > 0) {
    raise--;

    // Height half
    height = height/2;
    width = width/2;

    originals[raise] = new double[width*height*3]();
    for (int b=0; b<2; b++) {
      fbs[raise][b] = new double[width*height*3]();
      ps[raise][b] =  new double[width*height]();
    }
  
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int xdiff = 0; xdiff <= 1; xdiff++) {
          for (int ydiff = 0; ydiff <= 1; ydiff++) {
            for (int b=0; b<2; b++) {
              ps[raise][b][y*width+x] +=
                ps[raise+1][b][(2*y+ydiff)*(2*width)+2*x+xdiff]/4;
            }
          }
        }

        for (int xdiff = 0; xdiff <= 1; xdiff++) {
          for (int ydiff = 0; ydiff <= 1; ydiff++) {
            for (int c=0; c<3; c++) {
              originals[raise][(y*width+x)*3+c] +=
                originals[raise+1][((2*y+ydiff)*(2*width)+2*x+xdiff)*3+c]/4;

              for (int b=0; b<2; b++) {
                double new_ps = ps[raise][b][y*width+x];
                if (new_ps > 0)
                  fbs[raise][b][(y*width+x)*3+c] +=
                    fbs[raise+1][b][((2*y+ydiff)*(2*width)+2*x+xdiff)*3+c]
                    *ps[raise+1][b][(2*y+ydiff)*(2*width)+2*x+xdiff]
                    /new_ps/4;
              }
            }
          }
        }
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "originals_%i.ppm", raise);
    save_image(buffer, width, height, 3, originals[raise]);
    
    snprintf(buffer, 100, RESULTS "foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, fbs[raise][0]);
    
    snprintf(buffer, 100, RESULTS "backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, fbs[raise][1]);
  }
  

  // Upscaling
  raise = 0;
  width = 1;
  height = 1;

  new_alphas[0] = new double[1];
  for (int b=0; b<2; b++) {
    new_colors[0][b] = new double[3];
    for (int c=0; c<3; c++) {
      new_colors[0][b][c] = fbs[0][b][c];
    }
  }

  // Calculate alpha
  double merged_point[3];
  for (int c=0; c<3; c++) {
    merged_point[c] = (originals[0][c]
      -ps[0][0][0]*fbs[0][0][c]
      -ps[0][1][0]*fbs[0][1][c])
      /(1-ps[0][0][0]-ps[0][1][0]);
  }

  projection(new_colors[0][0], new_colors[0][1],
      merged_point, &(new_alphas[0][0]));

  save_image(RESULTS "new_foregrounds_0.ppm", 1, 1, 3, new_colors[0][0]);
  save_image(RESULTS "new_backgrounds_0.ppm", 1, 1, 3, new_colors[0][1]);
  save_image(RESULTS "new_alphas_0.ppm", 1, 1, 1, new_alphas[0]);


  double final[3];
  for (int c=0; c<3; c++) {
    final[c] = ps[0][0][0]*fbs[0][0][c]
      +new_colors[0][0][c]*new_alphas[0][0]
      *(1-ps[0][0][0]-ps[0][1][0]);
  }
  save_image(RESULTS "final_0.ppm", 1, 1, 3, final);

  while (raise < 9) {
    cout << "============= raise: " << raise << "=============" << endl;
    width *= 2;
    height *= 2;
    raise++;

    for (int b=0; b<2; b++) {
      new_colors[raise][b] = new double[width*height*3];
    }
    new_alphas[raise] = new double[width*height];

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

            fbt[0] = &(new_colors[raise-1][1][((y/2+bydiff)*width/2+(x/2+bxdiff))*3]);
            for (int fxdiff=-3; fxdiff<=3; fxdiff++) {
              for (int fydiff=-3; fydiff<=3; fydiff++) {
                if (!(x/2 + fxdiff >= 0 && y + fydiff >= 0
                      && x/2 + fxdiff < width/2 && y + fydiff < height/2))
                  continue;

                fbt[1] = &(new_colors[raise-1][0]
                    [((y/2+fydiff)*width/2+(x/2+fxdiff))*3]);

                for (int nx=0; nx<2; nx++) {
                  for (int ny=0; ny<2; ny++) {
                    double fbs_merged[2][3];
                    double proj_a;

                    for (int c=0; c<3; c++) {
                      for (int b=0; b<2; b++){
                        // TODO fix verhältnis?
                        double cur_ps = ps[raise][b][(y*ny)*width+x+nx];
                        fbs_merged[b][c] = fbs[raise][b][((y+ny)*width+x+nx)*3+c]
                          * cur_ps + (1 - cur_ps)*fbt[b][c];
                      }
                    }

                    double score = projection(fbs_merged[0], fbs_merged[1],
                        &(originals[raise][((y+ny)*width+x+nx)*3]), &proj_a);

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
                double cur_ps = ps[raise][b][(y+ny)*width+x+nx];
                new_colors[raise][b][((y+ny)*width+x+nx)*3+c] =
                  ((best_fbs[nx][ny][b][c]
                  + fbs[raise-1][b][(y/2*width/2+x/2)*3+c])/2)*(1-cur_ps)
                  + fbs[raise][b][((y+ny)*width+x+nx)]*cur_ps;
              }
              projection(&(new_colors[raise][0][((y+ny)*width+x+nx)*3]),
                  &(new_colors[raise][1][((y+ny)*width+x+nx)*3]),
                  &(originals[raise][((y*ny)*width+x+nx)*3]),
                  &(new_alphas[raise][((y*ny)*width+x+nx)*3]));
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
          tmp[(y*width+x)*3+c] = ps[raise][0][y*width+1]
            *fbs[raise][0][(y*width+1)*3+c]
            +new_colors[raise][0][(y*width+x)*3+c]*new_alphas[raise][(y*width+x)]
            *(1-ps[raise][0][y*width+x]-ps[raise][1][y*width+x]);
        }
      }
    }
    save_image(buffer, width, height, 3, tmp);
    delete[] tmp;

    snprintf(buffer, 100, RESULTS "new_foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, new_colors[raise][0]);
    
    snprintf(buffer, 100, RESULTS "new_backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, new_colors[raise][1]);
    
    snprintf(buffer, 100, RESULTS "new_alphas_%i.ppm", raise);
    save_image(buffer, width, height, 1, new_alphas[raise]);
    
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
