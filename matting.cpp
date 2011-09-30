#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <limits.h>

using namespace std;

//#define DEBUG_PROJECTION
//#define DEBUG_GET_VALUES
//#define DEBUG_SOLVE_EQUATIONS
//#define DEBUG_OPTIMIZE_LOOP

unsigned char*
load_image (const char* filename, int dimx, int dimy, int num_colors) {
  unsigned char* data = new unsigned char[dimx*dimy*num_colors];

  FILE *fp = fopen (filename, "rb");

  fread(data, 1, dimx*dimy*num_colors, fp);
  
  fclose (fp);

  return data;
}


void
save_image (const char* filename, int dimx, int dimy, int num_colors, unsigned char* data) {
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

  fwrite(data, 1, dimx*dimy*num_colors, fp);
  
  fclose (fp);
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

double projection (unsigned char F[3], unsigned char B[3], unsigned char C[3], unsigned char* alpha_pointer) {
#ifdef DEBUG_PROJECTION
  cout << "F\tB\tC\n";
  for (int c=0; c<3; c++) {
    cout << (int)F[c] << "\t" << (int)B[c] << "\t" << (int)C[c] << "\n";
  }
#endif

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

#ifdef DEBUG_PROJECTION
  cout << "=> alpha = " << alpha << endl;
#endif
  
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

/*
  double new_B[3];
  double new_F[3];

  for (int c=0; c<3; c++) {
    new_B[c] = B[c];
    new_F[c] = F[c];
    new_B[c] = B[c] + CC[c];
    new_F[c] = F[c] + CC[c];
  }

#ifdef DEBUG_PROJECTION
  cout << "\nF\t=>\tnew_F\n";
  for (int c=0; c<3; c++) {
    cout << (int)F[c] << "\t" << "" << "\t" << new_F[c] << "\n";
  }
  cout << "\nB\t=>\tnew_B\n";
  for (int c=0; c<3; c++) {
    cout << (int)B[c] << "\t" << "" << "\t" << new_B[c] << "\n";
  }
#endif*/
/*
  for (int c=0; c<3; c++) {
    if (new_B[c] < 0) {
      double factor = (0-C[c])/(new_B[c]-C[c]);
      for (int new_c=0; new_c<3; new_c++) {
        new_B[new_c] = C[new_c]+factor*(new_B[new_c]-C[new_c]);
      }
    }
    if (new_B[c] > 255) {
      double factor = (255-C[c])/(new_B[c]-C[c]);
      for (int new_c=0; new_c<3; new_c++) {
        new_B[new_c] = C[new_c]+factor*(new_B[new_c]-C[new_c]);
      }
    }
    if (new_F[c] < 0) {
      double factor = (0-C[c])/(new_F[c]-C[c]);
      for (int new_c=0; new_c<3; new_c++) {
        new_F[new_c] = C[new_c]+factor*(new_F[new_c]-C[new_c]);
      }
    }
    if (new_F[c] > 255) {
      double factor = (255-C[c])/(new_F[c]-C[c]);
      for (int new_c=0; new_c<3; new_c++) {
        new_F[new_c] = C[new_c]+factor*(new_F[new_c]-C[new_c]);
      }
    }
  }
*/
#ifdef DEBUG_PROJECTION
  cout << "\nNormalize:\n";
  cout << "new_F\t\tnew_B\n";
  for (int c=0; c<3; c++) {
    cout << new_F[c] << "\t\t" << new_B[c] << "\n";
  }
#endif
  
//  double new_alpha = (C[0] - new_B[0])/(new_F[0] - new_B[0]);
#ifdef DEBUG_PROJECTION
  cout << "=> new_alpha = " << new_alpha << endl;
#endif

  *alpha_pointer = (unsigned char)(alpha*255);
/*
  for (int c=0; c<3; c++) {
    B[c] = new_B[c];
    F[c] = new_F[c];
  }*/
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

  unsigned char* mask = load_image("trimap.pnm", width, height, 1);

  unsigned char* originals[raise+1][2];
  unsigned char* new_alphas[raise+1][2];
  unsigned char* foregrounds[raise+1][2];
  unsigned char* new_foregrounds[raise+1][2];
  double* foreground_ps[raise+1][2];
  unsigned char* backgrounds[raise+1][2];
  unsigned char* new_backgrounds[raise+1][2];
  double* background_ps[raise+1][2];

  originals[9][0] = load_image("test.ppm", width, height, 3);
  foregrounds[9][0] = new unsigned char[width*height*3];
  foreground_ps[9][0] = new double[width*height];
  backgrounds[9][0] = new unsigned char[width*height*3];
  background_ps[9][0] = new double[width*height];

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      if (mask[width*y+x] == 255) {
        foreground_ps[9][0][width*y+x] = 1;
        background_ps[9][0][width*y+x] = 0;
        for (int c=0; c<3; c++) {
          backgrounds[9][0][(width*y+x)*3+c] = 0;
          foregrounds[9][0][(width*y+x)*3+c] =
            originals[9][0][(width*y+x)*3+c];
        }
      } else if (mask[width*y+x] == 0) {
        foreground_ps[9][0][width*y+x] = 0;
        background_ps[9][0][width*y+x] = 1;
        for (int c=0; c<3; c++) {
          foregrounds[9][0][(width*y+x)*3+c] = 0;
          backgrounds[9][0][(width*y+x)*3+c] =
            originals[9][0][(width*y+x)*3+c];
        }
      } else {
        foreground_ps[9][0][width*y+x] = 0;
        background_ps[9][0][width*y+x] = 0;
        for (int c=0; c<3; c++) {
          foregrounds[9][0][(width*y+x)*3+c] = 0;
          backgrounds[9][0][(width*y+x)*3+c] = 0;
        }
      }
    }
  }

  save_image(RESULTS "originals_9.ppm", width, height, 3, originals[9][0]);
  save_image(RESULTS "foregrounds_9.ppm", width, height, 3, foregrounds[9][0]);
  save_image(RESULTS "backgrounds_9.ppm", width, height, 3, backgrounds[9][0]);

  while (raise > 0) {
    raise--;

    // Height half
    height = height/2;
    originals[raise][1] = new unsigned char[width*height*3];
    foregrounds[raise][1] = new unsigned char[width*height*3];
    foreground_ps[raise][1] = new double[width*height];
    backgrounds[raise][1] = new unsigned char[width*height*3];
    background_ps[raise][1] = new double[width*height];
  
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        foreground_ps[raise][1][y*width+x] =
          (foreground_ps[raise+1][0][2*y*width+x] +
           foreground_ps[raise+1][0][(2*y+1)*width+x])/2;
        background_ps[raise][1][y*width+x] =
          (background_ps[raise+1][0][2*y*width+x] +
           background_ps[raise+1][0][(2*y+1)*width+x])/2;

        for (int c=0; c<3; c++) {
          originals[raise][1][(y*width+x)*3+c] =
            (originals[raise+1][0][(2*y*width+x)*3+c] +
             originals[raise+1][0][((2*y+1)*width+x)*3+c])/2;

          foregrounds[raise][1][(y*width+x)*3+c] =
            (foregrounds[raise+1][0][(2*y*width+x)*3+c]*
             foreground_ps[raise+1][0][2*y*width+x] +
             foregrounds[raise+1][0][((2*y+1)*width+x)*3+c]*
             foreground_ps[raise+1][0][(2*y+1)*width+x])/
             foreground_ps[raise][1][y*width+x]/2;
          
          backgrounds[raise][1][(y*width+x)*3+c] =
            (backgrounds[raise+1][0][(2*y*width+x)*3+c]*
             background_ps[raise+1][0][2*y*width+x] +
             backgrounds[raise+1][0][((2*y+1)*width+x)*3+c]*
             background_ps[raise+1][0][(2*y+1)*width+x])/
             background_ps[raise][1][y*width+x]/2;
        }
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "originals_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, originals[raise][1]);
    
    snprintf(buffer, 100, RESULTS "foregrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, foregrounds[raise][1]);
    
    snprintf(buffer, 100, RESULTS "backgrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, backgrounds[raise][1]);

    // Width half
    width = width/2;
    originals[raise][0] = new unsigned char[width*height*3];
    foregrounds[raise][0] = new unsigned char[width*height*3];
    foreground_ps[raise][0] = new double[width*height];
    
    backgrounds[raise][0] = new unsigned char[width*height*3];
    background_ps[raise][0] = new double[width*height];

    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        foreground_ps[raise][0][y*width+x] =
          (foreground_ps[raise][1][y*width*2+x*2] +
           foreground_ps[raise][1][y*width*2+x*2+1])/2;

        background_ps[raise][0][y*width+x] =
          (background_ps[raise][1][y*width*2+x*2] +
           background_ps[raise][1][y*width*2+x*2+1])/2;

        for (int c=0; c<3; c++) {
          originals[raise][0][(y*width+x)*3+c] =
            (originals[raise][1][(y*width*2+2*x)*3+c] +
             originals[raise][1][(y*width*2+2*x+1)*3+c])/2;

          foregrounds[raise][0][(y*width+x)*3+c] =
            (foregrounds[raise][1][(y*width*2+x*2)*3+c]*
             foreground_ps[raise][1][y*width*2+x*2] +
             foregrounds[raise][1][(y*width*2+x*2+1)*3+c]*
             foreground_ps[raise][1][y*width*2+x*2+1])/
             foreground_ps[raise][0][y*width+x]/2;
          
          backgrounds[raise][0][(y*width+x)*3+c] =
            (backgrounds[raise][1][(y*width*2+x*2)*3+c]*
             background_ps[raise][1][y*width*2+x*2] +
             backgrounds[raise][1][(y*width*2+x*2+1)*3+c]*
             background_ps[raise][1][y*width*2+x*2+1])/
             background_ps[raise][0][y*width+x]/2;
        }
      }
    }

    snprintf(buffer, 100, RESULTS "originals_%i.ppm", raise);
    save_image(buffer, width, height, 3, originals[raise][0]);
    
    snprintf(buffer, 100, RESULTS "foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, foregrounds[raise][0]);
    
    snprintf(buffer, 100, RESULTS "backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, backgrounds[raise][0]);
  }

  // Upscaling
  raise = 0;
  width = 1;
  height = 1;

  new_foregrounds[0][0] = new unsigned char[3];
  new_backgrounds[0][0] = new unsigned char[3];
  new_alphas[0][0] = new unsigned char[1];

  for (int c=0; c<3; c++) {
    new_foregrounds[0][0][c] = foregrounds[0][0][c];
    new_backgrounds[0][0][c] = backgrounds[0][0][c];
  }

  // Calculate alpha
  unsigned char merged_point[3];
  for (int c=0; c<3; c++) {
    merged_point[c] = (originals[0][0][c]
      -foreground_ps[0][0][0]*foregrounds[0][0][c]
      -background_ps[0][0][0]*backgrounds[0][0][c])
      /(1-foreground_ps[0][0][0]-background_ps[0][0][0]);
  }

  projection(new_foregrounds[0][0], new_backgrounds[0][0], merged_point, &(new_alphas[0][0][0]));

  save_image(RESULTS "new_foregrounds_0.ppm", 1, 1, 3, new_foregrounds[0][0]);
  save_image(RESULTS "new_backgrounds_0.ppm", 1, 1, 3, new_backgrounds[0][0]);
  save_image(RESULTS "new_alphas_0.ppm", 1, 1, 1, new_alphas[0][0]);
  unsigned char final[3];
  for (int c=0; c<3; c++) {
    final[c] = foreground_ps[0][0][0]*foregrounds[0][0][c]
      +new_foregrounds[0][0][c]*new_alphas[0][0][0]/255.
      *(1-foreground_ps[0][0][0]-background_ps[0][0][0]);
  }
  save_image(RESULTS "final_0.ppm", 1, 1, 3, final);

  while (raise < 9) {
    cout << "============= size: " << pow(2, raise) << "=============" << endl;
    width *= 2;
    new_foregrounds[raise][1] = new unsigned char[width*height*3];
    new_backgrounds[raise][1] = new unsigned char[width*height*3];
    new_alphas[raise][1] = new unsigned char[width*height];

    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x+=2) {
        unsigned char* f0 = &(new_foregrounds[raise][0][(y*width/2+x/2)*3]);
        unsigned char* b0 = &(new_backgrounds[raise][0][(y*width/2+x/2)*3]);
        unsigned char* a0 = &(new_alphas[raise][0][(y*width/2+x/2)]);

        unsigned char* f1 = &(new_foregrounds[raise][1][(y*width+x)*3]);
        unsigned char* f2 = &(new_foregrounds[raise][1][(y*width+x+1)*3]);
        
        unsigned char* b1 = &(new_backgrounds[raise][1][(y*width+x)*3]);
        unsigned char* b2 = &(new_backgrounds[raise][1][(y*width+x+1)*3]);
        
        unsigned char* a1 = &(new_alphas[raise][1][y*width+x]);
        unsigned char* a2 = &(new_alphas[raise][1][y*width+x+1]);

        // Subtract known foreground / background 
        unsigned char* original1 = &(originals[raise][1][(y*width+x)*3]);
        unsigned char* original2 = &(originals[raise][1][(y*width+x+1)*3]);

        unsigned char* foreground1 = &(foregrounds[raise][1][(y*width+x)*3]);
        unsigned char* foreground2 = &(foregrounds[raise][1][(y*width+x+1)*3]);
        
        unsigned char* background1 = &(backgrounds[raise][1][(y*width+x)*3]);
        unsigned char* background2 = &(backgrounds[raise][1][(y*width+x+1)*3]);

        double* foreground_ps1 = &(foreground_ps[raise][1][(y*width+x)]);
        double* foreground_ps2 = &(foreground_ps[raise][1][(y*width+x+1)]);
        
        double* background_ps1 = &(background_ps[raise][1][(y*width+x)]);
        double* background_ps2 = &(background_ps[raise][1][(y*width+x+1)]);

#ifdef DEBUG_GET_VALUES
        cout << "a0:\t" << (int)a0[0] << endl;
        cout << "fg_ps1:\t" << foreground_ps1[0]
          << "\nfg_ps2:\t" << foreground_ps2[0]
          << "\nbg_ps1:\t" << background_ps1[0]
          << "\nbg_ps2:\t" << background_ps2[0] << "\n";
        cout << "f0\tb0\torig1\torig2\tbg1\tbg2\tfg1\tfg2\t\n";
        for (int c=0; c<3; c++) {
          cout << (int)f0[c] << "\t" << (int)b0[c] << "\t"
            << (int)original1[c] << "\t" << (int)original2[c] << "\t"
            << (int)background1[c] << "\t" << (int)background2[c] << "\t"
            << (int)foreground1[c] << "\t" << (int)foreground2[c] << "\n";
        }
#endif

        unsigned char c1[3];
        unsigned char c2[3];

        for (int c=0; c<3; c++) {
          c1[c] = (original1[c] - foreground_ps1[0]*foreground1[c]
              -background_ps1[0]*background1[c])
            /(1-foreground_ps1[0]-background_ps1[0]);
          c2[c] = (original2[c] - foreground_ps2[0]*foreground2[c]
              -background_ps2[0]*background2[c])
            /(1-foreground_ps2[0]-background_ps2[0]);
        }

#ifdef DEBUG_GET_VALUES
        cout << "\nc1\tc2\n";
        for (int c=0; c<3; c++) {
          cout << (int)c1[c] << "\t" << (int)c2[c] << "\n";
        }
#endif
        unsigned char best_fg1[3];
        unsigned char best_bg1[3];
        unsigned char best_fg2[3];
        unsigned char best_bg2[3];
        unsigned char best_a1;
        unsigned char best_a2;
        double best_1 = -1;
        double best_2 = -1;
        int count;
        for (int bxdiff=-1; bxdiff<=1; bxdiff++) {
          for (int bydiff=-1; bydiff<=1; bydiff++) {
            if (!(x/2+bxdiff>=0 && y+bydiff>=0 && x/2+bxdiff < width/2 && y+bydiff < height))
              continue;

            unsigned char* bgt = &(new_backgrounds[raise][0][((y+bydiff)*width/2+(x/2+bxdiff))*3]);
            for (int fxdiff=-3; fxdiff<=3; fxdiff++) {
              for (int fydiff=-3; fydiff<=3; fydiff++) {
                if (!(x/2+fxdiff>=0 && y+fydiff>=0 && x/2+fxdiff < width/2 && y+fydiff < height))
                  continue;

                unsigned char* fgt = &(new_foregrounds[raise][0][((y+fydiff)*width/2+(x/2+fxdiff))*3]);
                for (int c=0; c<3; c++) {
                  f1[c]=foreground1[c]*foreground_ps1[0]+
                    (1-foreground_ps1[0])*(fgt[c]);
                  b1[c]=background1[c]*background_ps1[0]+
                    (1-background_ps1[0])*(bgt[c]);
                  f2[c]=foreground2[c]*foreground_ps2[0]
                    +(1-foreground_ps2[0])*(fgt[c]);
                  b2[c]=background2[c]*background_ps2[0]
                    +(1-background_ps2[0])*(bgt[c]);

                }
                double score1 = projection(f1, b1, original1, a1);
                if (score1 < best_1 || best_1 == -1) {
                  best_1 = score1;
                  best_a1 = a1[0];
                  for (int c=0; c<3; c++) {
                    best_fg1[c] = f1[c];
                    best_bg1[c] = b1[c];
                  }
                }
                double score2 = projection(f2, b2, original2, a2);
                if (score2 < best_2 || best_2 == -1) {
                  best_2 = score2;
                  best_a2 = a2[0];
                  for (int c=0; c<3; c++) {
                    best_fg2[c] = f2[c];
                    best_bg2[c] = b2[c];
                  }
                }
              }
            }
          }
        }
        for (int c=0; c<3; c++) {
          f1[c] = (best_fg1[c]+f0[c])/2;
          f2[c] = (best_fg2[c]+f0[c])/2;
          b1[c] = (best_bg1[c]+b0[c])/2;
          b2[c] = (best_bg2[c]+b0[c])/2;
        }
        for (int c=0; c<3; c++) {
          f1[c]=foreground1[c]*foreground_ps1[0]+
            (1-foreground_ps1[0])*(f1[c]);
          b1[c]=background1[c]*background_ps1[0]+
            (1-background_ps1[0])*(b1[c]);
          f2[c]=foreground2[c]*foreground_ps2[0]
            +(1-foreground_ps2[0])*(f2[c]);
          b2[c]=background2[c]*background_ps2[0]
            +(1-background_ps2[0])*(b2[c]);

        }
        projection(f1, b1, original1, a1);
        projection(f2, b2, original2, a2);
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "final_%i_h.ppm", raise);
    unsigned char *tmp = new unsigned char[width*height*3];
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int c=0; c<3; c++) {
          tmp[(y*width+x)*3+c] =
            foregrounds[raise][1][(y*width+x)*3+c]*
            foreground_ps[raise][1][(y*width+x)]+

            (new_foregrounds[raise][1][(y*width+x)*3+c]
            *new_alphas[raise][1][y*width+x]/255.)*

            (1-foreground_ps[raise][1][(y*width+x)]-
            background_ps[raise][1][(y*width+x)]);
        }
      }
    }
    save_image(buffer, width, height, 3, tmp);
    delete[] tmp;

    snprintf(buffer, 100, RESULTS "new_foregrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, new_foregrounds[raise][1]);
    
    snprintf(buffer, 100, RESULTS "new_backgrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, new_backgrounds[raise][1]);
    
    snprintf(buffer, 100, RESULTS "new_alphas_%i_h.ppm", raise);
    save_image(buffer, width, height, 1, new_alphas[raise][1]);
    
    // Stretch height
    cout << "============= size(h): " << pow(2, raise) << "=============" << endl;
    height *= 2;
    new_foregrounds[raise+1][0] = new unsigned char[width*height*3];
    new_backgrounds[raise+1][0] = new unsigned char[width*height*3];
    new_alphas[raise+1][0] = new unsigned char[width*height];

    for (int y=0; y<height; y+=2) {
      for (int x=0; x<width; x++) {
        unsigned char* f0 = &(new_foregrounds[raise][1][((y/2)*width+x)*3]);
        unsigned char* b0 = &(new_backgrounds[raise][1][((y/2)*width+x)*3]);
        unsigned char* a0 = &(new_alphas[raise][1][((y/2)*width+x)]);

        unsigned char* f1 = &(new_foregrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* f2 = &(new_foregrounds[raise+1][0][((y+1)*width+x)*3]);
        
        unsigned char* b1 = &(new_backgrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* b2 = &(new_backgrounds[raise+1][0][((y+1)*width+x)*3]);
        
        unsigned char* a1 = &(new_alphas[raise+1][0][y*width+x]);
        unsigned char* a2 = &(new_alphas[raise+1][0][(y+1)*width+x]);

        // Subtract known foreground / background 
        unsigned char* original1 = &(originals[raise+1][0][(y*width+x)*3]);
        unsigned char* original2 = &(originals[raise+1][0][((y+1)*width+x)*3]);
       
        unsigned char* foreground1 = &(foregrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* foreground2 = &(foregrounds[raise+1][0][((y+1)*width+x)*3]);
        
        unsigned char* background1 = &(backgrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* background2 = &(backgrounds[raise+1][0][((y+1)*width+x)*3]);

        double* foreground_ps1 = &(foreground_ps[raise+1][0][(y*width+x)]);
        double* foreground_ps2 = &(foreground_ps[raise+1][0][((y+1)*width+x)]);
        
        double* background_ps1 = &(background_ps[raise+1][0][(y*width+x)]);
        double* background_ps2 = &(background_ps[raise+1][0][((y+1)*width+x)]);

#ifdef DEBUG_GET_VALUES
        cout << "a0:\t" << (int)a0[0] << endl;
        cout << "fg_ps1:\t" << foreground_ps1[0]
          << "\nfg_ps2:\t" << foreground_ps2[0]
          << "\nbg_ps1:\t" << background_ps1[0]
          << "\nbg_ps2:\t" << background_ps2[0] << "\n";
        cout << "f0\tb0\torig1\torig2\tbg1\tbg2\tfg1\tfg2\t\n";
        for (int c=0; c<3; c++) {
          cout << (int)f0[c] << "\t" << (int)b0[c] << "\t"
            << (int)original1[c] << "\t" << (int)original2[c] << "\t"
            << (int)background1[c] << "\t" << (int)background2[c] << "\t"
            << (int)foreground1[c] << "\t" << (int)foreground2[c] << "\n";
        }
#endif

        unsigned char c1[3];
        unsigned char c2[3];

        for (int c=0; c<3; c++) {
          c1[c] = (original1[c] - foreground_ps1[0]*foreground1[c]
              -background_ps1[0]*background1[c])
            /(1-foreground_ps1[0]-background_ps1[0]);
          c2[c] = (original2[c] - foreground_ps2[0]*foreground2[c]
              -background_ps2[0]*background2[c])
            /(1-foreground_ps2[0]-background_ps2[0]);
        }

#ifdef DEBUG_GET_VALUES
        cout << "\nc1\tc2\n";
        for (int c=0; c<3; c++) {
          cout << (int)c1[c] << "\t" << (int)c2[c] << "\n";
        }
#endif
        unsigned char best_fg1[3];
        unsigned char best_bg1[3];
        unsigned char best_fg2[3];
        unsigned char best_bg2[3];
        unsigned char best_a1;
        unsigned char best_a2;
        double best_1 = -1;
        double best_2 = -1;
        int count;
        for (int bxdiff=-1; bxdiff<=1; bxdiff++) {
          for (int bydiff=-1; bydiff<=1; bydiff++) {
            if (!(x+bxdiff>=0 && y/2+bydiff>=0 && x+bxdiff < width && y/2+bydiff < height/2))
              continue;

            unsigned char* bgt = &(new_backgrounds[raise][1][((y/2+bydiff)*width+(x+bxdiff))*3]);
            for (int fxdiff=-3; fxdiff<=3; fxdiff++) {
              for (int fydiff=-3; fydiff<=3; fydiff++) {
                if (!(x+fxdiff>=0 && y/2+fydiff>=0 && x+fxdiff < width && y/2+fydiff < height/2))
                  continue;

                unsigned char* fgt = &(new_foregrounds[raise][1][((y/2+fydiff)*width+(x+fxdiff))*3]);
                for (int c=0; c<3; c++) {
                  f1[c]=foreground1[c]*foreground_ps1[0]+
                    (1-foreground_ps1[0])*(fgt[c]);
                  b1[c]=background1[c]*background_ps1[0]+
                    (1-background_ps1[0])*(bgt[c]);
                  f2[c]=foreground2[c]*foreground_ps2[0]
                    +(1-foreground_ps2[0])*(fgt[c]);
                  b2[c]=background2[c]*background_ps2[0]
                    +(1-background_ps2[0])*(bgt[c]);

                }
                double score1 = projection(f1, b1, original1, a1);
                if (score1 < best_1 || best_1 == -1) {
                  best_1 = score1;
                  best_a1 = a1[0];
                  for (int c=0; c<3; c++) {
                    best_fg1[c] = f1[c];
                    best_bg1[c] = b1[c];
                  }
                }
                double score2 = projection(f2, b2, original2, a2);
                if (score2 < best_2 || best_2 == -1) {
                  best_2 = score2;
                  best_a2 = a2[0];
                  for (int c=0; c<3; c++) {
                    best_fg2[c] = f2[c];
                    best_bg2[c] = b2[c];
                  }
                }
              }
            }
          }
        }

        for (int c=0; c<3; c++) {
          f1[c] = (best_fg1[c]+f0[c])/2;
          f2[c] = (best_fg2[c]+f0[c])/2;
          b1[c] = (best_bg1[c]+b0[c])/2;
          b2[c] = (best_bg2[c]+b0[c])/2;
        }
        for (int c=0; c<3; c++) {
          f1[c]=foreground1[c]*foreground_ps1[0]+
            (1-foreground_ps1[0])*(f1[c]);
          b1[c]=background1[c]*background_ps1[0]+
            (1-background_ps1[0])*(b1[c]);
          f2[c]=foreground2[c]*foreground_ps2[0]
            +(1-foreground_ps2[0])*(f2[c]);
          b2[c]=background2[c]*background_ps2[0]
            +(1-background_ps2[0])*(b2[c]);

        }
        projection(f1, b1, original1, a1);
        projection(f2, b2, original2, a2);
      }
    }

    snprintf(buffer, 100, RESULTS "final_%i.ppm", raise+1);
    tmp = new unsigned char[width*height*3];
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int c=0; c<3; c++) {
          tmp[(y*width+x)*3+c] =
            foregrounds[raise+1][0][(y*width+x)*3+c]*
            foreground_ps[raise+1][0][(y*width+x)]+

            (new_foregrounds[raise+1][0][(y*width+x)*3+c]
            *new_alphas[raise+1][0][y*width+x]/255.)*

            (1-foreground_ps[raise+1][0][(y*width+x)]-
            background_ps[raise+1][0][(y*width+x)]);
        }
      }
    }
    save_image(buffer, width, height, 3, tmp);
    delete[] tmp;

    snprintf(buffer, 100, RESULTS "new_foregrounds_%i.ppm", raise+1);
    save_image(buffer, width, height, 3, new_foregrounds[raise+1][0]);
    
    snprintf(buffer, 100, RESULTS "new_backgrounds_%i.ppm", raise+1);
    save_image(buffer, width, height, 3, new_backgrounds[raise+1][0]);
    
    snprintf(buffer, 100, RESULTS "new_alphas_%i.ppm", raise+1);
    save_image(buffer, width, height, 1, new_alphas[raise+1][0]);

    raise++;
  }

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
  delete[] new_alphas[9][0];
}
