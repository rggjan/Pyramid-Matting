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
#define DEBUG_TREE

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

  if (alpha_pointer != NULL)
    *alpha_pointer = (unsigned char)(alpha*255);
/*
  for (int c=0; c<3; c++) {
    B[c] = new_B[c];
    F[c] = new_F[c];
  }*/
  return quality;
}


inline void find_best_combination(unsigned char* original,
    unsigned char* foregrounds[10][2],
    unsigned char* backgrounds[10][2],
    double* foreground_ps[10][2],
    double* background_ps[10][2],
    double fg_ps,
    double bg_ps,
    unsigned char* foreground_sure,
    unsigned char* background_sure,
    int max_raise,
    int max_stretched,
    unsigned char* final_f,
    unsigned char* final_b,
    unsigned char* final_a){
  static int counter = 0;
  counter++;
  if (counter == 295731) {
    cout << "Original: " << endl;
    for (int c=0; c<3; c++) {
      cout << (int)original[c] << endl;
    }
  }

  int width=1;
  int height=1;
  int fg_x = 0;
  int fg_y = 0;
  int bg_x = 0;
  int bg_y = 0;

  int global_best_score = -1;

  for (int raise=0; raise<=max_raise; raise++) {
#ifdef DEBUG_TREE
    if (counter == 295731) {
      cout << "\n==============\nglobal fg: (" << fg_x << "/" << fg_y << "), bg: (" << bg_x << "/" << bg_y << ")\n";
    }
#endif

    width *= 2;
    fg_x *= 2;
    bg_x *= 2;

    int local_best_score = -1;
    int best_fxdiff = 0;
    int best_bxdiff = 0;
    for (int fxdiff = 0; fxdiff<=1; fxdiff++) {
      for (int bxdiff = 0; bxdiff <= 1; bxdiff++) {
        unsigned char alpha;
        unsigned char* ft = &(foregrounds[raise][1][(fg_y*width+fg_x+fxdiff)*3]);
        unsigned char* bt = &(backgrounds[raise][1][(bg_y*width+bg_x+bxdiff)*3]);
        
        unsigned char test_f[3];
        unsigned char test_b[3];
        for (int c=0; c<3; c++) {
          test_f[c] = foreground_sure[c]*fg_ps + (1-fg_ps)*ft[c];
          test_b[c] = background_sure[c]*bg_ps + (1-bg_ps)*bt[c];
        }
        double score = projection(test_f, test_b, original, &alpha);
#ifdef DEBUG_TREE
        if (counter == 295731) {
          cout << "\nfg: (" << fg_x+fxdiff << "/" << fg_y << "), bg: (" << bg_x+bxdiff << "/" << bg_y << ")\n";
          cout << "test_f\ttest_b" << endl;
          for (int c=0; c<3; c++) {
            cout << (int)test_f[c] << "\t" << (int)test_b[c] << endl;
          }
          cout << "alpha/score: " << (int)alpha << "/" << score << endl;
        }
#endif
        if (score < local_best_score || local_best_score == -1) {
          local_best_score = score;
          best_fxdiff = fxdiff;
          best_bxdiff = bxdiff;

          if (score < global_best_score || global_best_score == -1) {
            global_best_score = score;
            *final_a = alpha;
            for (int c=0; c<3; c++) {
              final_f[c] = test_f[c];
              final_b[c] = test_b[c];
            }
          }
        }
      }
    }

    fg_x += best_fxdiff;
    bg_x += best_bxdiff;
    if (raise < max_raise || max_stretched > 0) {
#ifdef DEBUG_TREE
      if (counter == 295731) {
        cout << "\n==============\nglobal fg: (" << fg_x << "/" << fg_y << "), bg: (" << bg_x << "/" << bg_y << ")\n";
      }
#endif

      // Stretch other direction
      height *= 2;
      bg_y *= 2;
      fg_y *= 2;

      int local_best_score = -1;
      int best_fydiff;
      int best_bydiff;
      for (int fydiff = 0; fydiff<=1; fydiff++) {
        for (int bydiff = 0; bydiff <= 1; bydiff++) {
          unsigned char alpha;
          unsigned char* ft = &(foregrounds[raise+1][0][((fg_y+fydiff)*width+fg_x)*3]);
          unsigned char* bt = &(backgrounds[raise+1][0][((bg_y+bydiff)*width+bg_x)*3]);

          unsigned char test_f[3];
          unsigned char test_b[3];
          for (int c=0; c<3; c++) {
            test_f[c] = foreground_sure[c]*fg_ps + (1-fg_ps)*ft[c];
            test_b[c] = background_sure[c]*bg_ps + (1-bg_ps)*bt[c];
          }
          double score = projection(test_f, test_b, original, &alpha);
#ifdef DEBUG_TREE
        if (counter == 295731) {
          cout << "\nfg: (" << fg_x << "/" << fg_y+fydiff << "), bg: (" << bg_x << "/" << bg_y+bydiff << ")\n";
          cout << "test_f\ttest_b" << endl;
          for (int c=0; c<3; c++) {
            cout << (int)test_f[c] << "\t" << (int)test_b[c] << endl;
          }
          cout << "alpha/score: " << (int)alpha << "/" << score << endl;
        }
#endif
          if (score < local_best_score || local_best_score == -1) {
            local_best_score = score;
            best_fydiff = fydiff;
            best_bydiff = bydiff;

            if (score < global_best_score || global_best_score == -1) {
              global_best_score = score;
              *final_a = alpha;
              for (int c=0; c<3; c++) {
                final_f[c] = test_f[c];
                final_b[c] = test_b[c];
              }
            }
          }
        }
      }
      
      fg_y += best_fydiff;
      bg_y += best_bydiff;
    }
  }
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

  unsigned char* originals[9+1][2];
  unsigned char* new_alphas[9+1][2];
  unsigned char* foregrounds[9+1][2];
  unsigned char* new_foregrounds[9+1][2];
  double* foreground_ps[9+1][2];
  unsigned char* backgrounds[9+1][2];
  unsigned char* new_backgrounds[9+1][2];
  double* background_ps[9+1][2];

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

        find_best_combination(original1,
            new_foregrounds,
            new_backgrounds,
            foreground_ps,
            background_ps,
            foreground_ps1[0],
            background_ps1[0],
            foreground1,
            background1,
            raise, 0,
            f1, b1, a1);
        
        find_best_combination(original2,
            new_foregrounds,
            new_backgrounds,
            foreground_ps,
            background_ps,
            foreground_ps2[0],
            background_ps2[0],
            foreground2,
            background2,
            raise, 0,
            f2, b2, a2);

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


        find_best_combination(original1,
            new_foregrounds,
            new_backgrounds,
            foreground_ps,
            background_ps,
            foreground_ps1[0],
            background_ps1[0],
            foreground1,
            background1,
            raise, 1,
            f1, b1, a1);
        if (raise == 8 && x == 410 && y == 64) {
          exit(0);
        }
        
        find_best_combination(original2,
            new_foregrounds,
            new_backgrounds,
            foreground_ps,
            background_ps,
            foreground_ps2[0],
            background_ps2[0],
            foreground2,
            background2,
            raise, 1,
            f2, b2, a2);
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
