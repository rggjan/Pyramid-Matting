#include <stdio.h>
#include <stdlib.h>


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

int main() {
  int width = 512;
  int height = 512;

  unsigned char* original = load_image("test.ppm", width, height, 3);
  unsigned char* mask = load_image("trimap.pnm", width, height, 1);
  unsigned char* final = new unsigned char[width*height*3];

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      for (int color=0; color<3; color++) {
        int index = y*width*3+x*3+color;
        if (mask[y*width+x] == 0)
          final[index] = original[index];
        else
          final[index] = 0;
      }
    }
  }

  save_image("final.ppm", width, height, 3, final);

  delete[] original;
  delete[] mask;
  delete[] final;
}
