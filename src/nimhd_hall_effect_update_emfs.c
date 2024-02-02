//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void HallEffect_UpdateEmfs_cpu(){

//<USER_DEFINED>
  INPUT(EmfxH);
  INPUT(EmfyH);
  INPUT(EmfzH);
  INPUT(Emfx);
  INPUT(Emfy);
  INPUT(Emfz);
  OUTPUT(Emfx);
  OUTPUT(Emfy);
  OUTPUT(Emfz);
//<\USER_DEFINED>

//<EXTERNAL>
  real* emfx = Emfx->data->field_cpu;
  real* emfy = Emfy->data->field_cpu;
  real* emfz = Emfz->data->field_cpu;
  real* emfxH = EmfxH->data->field_cpu;
  real* emfyH = EmfyH->data->field_cpu;
  real* emfzH = EmfzH->data->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
//<\INTERNAL>


//<MAIN_LOOP>
  for (k=1; k<size_z; k++) {
    for (j=1; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
//<#>
	ll = l;

	emfx[ll] += emfxH[ll];
	emfy[ll] += emfyH[ll];
       	emfz[ll] += emfzH[ll];
//<\#>
      }
    }
  }
//<\MAIN_LOOP>
}
