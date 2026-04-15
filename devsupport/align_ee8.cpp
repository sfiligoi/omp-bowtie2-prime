/*
 * Simple exerciser of aligner_swsse_ee_u8.cpp
 */

// assumed CDNA where warp size is 64. 32 on Nvidia and RDNA
#ifndef WARP_SIZE
#define WARP_SIZE 64
#endif

#ifndef SSE_SCALAR
 #ifndef E_PER_WARP
  #define E_PER_WARP 2
 #endif
 #ifndef LANE_SIZE
  #define LANE_SIZE 32 // number of lanes used in this
 #endif
#else // enabled SSE_SCALAR
 #define E_PER_WARP 1
 #define LANE_SIZE 1
#endif // SSE_SCALAR

#define SWSSE_INLINE_ONLY

#include "../aligner_swsse_ee_u8.cpp"

#define MAX_PB_EL   128*(32/NBYTES_PER_REG)
#define MAX_RF_EL   192*(32/NBYTES_PER_REG)
#define MAX_MAT_EL  4096*(32/NBYTES_PER_REG)

#ifdef HIP_KERNELS
#include <hip/hip_runtime.h>

#define HIP_CHECK(expression)                  \
{                                              \
    const hipError_t status = expression;      \
    if(status != hipSuccess){                  \
        std::cerr << "HIP error "              \
                  << status << ": "            \
                  << hipGetErrorString(status) \
                  << " at " << __FILE__ << ":" \
                  << __LINE__ << std::endl;    \
    }                                          \
}
#endif

#define tread(data,size,n, file) if (fread(data,size,n,file)!=n) {fprintf(stderr, "Read failed\n"); return 3;}

int load_data(int nels,
		size_t nrow[],
		size_t iter[],
		size_t colstride[],
		size_t lastWordIdx[],
		size_t minsc[],
		size_t rfd[],
                SSERegI profbuf[],
		char    rf[],
		uint8_t gaps[]
		) {
   char fname[256];
   sprintf(fname,"EEU8_an_in_%i.bin",int(sizeof(SSERegI)));
   FILE *file = fopen(fname, "rb");
   if (file==NULL) {
      fprintf(stderr, "Could not open %s\n",fname);
      return 2;
   }
   
   for (int i=0; i<nels; i++) {
      uint32_t magic = 0;
      uint32_t isize = 0;
      uint32_t cnt = 9999999;
      size_t mat_size = 0;
      size_t profbuf_size = 0;
      tread(&magic,sizeof(uint32_t), 1, file);
      if (magic!=0x12943ef2) {fprintf(stderr, "Wrong magic\n"); return 4;}
      tread(&isize,sizeof(uint32_t), 1, file);
      if (isize != sizeof(SSERegI)) {fprintf(stderr, "Wrong isize\n"); return 4;}
      tread(&cnt,sizeof(uint32_t), 1, file);
      if (cnt != i) {fprintf(stderr, "Wrong cnt\n"); return 4;}
      tread(nrow+i,sizeof(size_t), 1, file);
      tread(iter+i,sizeof(size_t), 1, file);
      tread(colstride+i,sizeof(size_t), 1, file);
      tread(lastWordIdx+i,sizeof(size_t), 1, file);
      tread(minsc+i,sizeof(int32_t), 1, file);
      tread(&mat_size,sizeof(size_t), 1, file);
      if (mat_size>MAX_MAT_EL) {fprintf(stderr, "Mat size too big %i\n", int(mat_size)); return 4;}
      tread(&profbuf_size,sizeof(size_t), 1, file);
      if (profbuf_size>MAX_PB_EL) {fprintf(stderr, "Profbuf size too big %i\n", int(profbuf_size)); return 4;}
      tread(profbuf+(i*size_t(MAX_PB_EL)), sizeof(SSERegI), profbuf_size, file);
      tread(rfd+i,sizeof(size_t), 1, file);
      if (rfd[i]>MAX_RF_EL) {fprintf(stderr, "RF size too big %i\n", int(rfd[i])); return 4;}
      tread(rf+(i*size_t(MAX_RF_EL)), sizeof(char), rfd[i], file);
      tread(gaps+i*4, sizeof(int8_t), 4, file);
   }

   fclose(file);
   return 0;
}

int load_results(int nels,
		int32_t lrmax[],
		int32_t btnfilled[]) {
   char fname[256];
   sprintf(fname,"EEU8_an_res_%i.bin",int(sizeof(SSERegI)));
   FILE *file = fopen(fname, "rb");
   if (file==NULL) {
      fprintf(stderr, "Could not open %s\n",fname);
      return 2;
   }
   
   for (int i=0; i<nels; i++) {
      uint32_t cnt = 9999999;
      tread(&cnt,sizeof(uint32_t), 1, file);
      if (cnt != i) {fprintf(stderr, "Wrong cnt\n"); return 4;}
      tread(lrmax+i,sizeof(int32_t), 1, file);
      tread(btnfilled+i,sizeof(int32_t), 1, file);
   }

   fclose(file);
   return 0;
}

#ifdef HIP_KERNELS
/*
 * Wrapper for EEU8_alignNucleotides_HIP
 * But handles the offsets and batch sizes
 *
 */
__global__
void EEU8_alignNucleotidesBatch_HIP(int npar, int nels,
		const size_t nrow_[],
		const size_t iter_[],
		const size_t colstride_[],
		const size_t lastWordIdx_[],
		const size_t minsc_[],
		const size_t rfd_[],
		const SSERegI profbuf_[],
		const char    rf_[],
		const uint8_t gaps_[],
		SSERegI          mat_[],
		DpBtCandidate    btncand_[],
		int32_t  lrmax_[]) {

    const int p = blockIdx.x; // parallel task (1 per gpu block)
    const int elp_pp = (nels+ (npar-1))/npar; // round up
    int lane_id_b = (threadIdx.x / LANE_SIZE);
    const int lane_id = threadIdx.x % LANE_SIZE;// used for output and debugging in this wrapper function.

    // Equivalently the grid dimension

    const int iend = std::min((p+1)*elp_pp,nels);
    for (int i=p*elp_pp+lane_id_b; i<iend; i+=E_PER_WARP) {
       // leave each block to parallelize
       // where each block does 32 lanes.

       // Only curBatchSize amount indexed, but this is to compile and is calculated in the alignNucleotides kernel
       uint8_t btnfilled;
       uint32_t lrmax;


       const size_t nrow        = nrow_[i];
       const size_t iter        = iter_[i];
       const size_t colstride   = colstride_[i];
       const size_t lastWordIdx = lastWordIdx_[i];
       const size_t minsc       = minsc_[i];
       const size_t rfd         = rfd_[i];

       const SSERegI* profbuf   = profbuf_ + (i * size_t(MAX_PB_EL)); // why is this needed here?
       const char* rf           = rf_      + (i * size_t(MAX_RF_EL));
       const uint8_t* gaps      = gaps_    + (i * 4);

       // NOTE: Use globalIdx here too, otherwise every block in the batch 
       // writes to the same 'p'... not used yet
       SSERegI* mat             = mat_ + ((E_PER_WARP*p+lane_id_b) * size_t(MAX_MAT_EL));
       DpBtCandidate* btncand   = btncand_ + ((E_PER_WARP*p+lane_id_b) * size_t(MAX_RF_EL));


       // Updates btnfilled and lrmax
       EEU8_alignNucleotides_HIP(
	     (uint8_t*)profbuf, rf, rfd,
	     (uint8_t*)mat,
	     iter, colstride, lastWordIdx,
	     minsc, nrow,
	     btncand, &btnfilled,
	     gaps[0], gaps[1], gaps[2], gaps[3],
	     &lrmax
	     );

       // write back local lrmax to global memory
       __syncthreads();
       if(lane_id == 0) {
	  lrmax_[i] = lrmax;
       }
       __syncthreads();

    }
}


/* Function exists to count errors and handle some of the separation on batching
 * @returns number of errors; zero if checking is disabled
 */
int align_ee8_batchLaunch(const int npar, const int nels,
		const size_t nrow_[],
		const size_t iter_[],
		const size_t colstride_[],
		const size_t lastWordIdx_[],
		const size_t minsc_[],
		const size_t rfd_[],
                const SSERegI profbuf_[],
		const char    rf_[],
		const uint8_t gaps_[],
                SSERegI          mat_[],
		DpBtCandidate    btncand_[],
                const int32_t ref_lrmax_[],
                const int32_t ref_btnfilled_[]) {

      static_assert(WARP_SIZE % E_PER_WARP == 0, "WARP_SIZE must be divisible by E_PER_WARP");
      constexpr int lanes = LANE_SIZE*E_PER_WARP;
      static_assert(lanes <= WARP_SIZE, "Number of total names used in a block must be <= Hardware Warp Size");

      int32_t* lrmax = (int32_t*)malloc(nels*sizeof(int32_t));
      int nerrs = 0;

      // Handle all the proper indexing in this call.
      EEU8_alignNucleotidesBatch_HIP<<<npar, lanes>>>(npar, nels,
	    nrow_, iter_, colstride_, lastWordIdx_, minsc_, rfd_,
	    profbuf_, rf_, gaps_,
	    mat_, btncand_,
	    lrmax);
      HIP_CHECK(hipDeviceSynchronize());

#ifndef NO_CHECK_PRINT
      for(int i = 0; i<nels; i++){
	 if (int(ref_lrmax_[i]) != lrmax[i]) {
	    printf("[%i] MISMATCH in lrmax (%i != %i)\n",i,int(lrmax[i]), int(ref_lrmax_[i]));
	    nerrs++;
         } 
      }
#endif
      free(lrmax);
      return nerrs;
}
#endif

int align_ee8_one(const int el, // for debuggging purpose
		const size_t nrow,
		const size_t iter,
		const size_t colstride,
		const size_t lastWordIdx,
		const size_t minsc,
		const size_t rfd,
                const SSERegI *profbuf,
		const char    *rf,
		const uint8_t gaps[],
                SSERegI          *mat,
		DpBtCandidate    *btncand,
                const int32_t ref_lrmax,
                const int32_t ref_btnfilled) {

	uint16_t btnfilled = 0;
	const EEU8_TCScore lrmax = EEU8_alignNucleotides<uint16_t>(profbuf, rf, rfd,
					mat,
                                        iter, colstride, lastWordIdx,
					minsc, nrow,
					btncand, btnfilled,
					gaps[0],gaps[1],gaps[2],gaps[3]);

	int nerrs = 0;
	if (int(ref_lrmax) != int(lrmax)) nerrs++;
#ifndef NO_CHECK_PRINT
	if (int(ref_lrmax) != int(lrmax)) fprintf(stderr, "[%i] MISMATCH in lrmax (%i != %i)\n",el,int(lrmax), int(ref_lrmax));
#endif
	return nerrs;
}


void align_ee8(const int npar, const int nels,
		const size_t nrow[],
		const size_t iter[],
		const size_t colstride[],
		const size_t lastWordIdx[],
		const size_t minsc[],
		const size_t rfd[],
                const SSERegI profbuf[],
		const char    rf[],
		const uint8_t gaps[],
                SSERegI          mat[], 
		DpBtCandidate    btncand[],
                const int32_t ref_lrmax[],
                const int32_t ref_btnfilled[]) {
   int nerrs = 0;

#ifndef HIP_KERNELS
   const int elp_pp = (nels+ (npar-1))/npar; // round up
#ifdef OMPGPU
#pragma omp target teams distribute parallel for reduction(+:nerrs)
#else
#pragma omp parallel for schedule(dynamic) reduction(+:nerrs)
#endif
   for (int p=0; p<npar; p++) {
      const int iend = std::min((p+1)*elp_pp,nels);
      for (int i=p*elp_pp; i<iend; i++) {
         nerrs+= align_ee8_one(i,
		nrow[i], iter[i], colstride[i], lastWordIdx[i],minsc[i], rfd[i],
                profbuf+i*size_t(MAX_PB_EL),rf+i*size_t(MAX_RF_EL),gaps+4*i,
		mat+p*size_t(MAX_MAT_EL),btncand+p*size_t(MAX_RF_EL),
		ref_lrmax[i],ref_btnfilled[i]);
      }
   }

#else
   // while previously npar meant splitting into nels/npar processes, this is SIMT
   // so a smaller batch size means running once on gpu and letting the GPU kernel launch
   // handle the indexing for SIMT
   nerrs = align_ee8_batchLaunch(npar, nels,
	 nrow, iter, colstride, lastWordIdx,minsc, rfd,
	 profbuf,rf,gaps,
	 mat,btncand,
	 ref_lrmax,ref_btnfilled);
#endif

   if (nerrs!=0) {
      fprintf(stderr, "FAILED matching results %i times!\n", nerrs);
   } else {
      fprintf(stderr, "SUCCESS, all results matched.\n");
   }
}

int load_and_align_ee8(const int npar, const int nels) {
   int32_t *ref_lrmax = new int32_t[nels];
   int32_t *ref_btnfilled = new int32_t[nels];
   SSERegI *profbuf = new SSERegI[size_t(nels)*MAX_PB_EL];
   SSERegI *mat = new SSERegI[size_t(npar)*MAX_MAT_EL*E_PER_WARP];
   char    *rf = new char[size_t(MAX_RF_EL)*nels];
   size_t  *nrow = new size_t[nels];
   size_t  *iter = new size_t[nels];
   size_t  *colstride = new size_t[nels];
   size_t  *lastWordIdx = new size_t[nels];
   size_t  *minsc = new size_t[nels];
   size_t  *rfd = new size_t[nels];
   uint8_t *gaps = new uint8_t[4*nels];
   DpBtCandidate    *btncand = new DpBtCandidate[size_t(MAX_RF_EL)*npar*E_PER_WARP];

   auto t1 = std::chrono::high_resolution_clock::now();
   // test tool, don't worry about perfect cleanup
   if (load_data(nels,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		profbuf,rf,gaps) !=0 ) return 1;
   if (load_results(nels,
		ref_lrmax, ref_btnfilled) !=0 ) return 1;
   auto t2 = std::chrono::high_resolution_clock::now();
   align_ee8(npar, nels,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		profbuf,rf,gaps,mat,btncand,
		ref_lrmax,ref_btnfilled);
   auto t3 = std::chrono::high_resolution_clock::now();
   align_ee8(npar, nels,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		profbuf,rf,gaps,mat,btncand,
		ref_lrmax, ref_btnfilled);
   auto t4 = std::chrono::high_resolution_clock::now();

   auto time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
   auto time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t3 - t2);
   auto time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
   fprintf(stderr, "load %.4f s cold align %.4f s hot align %.4f s\n", time_span1.count(), time_span2.count(), time_span3.count());

   delete[] btncand;
   delete[] gaps;
   delete[] rfd;
   delete[] minsc;
   delete[] lastWordIdx;
   delete[] colstride;
   delete[] iter;
   delete[] nrow;
   delete[] rf;
   delete[] mat;
   delete[] profbuf;
   delete[] ref_btnfilled;
   delete[] ref_lrmax;
   return 0;
}



int main(int argv, const char* argc[]) {
   if (argv<2) {
     fprintf(stderr, "Usage: test_align_ee8 <npar> <nels>\n");
     return 1;
   }
   int npar = atoi(argc[1]);
   int nels = atoi(argc[2]);

   return load_and_align_ee8(npar,nels);
}
