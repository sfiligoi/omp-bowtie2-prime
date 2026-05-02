/*
 * Simple exerciser of aligner_swsse_ee_u8.cpp
 */

#define SWSSE_INLINE_ONLY

#include <omp.h>

#include "../aligner_swsse_ee_u8.cpp"

#define MAX_PB_EL   128*(32/NBYTES_PER_REG)
#define MAX_RF_EL   192*(32/NBYTES_PER_REG)
#define MAX_MAT_EL  4096*(32/NBYTES_PER_REG)

#ifndef VEC_WIDTH
#define VEC_WIDTH 64
#endif

#define tread(data,size,n, file) if (fread(data,size,n,file)!=n) {fprintf(stderr, "Read failed\n"); return 3;}

typedef struct {
  uint32_t vec_els[VEC_WIDTH];
} TestScoreVec;

typedef struct {
  char vec_els[VEC_WIDTH];
} TestRFVec;

int load_data(int nels,
		size_t nrow[],
		size_t iter[],
		size_t colstride[],
		size_t lastWordIdx[],
		int32_t minsc[],
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

void transpose_data(int nels,
		const size_t  iter[],
                const uint8_t profbuf[], // the semantics only make sense if SSERegI==uint8_t
		const char    rf[],
		TestScoreVec vecscorebuf[],
		TestRFVec    vecrf[]
		) {
#pragma omp parallel for
	for (int i=0; i<nels; i++) {
		int64_t nblock = i/VEC_WIDTH;
		int vec_id = i%VEC_WIDTH;

		// use packed representation

		const int my_iter = iter[i];
		// not all of the PB is used... only in multiples of iter
		const int MAX_PB_ITER_EL = (MAX_PB_EL/2)/my_iter;
		const int MAX_PB_VEC_EL = MAX_PB_EL/4;

		const uint16_t *my_profbuf2 = (const uint16_t *) (profbuf+i*size_t(MAX_PB_EL));  // always pairs
		TestScoreVec* my_vecscorebuf = vecscorebuf+nblock*size_t(MAX_PB_VEC_EL);

		size_t ve=0;
		for (int r=0; r<MAX_PB_ITER_EL; r++) {
		   for (int j=0; j<(my_iter-1); j+=2) {
			const int e=r*my_iter+j;
			my_vecscorebuf[ve].vec_els[vec_id] = ((const uint32_t *)(my_profbuf2+e))[0]; // access both elements
			ve++;
		   }
		   if ((my_iter%2)!=0) {
			const int e=r*my_iter+(my_iter-1);
			my_vecscorebuf[ve].vec_els[vec_id] = ((const uint32_t *)(my_profbuf2+e))[0]; // even if I read past the end, it's ok here
			ve++;
		   }
		}

		for (int e=0; e<MAX_RF_EL; e++) {
			vecrf[nblock*size_t(MAX_RF_EL)+e].vec_els[vec_id] = rf[i*size_t(MAX_RF_EL)+e];
		}
	}
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

template<typename TP, typename TR>
int align_ee8_one(const int el, // for debuggging purpose only
		const int    vec_id,
		const size_t nrow,
		const size_t iter,
		const size_t colstride,
		const size_t lastWordIdx,
		const int32_t minsc,
		const size_t rfd,
                const TP     *profbuf,
		const TR     *rf,
		const uint8_t gaps[],
                SSERegI          *mat,
		DpBtCandidate    *btncand,
                const int32_t ref_lrmax,
                const int32_t ref_btnfilled) {
#ifndef SSE_FAST_SCALAR
	uint16_t btnfilled = 0;
	const EEU8_TCScore lrmax = EEU8_alignNucleotides<uint16_t>(profbuf, rf, rfd,
					mat,
                                        iter, colstride, lastWordIdx,
					minsc, nrow,
					btncand, btnfilled,
					gaps[0],gaps[1],gaps[2],gaps[3]);
#else
	// Note: The vast majority of reads have iter==151
	if (iter!=151) return 0;  // TODO: We may properly use the right template function, or call slow align directly
	const EEU8_TCScore lrmax = EEU8_alignNucleotidesScalar<uint16_t,151>(vec_id, profbuf, rf, rfd,
					nrow,
					gaps[0],gaps[1],gaps[2],gaps[3]);
#endif
	int nerrs = 0;
	if (int(ref_lrmax) != int(lrmax)) nerrs++;
#ifndef SSE_FAST_SCALAR
	if (int(ref_btnfilled) != int(btnfilled)) nerrs++;
#else
	TAlScore sc = (TAlScore)(lrmax - 0xff);
	if ((ref_btnfilled!=0) != (sc>=minsc)) nerrs++;
#endif
#ifndef NO_CHECK_PRINT
	if (int(ref_lrmax) != int(lrmax)) fprintf(stderr, "[%i] MISMATCH in lrmax (%i != %i)\n",el,int(lrmax), int(ref_lrmax));
#ifndef SSE_FAST_SCALAR
	if (int(ref_btnfilled) != int(btnfilled)) fprintf(stderr, "[%i] MISMATCH in ref_btnfilled (%i != %i)\n",el,int(btnfilled), int(ref_btnfilled));
#else
	if ((ref_btnfilled!=0) != (sc>=minsc)) fprintf(stderr, "[%i] MISMATCH in ref_btnfilled (%i) vs sc %i minsc %i)\n",el,int(ref_btnfilled), int(sc), int(minsc));
#endif
#endif
	// TODO: In case of SSE_FAST_SCALAR, we should also likely want to test the complete align on the elements that need it
	//       But that's probably better done in a separate test
	return nerrs;
}


template<typename TP, typename TR>
void align_ee8(const int npar, const int nels,
		const size_t nrow[],
		const size_t iter[],
		const size_t colstride[],
		const size_t lastWordIdx[],
		const int32_t minsc[],
		const size_t rfd[],
                const TP     profbuf[],
		const TR     rf[],
		const uint8_t gaps[],
                SSERegI          mat[], 
		DpBtCandidate    btncand[],
                const int32_t ref_lrmax[],
                const int32_t ref_btnfilled[]) {
   int nerrs = 0;
#ifdef OMPGPU
#pragma omp target teams distribute parallel for reduction(+:nerrs) num_teams(npar) num_threads(VEC_WIDTH*8)
#else
#pragma omp parallel for reduction(+:nerrs) num_threads(npar) schedule(static,VEC_WIDTH)
#endif
   for (int i=0; i<nels; i++) {
	 int64_t nblock = i/VEC_WIDTH;
	 int vec_id = i%VEC_WIDTH;
#ifndef SSE_FAST_SCALAR 
	 // unpacked/raw
	 constexpr uint64_t maxPB = MAX_PB_EL;
	 const int p = omp_get_thread_num(); // only supported on CPU
#else
	 // packed
	 constexpr uint64_t maxPB = MAX_PB_EL/4;
	 // unusued
	 const int p = 0;
#endif
         nerrs+= align_ee8_one(i, vec_id,
		nrow[i], iter[i], colstride[i], lastWordIdx[i],minsc[i], rfd[i],
                profbuf+nblock*maxPB,rf+nblock*size_t(MAX_RF_EL),gaps+4*i,
		mat+p*size_t(MAX_MAT_EL),btncand+p*size_t(MAX_RF_EL),
		ref_lrmax[i],ref_btnfilled[i]);
   }

   if (nerrs!=0) {
      fprintf(stderr, "FAILED matching results %i times!\n", nerrs);
   } else {
      fprintf(stderr, "SUCCESS, all results matched.\n");
   }

}

int load_and_align_ee8(const int npar, const int nels) {
   int32_t *ref_lrmax = new int32_t[nels];
   int32_t *ref_btnfilled = new int32_t[nels];
   SSERegI *raw_profbuf = new SSERegI[size_t(nels)*MAX_PB_EL];
#ifndef SSE_FAST_SCALAR 
   SSERegI *mat = new SSERegI[size_t(npar)*MAX_MAT_EL];
#else
   // matrix won't be used in fast scalar
   SSERegI *mat = nullptr;
#endif
   char    *raw_rf = new char[size_t(MAX_RF_EL)*nels];
   size_t  *nrow = new size_t[nels];
   size_t  *iter = new size_t[nels];
   size_t  *colstride = new size_t[nels];
   size_t  *lastWordIdx = new size_t[nels];
   int32_t  *minsc = new int32_t[nels];
   size_t  *rfd = new size_t[nels];
   uint8_t *gaps = new uint8_t[4*nels];
   DpBtCandidate    *btncand = new DpBtCandidate[size_t(MAX_RF_EL)*npar];

#ifdef SSE_FAST_SCALAR
   TestScoreVec *vecscorebuf = new TestScoreVec[size_t((nels+VEC_WIDTH-1)/VEC_WIDTH)*(MAX_PB_EL/4)];
   TestRFVec  *vecrf = new TestRFVec[size_t(MAX_RF_EL)*((nels+VEC_WIDTH-1)/VEC_WIDTH)];
#endif


   auto t1 = std::chrono::high_resolution_clock::now();
   // test tool, don't worry about perfect cleanup
   if (load_data(nels,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		raw_profbuf,raw_rf,gaps) !=0 ) return 1;
   if (load_results(nels,
		ref_lrmax, ref_btnfilled) !=0 ) return 1;
#ifdef SSE_FAST_SCALAR
   // pass transoposed version into the benchmarking section
   transpose_data(nels,
		iter,
		raw_profbuf,raw_rf, vecscorebuf,vecrf);
   const auto *rf = vecrf;
   const auto *profbuf = vecscorebuf;
#else
   // keep the original versions everywhere
   const auto *rf = raw_rf;
   const auto *profbuf = raw_profbuf;
#endif
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
