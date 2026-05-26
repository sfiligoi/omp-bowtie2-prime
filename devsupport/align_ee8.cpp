/*
 * Simple exerciser of aligner_swsse_ee_u8.cpp
 */

#define SWSSE_INLINE_ONLY

#include "../aligner_swsse_ee_u8.cpp"

#define MAX_PB_EL   128*(32/NBYTES_PER_REG)
#define MAX_RF_EL   192*(32/NBYTES_PER_REG)
#define MAX_MAT_EL  4096*(32/NBYTES_PER_REG)

#define tread(data,size,n, file) if (fread(data,size,n,file)!=n) {fprintf(stderr, "Read failed\n"); return 3;}

#define SKIPPED_LRMAX -2000

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

int align_ee8_one(const int el, // for debuggging purpose
		const size_t nrow,
		const size_t iter,
		const size_t colstride,
		const size_t lastWordIdx,
		const int32_t minsc,
		const size_t rfd,
                const SSERegI *profbuf,
		const char    *rf,
		const uint8_t gaps[],
                SSERegI       *mat,
		DpBtCandidate *btncand,
                int32_t       &computed_lrmax, // out
                const int32_t ref_lrmax,
                const int32_t ref_btnfilled) {
#ifndef PRE_LR_SCALAR
	uint16_t btnfilled = 0;
	const EEU8_TCScore lrmax = EEU8_alignNucleotides<uint16_t>(profbuf, rf, rfd,
					mat,
                                        iter, colstride, lastWordIdx,
					minsc, nrow,
					btncand, btnfilled,
					gaps[0],gaps[1],gaps[2],gaps[3]);
#else
	// Note: The vast majority of reads have iter==151
	if (iter!=151) {
		computed_lrmax = SKIPPED_LRMAX; // special value to signal we did not actually compute it
		return 0;  // TODO: We may properly use the right template function, or call slow align directly
	}
#ifndef PBMASK
	const EEU8_TCScore lrmax = EEU8_alignNucleotidesLRScalar<uint16_t,151>(profbuf, rf, rfd,
					nrow,
					gaps[0],gaps[1],gaps[2],gaps[3]);

#else
	const EEU8_TCScore lrmax = EEU8_alignNucleotidesLRM11Scalar<uint16_t,151>(profbuf, rf, rfd,
					nrow,
					gaps[0],gaps[1],gaps[2],gaps[3]);
#endif // PBMASK
#endif // PR_LR_SCALAR
	int nerrs = 0;
	computed_lrmax = lrmax;
	if (int(ref_lrmax) != int(lrmax)) nerrs++;
#ifndef PRE_LR_SCALAR
	if (int(ref_btnfilled) != int(btnfilled)) nerrs++;
#else
	TAlScore sc = (TAlScore)(lrmax - 0xff);
	if ((ref_btnfilled!=0) != (sc>=minsc)) nerrs++;
#endif
#ifndef NO_CHECK_PRINT
	if (int(ref_lrmax) != int(lrmax)) fprintf(stderr, "[%i] MISMATCH in lrmax (%i != %i)\n",el,int(lrmax), int(ref_lrmax));
#ifndef PRE_LR_SCALAR
	if (int(ref_btnfilled) != int(btnfilled)) fprintf(stderr, "[%i] MISMATCH in ref_btnfilled (%i != %i)\n",el,int(btnfilled), int(ref_btnfilled));
#else
	if ((ref_btnfilled!=0) != (sc>=minsc)) fprintf(stderr, "[%i] MISMATCH in ref_btnfilled (%i) vs sc %i minsc %i)\n",el,int(ref_btnfilled), int(sc), int(minsc));
#endif
#endif
	// TODO: In case of PRE_LR_SCALAR, we should also likely want to test the complete align on the elements that need it
	//       But that's probably better done in a separate test
	return nerrs;
}


void align_ee8(const int npar, const int nels,
		const size_t nrow[],
		const size_t iter[],
		const size_t colstride[],
		const size_t lastWordIdx[],
		const int32_t minsc[],
		const size_t rfd[],
                const SSERegI profbuf[],
		const char    rf[],
		const uint8_t gaps[],
                SSERegI       mat[], 
		DpBtCandidate btncand[],
                int32_t       computed_lrmax[], // out
                const int32_t ref_lrmax[],
                const int32_t ref_btnfilled[]) {
   const int elp_pp = (nels+ (npar-1))/npar; // round up
   int nerrs = 0;
#ifdef OMPGPU
#pragma omp target teams distribute reduction(+:nerrs)
#else
#pragma omp parallel for schedule(dynamic) reduction(+:nerrs)
#endif
   for (int p=0; p<npar; p++) {
      const int iend = std::min((p+1)*elp_pp,nels);
#if defined(PRE_LR_SCALAR)
	// no mat or btncand, just pass NULLs around 
	SSERegI       *my_mat = nullptr;
	DpBtCandidate *my_btncand = nullptr;
#elif !defined(OMPGPU)
	// the following loop is sequential, so can reuse the buffer
	SSERegI       *my_mat = mat+p*size_t(MAX_MAT_EL);
	DpBtCandidate *my_btncand = btncand+p*size_t(MAX_RF_EL);
#endif
#ifdef OMPGPU
#pragma omp parallel for reduction(+:nerrs)
#endif
      for (int i=p*elp_pp; i<iend; i++) {
#if defined(OMPGPU) && !defined(PRE_LR_SCALAR)
	// the loop is parallel, so we must use a different buffer for each element
	SSERegI       *my_mat = mat+i*size_t(MAX_MAT_EL);
	DpBtCandidate *my_btncand = btncand+i*size_t(MAX_RF_EL);
#endif
         nerrs+= align_ee8_one(i,
		nrow[i], iter[i], colstride[i], lastWordIdx[i],minsc[i], rfd[i],
                profbuf+i*size_t(MAX_PB_EL),rf+i*size_t(MAX_RF_EL),gaps+4*i,
		my_mat,my_btncand,
		computed_lrmax[i],ref_lrmax[i],ref_btnfilled[i]);
      }
   }

   if (nerrs!=0) {
      fprintf(stderr, "FAILED matching results %i times!\n", nerrs);
   } else {
      fprintf(stderr, "SUCCESS, all results matched.\n");
   }

}

int load_and_align_ee8(const int npar, const int nels) {
   int32_t *computed_lrmax = new int32_t[nels];
   int32_t *ref_lrmax = new int32_t[nels];
   int32_t *ref_btnfilled = new int32_t[nels];
   SSERegI *profbuf = new SSERegI[size_t(nels)*MAX_PB_EL];
#if defined(PRE_LR_SCALAR)
   // no mat or btncand, just pass NULLs around 
   SSERegI        *mat = nullptr;
   DpBtCandidate *btncand = nullptr;
#elif !defined(OMPGPU)
   // one per thread is enough on the CPU
   SSERegI       *mat = new SSERegI[size_t(npar)*MAX_MAT_EL];
   DpBtCandidate *btncand = new DpBtCandidate[size_t(MAX_RF_EL)*npar];
#else
   // no shortcuts on the GPU
   SSERegI       *mat = new SSERegI[size_t(nels)*MAX_MAT_EL];
   DpBtCandidate *btncand = new DpBtCandidate[size_t(MAX_RF_EL)*nels];
#endif
   char    *rf = new char[size_t(MAX_RF_EL)*nels];
   size_t  *nrow = new size_t[nels];
   size_t  *iter = new size_t[nels];
   size_t  *colstride = new size_t[nels];
   size_t  *lastWordIdx = new size_t[nels];
   int32_t  *minsc = new int32_t[nels];
   size_t  *rfd = new size_t[nels];
   uint8_t *gaps = new uint8_t[4*nels];

   auto t1 = std::chrono::high_resolution_clock::now();
   // test tool, don't worry about perfect cleanup
   if (load_data(nels,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		profbuf,rf,gaps) !=0 ) return 1;
   if (load_results(nels,
		ref_lrmax, ref_btnfilled) !=0 ) return 1;
   auto t2a = std::chrono::high_resolution_clock::now();
   for (int i=0; i<nels; i++) computed_lrmax[i] = -1; // invalidate
   auto t2b = std::chrono::high_resolution_clock::now();
   align_ee8(npar, nels,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		profbuf,rf,gaps,mat,btncand,
		computed_lrmax,ref_lrmax,ref_btnfilled);
   auto t3a = std::chrono::high_resolution_clock::now();
   {
     int nerrs = 0;
     int ncomputed = 0;
     for (int i=0; i<nels; i++) {
      if (computed_lrmax[i] != SKIPPED_LRMAX ) {
	ncomputed++;
	if (computed_lrmax[i] != ref_lrmax[i]) {
	   nerrs++;
#ifndef NO_CHECK_PRINT
	   fprintf(stderr, "[%i] MISMATCH in lrmax (%i != %i)\n",i,int(computed_lrmax[i]), int(ref_lrmax[i]));
#endif
	}
      }
     }
     fprintf(stderr, "INFO: Computed %i out of %i alignments, %i errors\n", ncomputed, nels, nerrs);
   }
   for (int i=0; i<nels; i++) computed_lrmax[i] = -1; // invalidate
   auto t3b = std::chrono::high_resolution_clock::now();
   align_ee8(npar, nels,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		profbuf,rf,gaps,mat,btncand,
		computed_lrmax, ref_lrmax, ref_btnfilled);
   auto t4 = std::chrono::high_resolution_clock::now();
   {
     int nerrs = 0;
     int ncomputed = 0;
     int npartials = 0;
     for (int i=0; i<nels; i++) {
      if (computed_lrmax[i] != SKIPPED_LRMAX ) {
	ncomputed++;
	if (computed_lrmax[i] != ref_lrmax[i]) {
	   nerrs++;
#ifndef NO_CHECK_PRINT
	   fprintf(stderr, "[%i] MISMATCH in lrmax (%i != %i)\n",i,int(computed_lrmax[i]), int(ref_lrmax[i]));
#endif
	}
#if defined(PRE_LR_SCALAR)
	// The fast method does not provide all info is sc too large
	TAlScore sc = (TAlScore)(computed_lrmax[i] - 0xff);
	if (sc>=minsc[i]) npartials++;
#endif
      }
     }
     fprintf(stderr, "INFO: Computed %i out of %i alignments, %i npartials, %i errors\n", ncomputed, nels, npartials, nerrs);
     if ((ncomputed-npartials+nerrs)!=nels) fprintf(stderr, "INFO: A 2nd round for %i els still needed\n",nels-(ncomputed-npartials+nerrs));
   }

   auto time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2a - t1);
   auto time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t3a - t2b);
   auto time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3b);
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
   delete[] computed_lrmax;
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
