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

#ifndef CPUVECT
typedef SSERegI TMatBuf;
#else
// TODO, make it more dynameic
#define VECT_SIZE 16
typedef uint8_t TMatBuf __attribute__((ext_vector_type(VECT_SIZE)));
#endif

int load_data(int nels, int skip_blocks,
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
 
   for (int b=0; b<=skip_blocks; b++) {
    // read all elements, keep only the last pass
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
      if (cnt != (i+b*nels)) {fprintf(stderr, "Wrong cnt\n"); return 4;}
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
   }

   fclose(file);
   return 0;
}

int load_results(int nels, int skip_blocks,
		int32_t lrmax[],
		int32_t btnfilled[]) {
   char fname[256];
   sprintf(fname,"EEU8_an_res_%i.bin",int(sizeof(SSERegI)));
   FILE *file = fopen(fname, "rb");
   if (file==NULL) {
      fprintf(stderr, "Could not open %s\n",fname);
      return 2;
   }
   
   for (int b=0; b<=skip_blocks; b++) {
    // read all elements, keep only the last pass
    for (int i=0; i<nels; i++) {
      uint32_t cnt = 9999999;
      tread(&cnt,sizeof(uint32_t), 1, file);
      if (cnt != (i+b*nels)) {fprintf(stderr, "Wrong cnt\n"); return 4;}
      tread(lrmax+i,sizeof(int32_t), 1, file);
      tread(btnfilled+i,sizeof(int32_t), 1, file);
    }
   }

   fclose(file);
   return 0;
}

// Swap all data for two elements
void swap_data(size_t my, size_t other,
                size_t nrow[],
                size_t iter[],
                size_t colstride[],
                size_t lastWordIdx[],
                int32_t minsc[],
                size_t rfd[],
                SSERegI profbuf[],
                char    rf[],
                uint8_t gaps[],
                int32_t lrmax[],
                int32_t btnfilled[]
                ) {
	if (my==other) return; // nothing to do
	std::swap(nrow[my],nrow[other]);
	std::swap(iter[my],iter[other]);
	std::swap(colstride[my],colstride[other]);
	std::swap(lastWordIdx[my],lastWordIdx[other]);
	std::swap(minsc[my],minsc[other]);
	std::swap(rfd[my],rfd[other]);
	for (size_t i=0; i<MAX_PB_EL; i++) std::swap(profbuf[my*size_t(MAX_PB_EL)+i],profbuf[other*size_t(MAX_PB_EL)+i]);
	for (size_t i=0; i<MAX_RF_EL; i++) std::swap(rf[my*size_t(MAX_RF_EL)+i],rf[other*size_t(MAX_RF_EL)+i]);
	for (size_t i=0; i<4; i++) std::swap(gaps[my*4+i],gaps[other*4+i]);
	std::swap(lrmax[my],lrmax[other]);
	std::swap(btnfilled[my],btnfilled[other]);
}

// The accelerated algoritms need homogeneous alignment sizes
// Move anything that does not fit to the end
// Return how many homogeneous elements we have
int filter(int nels,
		size_t nrow[],
		size_t iter[],
		size_t colstride[],
		size_t lastWordIdx[],
		int32_t minsc[],
		size_t rfd[],
                SSERegI profbuf[],
		char    rf[],
		uint8_t gaps[],
		int32_t lrmax[],
		int32_t btnfilled[],
		const bool pass2_only
		) {
   uint32_t *gaps4 = (uint32_t*)gaps;
   bool is_first = true;
   size_t nrow_first;
   size_t colstride_first;
   size_t lastWordIdx_first;
   int32_t minsc_first;
   size_t rfd_first;
   uint32_t gaps_first;

   int filter_pass2_n = 0;
   int filter_iter_n = 0;
   int filter_nrow_n = 0;
   int filter_colstride_n = 0;
   int filter_lastWordIdx_n = 0;
   int filter_minsc_n = 0;
   int filter_rfd_n = 0;
   int filter_gaps_n = 0;
   int good_els = nels;
   int i = 0;
   while (i<good_els) {
      if (pass2_only) {
	const int32_t minlr = 
		(minsc[i]<=(-255)) 
		? 0                 // underflow
		: ((minsc[i]>0) 
			? 0xff      // overflow
			: (minsc[i] + 0xff)
		  );
        if (lrmax[i]<minlr) {
    	   // this one does not need a 2nd pass
	   swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	   good_els--;
	   filter_pass2_n++;
	   continue;
	}
      }
      if (iter[i]!=((151+(NBYTES_PER_REG-1))/NBYTES_PER_REG)) {
	// keep the most relevant one
	swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	good_els--;
	filter_iter_n++;
	continue;
      }
      // for all other values, we will default to the first we see
      if (is_first) {
	nrow_first = nrow[i];
	colstride_first = colstride[i];
	lastWordIdx_first = lastWordIdx[i];
 	minsc_first = minsc[i];
 	rfd_first = rfd[i];
   	gaps_first = gaps4[i];
	is_first = false;
      }
      if (rfd[i]!=rfd_first) {
	// keep the most relevant one
	swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	good_els--;
	filter_rfd_n++;
	continue;
      }
      if (gaps4[i]!=gaps_first) {
	// keep the most relevant one
	swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	good_els--;
	filter_gaps_n++;
	continue;
      }
      if (nrow[i]!=nrow_first) {
	// keep the most relevant one
	swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	good_els--;
	filter_nrow_n++;
	continue;
      }
      if (colstride[i]!=colstride_first) {
	// keep the most relevant one
	swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	good_els--;
	filter_colstride_n++;
	continue;
      }
      if (lastWordIdx[i]!=lastWordIdx_first) {
	// keep the most relevant one
	swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	good_els--;
	filter_lastWordIdx_n++;
	continue;
      }
      if (minsc[i]!=minsc_first) {
	// keep the most relevant one
	swap_data(i, good_els-1,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, lrmax, btnfilled);
	good_els--;
	filter_minsc_n++;
	continue;
      }
      i++;
   }

   fprintf(stderr, "INFO: Filtered from %i down to %i (%3.2f%%)\n",int(nels), int(good_els),100.0*good_els/nels);
   fprintf(stderr, "INFO: Filter pass2: %i iter: %i rfd: %i gaps: %i nrow: %i colstride: %i lastWordIdx: %i minsc: %i\n",
		   filter_pass2_n, filter_iter_n,filter_rfd_n,filter_gaps_n,filter_nrow_n,filter_colstride_n,filter_lastWordIdx_n,filter_minsc_n);
   return good_els;
}

#ifdef CPUVECT
// The input buffers are meant for scalar execution
// If we have vectors, we must transpose some buffers accordingly
// Remonder: We groupped same rfd together, so use that info
//
// rf
//  Org:
//   el=0 0...rfd0..MAX_RF_EL-1
//   el=1                       MAX_RF_EL..+rfd1..2*MAX_RF_EL-1
//   ...
//   el=N                                                            N*MAX_RF_EL..+rfdN..(N-1)*MAX_RF_EL-1
//  New:
//   el=0   0      V           ...   V*(rfd-1)
//   el=1    1      V+1          ...         V*(rfd-1)+1
//   ...
//   el=V-1     V-1     2*V-1            ...              V*rfd-1
//   el=V                                                         V*rfd ...
//
void vect_transpose_rf(const int nels,
		const size_t  rfd,
		const char rf_in[],
		char       rf_out[]
		) {
    // Note: We are assuming rounded nels
#pragma omp parallel for collapse(2)
    for (int b=0; b<(nels/VECT_SIZE); b++) {
       for (size_t r=0; r<rfd; r++) {
	  for (int v=0; v<VECT_SIZE; v++) {
             uint64_t source_i = (b*VECT_SIZE+v)*uint64_t(MAX_RF_EL)+r;
             uint64_t dest_i = (b*uint64_t(rfd)+r)*uint64_t(VECT_SIZE)+v;
	     rf_out[dest_i] = rf_in[source_i];
	  }
       }
    }
}

// The input buffers are meant for scalar execution
// If we have vectors, we must transpose some buffers accordingly
// Remonder: We groupped same iter together, so use that info
//
// rf
//  Org: (vs0,vs1) always together, iter consecutive, x5 
//   el=0 0 1 ...itr5a..MAX_PB_EL-1
//   el=1                       MAX_PB_EL.+1..+itr5b..2*MAX_PB_EL-1
//   ...
//   el=N                                                            N*MAX_PB_EL.+1..+itr5n..(N-1)*MAX_PB_EL-1
//  New: (vs0x5) (vs1x5), repeated iter times
//   el=0   0      V           ...   V*(itr5-1)
//   el=1    1      V+1          ...         V*(itr5-1)+1
//   ...
//   el=V-1     V-1     2*V-1            ...              V*itr5-1
//   el=V                                                         V*itr5 ...
//
void vect_transpose_profbuf(const int nels,
		const size_t  iter,
		const uint8_t profbuf_in[],
		uint8_t       profbuf_out[]
		) {
    // Note: We are assuming rounded nels
#pragma omp parallel for collapse(3)
    for (int b=0; b<(nels/VECT_SIZE); b++) {
      for (size_t j=0; j<iter; j++) {
        for (size_t r=0; r<5; r++) {
	  for (int v=0; v<VECT_SIZE; v++) {
	     // [MAX_PB_EL,5,iter,2],
             // ((b,v),r,j,0)
             uint64_t source_vs0 = (b*VECT_SIZE+v)*uint64_t(MAX_PB_EL)+((r*iter+j)*2);
             uint64_t source_vs1 = source_vs0+1;
	     // [BLK,iter,2,5,v]
	     // (b,j,0,r,v)
             uint64_t dest_vs0 = ((b*uint64_t(iter)+j)*2*5+r)*uint64_t(VECT_SIZE)+v;
             uint64_t dest_vs1 = dest_vs0+VECT_SIZE*5;
	     profbuf_out[dest_vs0] = profbuf_in[source_vs0];
	     profbuf_out[dest_vs1] = profbuf_in[source_vs1];
	  }
	}
      }
    }
}

#endif

#ifndef CPUVECT
int align_ee8_one(const int el, // for debuggging purpose
		const size_t nrow,
		const size_t iter,
		const size_t colstride,
		const size_t lastWordIdx,
		const int32_t minsc,
		const size_t rfd,
		const uint8_t gaps[4],
                const SSERegI *profbuf,
		const char    *rf,
                TMatBuf       *mat,
		DpBtCandidate *btncand,
                int32_t       &computed_lrmax, // out
                const int32_t ref_lrmax,
                const int32_t ref_btnfilled) {
#ifndef PRE_LR_SCALAR
	uint16_t btnfilled = 0;
#ifndef PRE_M11_SCALAR
	const EEU8_TCScore lrmax = EEU8_alignNucleotides<uint16_t>(profbuf, rf, rfd,
					mat,
                                        iter, colstride, lastWordIdx,
					minsc, nrow,
					btncand, btnfilled,
					gaps[0],gaps[1],gaps[2],gaps[3]);
#else
	const EEU8_TCScore lrmax = EEU8_alignNucleotidesM11Scalar<uint16_t,151>(profbuf, rf, rfd,
					mat,
					minsc, nrow,
					btncand, btnfilled,
					gaps[0],gaps[1],gaps[2],gaps[3]);

#endif

#else

#if 0
	// The upstream filter takes care of this...
	// Note: The vast majority of reads have iter==151
	if (iter!=151) {
		computed_lrmax = SKIPPED_LRMAX; // special value to signal we did not actually compute it
		return 0;  // TODO: We may properly use the right template function, or call slow align directly
	}
#endif

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
#if  !defined(PRE_LR_SCALAR)
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

#else
template<int VLen>
int align_ee8_vect(const int el, // first, for debuggging purpose
		const size_t nrow,
		const size_t iter,
		const size_t colstride,
		const size_t lastWordIdx,
		const int32_t minsc,
		const size_t rfd,
		const uint8_t gaps[4],
                const SSERegI *profbuf,
		const char    *rf,
                TMatBuf       *mat,
		DpBtCandidate *btncand,
                int32_t       *computed_lrmax, // out
                const int32_t *ref_lrmax,
                const int32_t *ref_btnfilled) {
	typedef uint8_t vectu8_t __attribute__((ext_vector_type(VLen)));
	uint16_t btnfilled[VLen];
        //fprintf(stderr,"Using profbuf %p\n",profbuf);
        //fprintf(stderr,"Using rf %p\n",rf);
	const vectu8_t lrmax = EEU8_alignNucleotidesVect<VLen,uint16_t>(
					(const vectu8_t*) profbuf, (const vectu8_t*) rf, rfd,
					mat,
                                        iter, colstride, lastWordIdx,
					minsc, nrow,
					btncand, btnfilled,
					gaps[0],gaps[1],gaps[2],gaps[3]);
	int nerrs = 0;
	for(int i=0; i<VLen; i++) computed_lrmax[i] = lrmax[i];
	for(int i=0; i<VLen; i++) if (int(ref_lrmax[i]) != int(lrmax[i])) nerrs++;
#ifndef PRE_LR_SCALAR
	for(int i=0; i<VLen; i++) if (int(ref_btnfilled[i]) != int(btnfilled[i])) nerrs++;
#else
	//TAlScore sc = (TAlScore)(lrmax - 0xff);
	//if ((ref_btnfilled!=0) != (sc>=minsc)) nerrs++;
#endif

#if 0

#ifndef NO_CHECK_PRINT
	if (int(ref_lrmax) != int(lrmax)) fprintf(stderr, "[%i] MISMATCH in lrmax (%i != %i)\n",el,int(lrmax), int(ref_lrmax));
#ifndef PRE_LR_SCALAR
	if (int(ref_btnfilled) != int(btnfilled)) fprintf(stderr, "[%i] MISMATCH in ref_btnfilled (%i != %i)\n",el,int(btnfilled), int(ref_btnfilled));
#else
	if ((ref_btnfilled!=0) != (sc>=minsc)) fprintf(stderr, "[%i] MISMATCH in ref_btnfilled (%i) vs sc %i minsc %i)\n",el,int(ref_btnfilled), int(sc), int(minsc));
#endif
#endif
#endif
	// TODO: In case of PRE_LR_SCALAR, we should also likely want to test the complete align on the elements that need it
	//       But that's probably better done in a separate test
	return nerrs;
}

#endif

void align_ee8(const int npar, const int nels,
		const size_t nrow,
		const size_t iter,
		const size_t colstride,
		const size_t lastWordIdx,
		const int32_t minsc,
		const size_t  rfd,
		const uint8_t gaps[4],
                const SSERegI profbuf[],
		const char    rf[],
                TMatBuf       mat[], 
		DpBtCandidate btncand[],
                int32_t       computed_lrmax[], // out
                const int32_t ref_lrmax[],
                const int32_t ref_btnfilled[]) {
   int elp_pp = (nels+ (npar-1))/npar; // round up
   // Keep block multiple of 64
   // That's cache line size for bytes
   elp_pp = ((elp_pp + (64-1))/64)*64; // round up
   fprintf(stderr, "INFO: Starting %i parallel iterations of %i each (last %i).\n",int(npar), int(elp_pp), std::max(int(nels-(npar-1)*elp_pp),int(0)));
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
	TMatBuf       *my_mat = nullptr;
	DpBtCandidate *my_btncand = nullptr;
#elif !defined(OMPGPU)
	// the following loop is sequential, so can reuse the buffer
	TMatBuf       *my_mat = mat+p*size_t(MAX_MAT_EL);
	DpBtCandidate *my_btncand = btncand+p*size_t(MAX_RF_EL);
#endif
#ifndef CPUVECT

#ifdef OMPGPU
#pragma omp parallel for reduction(+:nerrs)
#endif
      for (int i=p*elp_pp; i<iend; i++) {
#if defined(OMPGPU) && !defined(PRE_LR_SCALAR)
	// the loop is parallel, so we must use a different buffer for each element
#ifndef PRE_M11_SCALAR
	TMatBuf       *my_mat = mat+i*size_t(MAX_MAT_EL);
#else
	// output is aligned on the 64-wide vector
	int b64 = i/64;
	int v64 = i%64;
	TMatBuf       *my_mat = mat+b64*64*ROWSTRIDE*iter*rfd + v64;
#endif
	DpBtCandidate *my_btncand = btncand+i*size_t(MAX_RF_EL);
#endif
         nerrs+= align_ee8_one(i,
		nrow, iter, colstride, lastWordIdx,minsc, rfd, gaps,
                profbuf+i*size_t(MAX_PB_EL),rf+i*size_t(MAX_RF_EL),
		my_mat,my_btncand,
		computed_lrmax[i],ref_lrmax[i],ref_btnfilled[i]);
      }

#else /* CPUVECT */
      // Note: Assuming nels was rounded
      for (int i=p*elp_pp; (i+VECT_SIZE)<=iend; i+=VECT_SIZE) {
         nerrs+= align_ee8_vect<VECT_SIZE>(i,
		nrow, iter, colstride, lastWordIdx,minsc, rfd, gaps,
                profbuf+i*size_t(iter*10), rf+i*size_t(rfd), // transposed, now different spacing
		my_mat,my_btncand,
		computed_lrmax+i,ref_lrmax+i,ref_btnfilled+i);
      }
#endif
   }

   if (nerrs!=0) {
      fprintf(stderr, "FAILED matching results %i times!\n", nerrs);
   } else {
      fprintf(stderr, "SUCCESS, all results matched.\n");
   }

}

int load_and_align_ee8(const int npar, int nels, bool pass2_only, int skip_blocks) {
   int32_t *computed_lrmax = new int32_t[nels];
   int32_t *ref_lrmax = new int32_t[nels];
   int32_t *ref_btnfilled = new int32_t[nels];
   SSERegI *profbuf = new (std::align_val_t(1024))  SSERegI[size_t(nels)*MAX_PB_EL];
   char    *rf = new  (std::align_val_t(1024)) char[size_t(MAX_RF_EL)*nels];
   size_t  *nrow = new size_t[nels];
   size_t  *iter = new size_t[nels];
   size_t  *colstride = new size_t[nels];
   size_t  *lastWordIdx = new size_t[nels];
   int32_t  *minsc = new int32_t[nels];
   size_t  *rfd = new size_t[nels];
   uint8_t *gaps = new uint8_t[4*nels];

   auto t1 = std::chrono::high_resolution_clock::now();
   // test tool, don't worry about perfect cleanup
   if (load_data(nels, skip_blocks,
		nrow, iter, colstride, lastWordIdx,minsc,rfd,
		profbuf,rf,gaps) !=0 ) return 1;
   if (load_results(nels, skip_blocks,
		ref_lrmax, ref_btnfilled) !=0 ) return 1;
   nels = filter(nels,
                nrow, iter, colstride, lastWordIdx, minsc,
                rfd, profbuf, rf, gaps, ref_lrmax, ref_btnfilled,
		pass2_only);
#ifdef CPUVECT
   // TODO: Fix the code to support non-multiple numbers
   if ((nels%VECT_SIZE)!=0) {
      nels = (nels/VECT_SIZE)*VECT_SIZE; // round down
      fprintf(stderr, "INFO: Rounded down to %i\n",int(nels));
   }
   // create support buffers
   char    *rf_alt = new  (std::align_val_t(1024)) char[size_t(rfd[0])*nels];
   uint8_t *profbuf_alt = new (std::align_val_t(1024)) uint8_t[size_t(iter[0])*10*nels];
   // need transposed values for vector operation
   vect_transpose_rf(nels,rfd[0],rf,rf_alt); // rfd homogeneous
   vect_transpose_profbuf(nels,iter[0],profbuf,profbuf_alt); // iter homogeneous
   std::swap(rf,rf_alt); // use the transposed version from now on... keeping the old just for debugging
   std::swap(profbuf,profbuf_alt); // use the transposed version from now on... keeping the old just for debugging
#endif
   auto t2a = std::chrono::high_resolution_clock::now();
   for (int i=0; i<nels; i++) computed_lrmax[i] = -1; // invalidate
#if defined(PRE_LR_SCALAR)
   // no mat or btncand, just pass NULLs around 
   TMatBuf        *mat = nullptr;
   DpBtCandidate *btncand = nullptr;
#elif !defined(OMPGPU)
   // one per thread is enough on the CPU
   TMatBuf       *mat = new  (std::align_val_t(1024))  TMatBuf[size_t(npar)*MAX_MAT_EL];
   DpBtCandidate *btncand = new  (std::align_val_t(1024)) DpBtCandidate[size_t(MAX_RF_EL)*npar];
#else
   // no shortcuts on the GPU
   TMatBuf       *mat = new  (std::align_val_t(1024)) TMatBuf[size_t(nels)*MAX_MAT_EL];
   DpBtCandidate *btncand = new  (std::align_val_t(1024)) DpBtCandidate[size_t(MAX_RF_EL)*nels];
#endif
   auto t2b = std::chrono::high_resolution_clock::now();
   // we filtered the eleemnts, we can now assume they are homogeneous
   fprintf(stderr, "INFO: #iter = %i #rfd = %i #row = %i #minsc = %i\n",int(iter[0]),int(rfd[0]),int(nrow[0]),int(minsc[0]));
   align_ee8(npar, nels,
		nrow[0], iter[0], colstride[0], lastWordIdx[0],minsc[0],rfd[0], // homogeneous
		gaps, // 4 elements, so passing by pointer, but still homogeneous
		profbuf,rf,
		mat,btncand,
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
		nrow[0], iter[0], colstride[0], lastWordIdx[0],minsc[0],rfd[0], // homogeneous
		gaps, // 4 elements, so passing by pointer, but still homogeneous
		profbuf,rf,
		mat,btncand,
		computed_lrmax,ref_lrmax,ref_btnfilled);
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

#ifdef CPUVECT
   delete[] profbuf_alt;
   delete[] rf_alt;
#endif
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
   if (argv<3) {
     fprintf(stderr, "Usage: test_align_ee8 <npar> <nels> [<pass2_only>] [<skip_blocks>]\n");
     return 1;
   }
   int npar = atoi(argc[1]);
   int nels = atoi(argc[2]);
   int pass2_only = 0;
   int skip_blocks = 0;
   if (argv>3) {
     pass2_only = atoi(argc[3]);
   }
   if (argv>4) {
     skip_blocks = atoi(argc[4]);
   }

#ifdef USE_BUILTIN_SUB_SAT
   {
     // make sure this compiler is not broken
     uint8_t a = 20;
     uint8_t b = 40;
     if (__builtin_elementwise_sub_sat(a,b)!=0) {
        fprintf(stderr, "ERROR: Broken compiler detected, builtin scalar sub_sat produces wrong result!\n");
	return 2;
     }
   }
   {
     // make sure this compiler is not broken
     typedef uint8_t uint8x16_t __attribute__((ext_vector_type(16)));
     uint8x16_t a = 20;
     uint8x16_t b = 40;
     uint8x16_t r = 0;
     if (__builtin_reduce_or(__builtin_elementwise_sub_sat(a,b)!=r)) { // enough is one is not equal
        fprintf(stderr, "ERROR: Broken compiler detected, builtin vector sub_sat produces wrong result!\n");
	return 2;
     }
   }
#endif
   return load_and_align_ee8(npar,nels,pass2_only!=0i,skip_blocks);
}
