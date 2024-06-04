/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <array>
#include <execution>
#include "aligner_cache.h"
#include "aligner_seed.h"
#include "search_globals.h"
#include "bt2_idx.h"

using namespace std;

/**
 * Construct a constraint with no edits of any kind allowed.
 */
Constraint Constraint::exact() {
	Constraint c;
	c.edits = c.mms = c.ins = c.dels = c.penalty = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::penaltyBased(int pen) {
	Constraint c;
	c.penalty = pen;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint related to the length of the read.
 */
Constraint Constraint::penaltyFuncBased(const SimpleFunc& f) {
	Constraint c;
	c.penFunc = f;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::mmBased(int mms) {
	Constraint c;
	c.mms = mms;
	c.edits = c.dels = c.ins = 0;
	return c;
}

/**
 * Construct a constraint where the only constraint is a total
 * penalty constraint.
 */
Constraint Constraint::editBased(int edits) {
	Constraint c;
	c.edits = edits;
	c.dels = c.ins = c.mms = 0;
	return c;
}

// Input to seachSeedBi
class SeedAlignerSearchParams {
public:
	class CacheAndSeed {
	public:
		CacheAndSeed()
		: seq(NULL), seq_len(0), n_seed_steps(0)  {}

		CacheAndSeed(
			const char *   _seq,             // sequence of the local seed alignment cache
			const uint32_t _seq_len,         // and its length
			const InstantiatedSeed& _seed    // current instantiated seed
		)
		: seq(NULL), seq_len(0), n_seed_steps(0)  // just set a default
		{ reset(_seq, _seq_len, _seed); }

		void reset(
			const char *   _seq,             // sequence of the local seed alignment cache
			const uint32_t _seq_len,         // and its length
			const InstantiatedSeed& _seed    // current instantiated seed
		)
		{
			seq = _seq;
			if (_seq_len>127) {printf("Unexpected seq_len %i\n",int(seq_len)); throw 1;}
			seq_len = _seq_len;
			if (abs(_seed.n_steps)>127) {printf("Unexpected n_seed_steps %i\n",int(_seed.n_steps)); throw 1;}
			n_seed_steps = _seed.n_steps;
		}

		CacheAndSeed(CacheAndSeed &other) = default;
		CacheAndSeed(CacheAndSeed &&other) = default;

		CacheAndSeed& operator=(const CacheAndSeed& other) = default;

		// Since we only suport exact matches, use -n_steps
		int seed_step_min() const {return -n_seed_steps;}

		// Maximum number of positions that the aligner may advance before
		// its first step.  This lets the aligner know whether it can use
		// the ftab or not.
		int maxjump() const {return n_seed_steps;}
	
		// Use pointers, so they can be changed 
		const char *seq;
		uint8_t seq_len;
		int8_t n_seed_steps;               // steps in the current instantiated seed
	};

	// create an empty bwt
	// and constratins from seed, for initial searchSeedBi invocation
	SeedAlignerSearchParams(
		const char *   _seq,             // sequence of the local seed alignment cache
		const uint32_t _seq_len,         // and its length
		const InstantiatedSeed& _seed    // current instantiated seed
	)
	: cs(_seq,_seq_len,_seed)
	{}

	// create an empty bwt
	// and constratins from seed, for initial searchSeedBi invocation
	SeedAlignerSearchParams()
	: cs()
	{}

	SeedAlignerSearchParams& operator=(const SeedAlignerSearchParams& other) = default;

	void reset(
		const char *   _seq,             // sequence of the local seed alignment cache
		const uint32_t _seq_len,         // and its length
		const InstantiatedSeed& _seed    // current instantiated seed
	)
	{
	  cs.reset(_seq,_seq_len,_seed);
	}

	// A sub-class for historical reasons
	CacheAndSeed cs;      // local seed alignment cache and associated instatiated seed
};


class SeedAlignerSearchData {
public:
	// create an empty bwt
	// and constratins from seed, for initial searchSeedBi invocation
	SeedAlignerSearchData()
	: bwt()
	, sak()
	, need_reporting(false)
	{}

	SeedAlignerSearchData& operator=(const SeedAlignerSearchData& other) = default;

	void resetData(
		const char *   seq,            // sequence of the local seed alignment cache
		const uint8_t seq_len          // and its length
	) {
	  bwt.set(0,0);
	  sak.init(seq, seq_len);
	  need_reporting = false;
	}

	void set_reporting() { need_reporting = true; }

	BwtTopBotFw bwt;      // The 2 BWT idxs
	SAKey       sak;      // seed key
	bool need_reporting;
};

class SeedAlignerSearchState {
public:

	SeedAlignerSearchState()
	: tloc()
	, bloc()
	, step(0)
	{}

	void reset() {
		tloc.invalidate();
		bloc.invalidate();
		step = 0;
	}

	SideLocus tloc;       // locus for top (perhaps unititialized)
	SideLocus bloc;       // locus for bot (perhaps unititialized)
	int step;             // depth into steps[] array, i.e. step_min+step
};



//
// Some static methods for constructing some standard SeedPolicies
//

/**
 * Given a read, depth and orientation, extract a seed data structure
 * from the read and fill in the steps & zones arrays.  The Seed
 * contains the sequence and quality values.
 */

bool
InstantiatedSeed::instantiateExact(
	const int seed_len,
	const char *seq) // seed read sequence
{
	int seedlen = seed_len;
	assert_gt(seedlen, 0);
	n_steps = seedlen;
	//step_min = -seedlen;
	
	// Take a sweep through the seed sequence.  Consider where the Ns
	// occur and how zones are laid out.  Calculate the maximum number
	// of positions we can jump over initially (e.g. with the ftab) and
	// perhaps set this function's return value to false, indicating
	// that the arrangements of Ns prevents the seed from aligning.
	for(int i = 0; i < seedlen; i++) {
		const int stepi = i-seedlen; //step_min+i;
		int off = abs(stepi)-1;
		int c = seq[off];  assert_range(0, 4, c);
		if(c == 4) {
			// cons.canN always false since mms and edit both ==0
			// Seed disqualified due to arrangement of Ns
			invalidate();
			return false;
		}
	}
	return true;
}

/**
 * Return a set consisting of 1 seed encapsulating an exact matching
 * strategy.
 */
void
Seed::zeroMmSeeds(int ln, EList<Seed>& pols) {
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_EXACT;
}


void
Seed::zeroMmSeed(int ln, Seed& pol) {
	pol.len = ln;
	pol.type = SEED_TYPE_EXACT;
}

/**
 * Types of actions that can be taken by the SeedAligner.
 */
enum {
	SA_ACTION_TYPE_RESET = 1,
	SA_ACTION_TYPE_SEARCH_SEED, // 2
	SA_ACTION_TYPE_FTAB,        // 3
	SA_ACTION_TYPE_FCHR,        // 4
	SA_ACTION_TYPE_MATCH,       // 5
	SA_ACTION_TYPE_EDIT         // 6
};

/**
 * Given a read and a few coordinates that describe a substring of the read (or
 * its reverse complement), fill in 'seq' and 'qual' objects with the seed
 * sequence and qualities.
 *
 * The seq field is filled with the sequence as it would align to the Watson
 * reference strand.  I.e. if fw is false, then the sequence that appears in
 * 'seq' is the reverse complement of the raw read substring.
 */
void
SeedAligner::instantiateSeq(
	const Read& read, // input read
	char        seq[], // output sequence
	int len,          // seed length
	int depth,        // seed's 0-based offset from 5' end
	bool fw)         // seed's orientation
{
	// Fill in 'seq'
	// If fw is false, we take characters starting at the 3' end of the
	// reverse complement of the read.
	for(int i = 0; i < len; i++) {
		seq[i] = read.patFw.windowGetDna(i, fw, depth, len);
	}
}

/**
 * We assume that all seeds are the same length.
 *
 * For each seed, instantiate the seed, retracting if necessary.
 */
void SeedAligner::prepareSeed(
	const int seed_len,        // search seed length
	size_t off,                // offset into read to start extracting
	int per,                   // interval between seeds
	const Read& read,          // read to align
	SeedResults& sr)           // holds all the seed hits
{
	assert_gt(read.length(), 0);
	// Check whether read has too many Ns
	const int len = seed_len;
	// Calc # seeds within read interval
	int nseeds = 1;
	if((int)read.length() - (int)off > len) {
		nseeds += ((int)read.length() - (int)off - len) / per;
	}

	const int min_len = std::min<int>(len, (int)read.length());
	sr.prepare(off, per, nseeds, min_len);
}

size_t SeedAligner::instantiateSeed(
	char*              seqBuf,      // content of all the seqs, no separators
	InstantiatedSeed*  seedsBuf,    // all the instantiated seeds, both fw and rc
	const Read& read,          // read to align
	bool nofw,                 // don't align forward read
	bool norc,                 // don't align revcomp read
	SeedResults& sr)           // holds all the seed hits
{
	int insts[3];
	insts[0] = insts[1] = insts[2] = 0;

	const int nseeds = sr.numOffs();
	const int min_len = sr.seqs_len();
	sr.reset(seqBuf, seedsBuf);

	const int off = sr.getOffset();
	const int per = sr.getInterval();
	// For each seed position
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			// Skip this orientation b/c user specified --nofw or --norc
			continue;
		}
		// For each seed position
		for(int i = 0; i < nseeds; i++) {
			int depth = i * per + off;
			// Extract the seed sequence at this offset
			// If fw == true, we extract the characters from i*per to
			// i*(per-1) (exclusive).  If fw == false, 
			instantiateSeq(
				read,
				sr.seqs(fw,i),
				min_len,
				depth,
				fw);
			// For each search strategy
			InstantiatedSeed& is = sr.instantiatedSeed(fw, i);
			{
				if(is.instantiateExact(
					min_len,
					sr.seqs(fw,i)))
				{
					// Can we fill this seed hit in from the cache?
					insts[0]++;
					insts[fw ? 1 : 2]++;
				}
			}
		}
	}

	return insts[0];
}

void MultiSeedResults::prepareOneSeed(
		uint32_t idx,              // index in the multi-seed array
	        size_t off,                // offset into read to start extracting
	        int per,                   // interval between seeds
	        const Read* pread)         // read to align
{
	SeedAligner::prepareSeed(
			_seed_len,    // length of a multiseed seed
			off,          // offset to begin extracting
			per,          // interval between seeds
			*pread,       // read to align
			_srs[idx]);
}

size_t MultiSeedResults::instantiateSeeds(Read const * const preads[])
{
	size_t seq_total_size = 0;
	size_t seeds_total_size = 0;

	for (uint32_t i=0; i<_n_sr; i++) {
		_bufOffs[i].first  = seq_total_size;
		_bufOffs[i].second = seeds_total_size;
		seq_total_size+=_srs[i].getSeqSize();
		seeds_total_size+=_srs[i].getSeedsSize();
	}

	if (_seqBuf_size<seq_total_size) {
		// need more space
		if (_seqBuf!=NULL) delete[] _seqBuf;
		_seqBuf = new char[seq_total_size];
		_seqBuf_size = seq_total_size;
	}

	if (_seedsBuf_size<seeds_total_size) {
		// need more space
		if (_seedsBuf!=NULL) delete[] _seedsBuf;
		_seedsBuf = new InstantiatedSeed[seeds_total_size];
		_seedsBuf_size = seeds_total_size;
	}

	size_t tries = 0;
#pragma omp parallel for reduction(+:tries) default(shared) schedule(dynamic,8)
	for (uint32_t i=0; i<_n_sr; i++) {
		tries += SeedAligner::instantiateSeed(
			_seqBuf + _bufOffs[i].first, _seedsBuf + _bufOffs[i].second,
			*preads[i], _nofw, _norc,
			_srs[i]);
	}

	return tries;
}

uint32_t SeedAligner::computeValidInstantiatedSeeds(
				const SeedResults& sr)
 {
	uint32_t seedsearches = 0;
	for(int fwi = 0; fwi < 2; fwi++) {
		const bool fw = (fwi == 0);
		// For each instantiated seed
		// start aligning and find list of seeds to search
		for(size_t i =0; i < sr.numOffs(); i++) {
			const InstantiatedSeed& is = sr.instantiatedSeed(fw, i);
			if(is.isValid()) seedsearches++;
			// else, Cache hit in an across-read cache
		}
	}
	seedsearches_ = seedsearches;
	return seedsearches;
}

/**
 * We assume that all seeds are the same length.
 *
 * For each seed:
 *
 * 1. Instantiate all seeds, retracting them if necessary.
 * 2. Calculate zone boundaries for each seed
 *
 * Return number of batches
 */
uint32_t SeedAligner::searchAllSeedsPrepare(
	const SeedResults& sr,
	const int ftabLen)           // forward index (BWT) value
{
	uint32_t seedsearches = 0;

	SeedAlignerSearchParams* paramVec = paramVec_;
	const auto seqs_len = sr.seqs_len();

	// Build the support structures
	// the order is arbirtrary
	for(int fwi = 0; fwi < 2; fwi++) {
		const bool fw = (fwi == 0);
		// For each instantiated seed
		// start aligning and find list of seeds to search
		for(size_t i =0; i < sr.numOffs(); i++) {
			const InstantiatedSeed& is = sr.instantiatedSeed(fw, i);
			if(!is.isValid()) {
				// Cache hit in an across-read cache
				continue;
			}
			const char *   seq = sr.seqs(fw,i);
			{
					assert_eq(fw, is.fw);
					assert_eq(i, (int)is.seedoffidx);
					paramVec[seedsearches].reset(seq, seqs_len, is);
			}
			seedsearches++;
		} // for i
	} // for fwi
	bwops_ = 0;
	assert_eq(seedsearches_,seedsearches);

	// return number of batches (rounded up)
	return (seedsearches+(ibatch_size-1))/ibatch_size;
}

inline void SeedAligner::searchAllSeedsDoAll(const Ebwt* ebwtFw)
{
	assert(ebwtFw != NULL);
	assert(ebwtFw->isInMemory());

	const SeedAlignerSearchParams* paramVec = paramVec_;
	SeedAlignerSearchData*         dataVec = dataVec_;
	const size_t seedsearches = seedsearches_;

	// do the searches in batches
	for (size_t mnr=0; mnr<seedsearches; mnr+=ibatch_size) {
		const size_t ibatch_max = std::min(mnr+ibatch_size,seedsearches);
		searchSeedBi<ibatch_size>(ebwtFw, bwops_, ibatch_max-mnr, &(paramVec[mnr]), &(dataVec[mnr]) );
	} // mnr loop

}

inline void SeedAligner::searchAllSeedsDoBatch(uint32_t ibatch, const Ebwt* ebwtFw)
{
	const size_t mnr = ibatch*ibatch_size;
	const size_t seedsearches = seedsearches_;
	const size_t ibatch_max = std::min(mnr+ibatch_size,seedsearches);
	if (mnr<ibatch_max) {
		assert(ebwtFw != NULL);
		assert(ebwtFw->isInMemory());

		// Note: Updates on bwops_ may not be atomi
		// But a small discrepancy is acceptable, as it is only rarely used diagnostics
		searchSeedBi<ibatch_size>(ebwtFw, bwops_, ibatch_max-mnr, &(paramVec_[mnr]), &(dataVec_[mnr]) );
	}
}


void SeedAligner::searchAllSeedsFinalize(
	AlignmentCache& cache,  // local cache for seed alignments
	SeedResults& sr)
{
	const SeedAlignerSearchData* dataVec = dataVec_;

	uint32_t seedsearches = 0;
	cache.clear();

	// finish aligning and add to SeedResult
	for(int fwi = 0; fwi < 2; fwi++) {
		const bool fw = (fwi == 0);
		// For each instantiated seed
		// start aligning and find list of seeds to search
		for(size_t i =0; i < sr.numOffs(); i++) {
			const InstantiatedSeed& is = sr.instantiatedSeed(fw, i);
			if(!is.isValid()) {
				// Cache hit in an across-read cache
				continue;
			}
			const SeedAlignerSearchData& sdata= dataVec[seedsearches];
			seedsearches++;

			if ( sdata.need_reporting ) {
				const SAKey &sak = sdata.sak;
				QVal  qv;
				cache.addOnTheFly(qv,sak, sdata.bwt.topf, sdata.bwt.botf);
				sr.add(qv,i,fw);
			}
		} // for i
	} // for fwi
	cache.finalize();
}

MultiSeedAligner::MultiSeedAligner(
		const Ebwt*       ebwtFw, // forward index (BWT)
		MultiSeedResults& srs)
 	: _ebwtFw(ebwtFw)
	, _srs(srs)
	, _als(new SeedAligner[srs.nSRs()])
	, _ftabLen(ebwtFw->eh().ftabChars()) // cache the value
	,  _paramVec(NULL), _dataVec(NULL)
	, _bufVec_size(0), _bufVec_filled(0)
{}

MultiSeedAligner::~MultiSeedAligner() {
	if (_dataVec!=NULL) delete[] _dataVec;
	if (_paramVec!=NULL) delete[] _paramVec;
	delete[] _als;
}

// Update buffer, based on content of _srs
void MultiSeedAligner::reserveBuffers()
{
	const uint32_t n_sr = _srs.nSRs();
	size_t buf_total_size = 0;

	// need to determine how much I actually need
	// also compute total size first to make sure buffers are big enough
#pragma omp parallel for reduction(+:buf_total_size) default(shared) schedule(static,8)
	for (uint32_t i=0; i<n_sr; i++) {
		// both update and use the size
		buf_total_size+=_als[i].computeValidInstantiatedSeeds(_srs.getSR(i));
	}

	if (_bufVec_size<buf_total_size) {
		// need bigger buffers
		if (_dataVec!=NULL) delete[] _dataVec;
		if (_paramVec!=NULL) delete[] _paramVec;
		_paramVec = new SeedAlignerSearchParams[buf_total_size];
		_dataVec = new SeedAlignerSearchData[buf_total_size];
		_bufVec_size = buf_total_size;
	}

	// now update buffers in SearchResults
	buf_total_size = 0;
	for (uint32_t i=0; i<n_sr; i++) {
		_als[i].setBufs(
			_paramVec+buf_total_size,
			_dataVec+buf_total_size); 
		buf_total_size+=_als[i].getBufsSize();
	}
	_bufVec_filled = buf_total_size;
	assert_leq(_bufVec_filled,_bufVec_size);
}

void MultiSeedAligner::searchAllSeedsDoAll()
{
	const Ebwt* ebwtFw= _ebwtFw;
	const SeedAlignerSearchParams* paramVec = _paramVec;
	SeedAlignerSearchData*         dataVec  = _dataVec;

	// do the searches in batches
	const uint64_t total_els  = _bufVec_filled;
	const uint32_t total_batches = (total_els+(ibatch_size-1))/ibatch_size; // round up
#ifdef FORCE_ALL_OMP
#pragma omp parallel for
	for (uint32_t gbatch=0; gbatch<total_batches; gbatch++) {
#else
	std::for_each_n(std::execution::par_unseq,
		thrust::counting_iterator(0), total_batches,
		[ebwtFw,paramVec,dataVec,total_els](uint32_t gbatch) mutable {
#endif
		const size_t start_el = gbatch*ibatch_size;
		size_t end_el = start_el+ibatch_size;
		if (end_el>total_els) end_el = total_els;
		uint64_t bwops; // just ignore the bwops for now, keep it local
		SeedAligner::searchSeedBi<ibatch_size>(ebwtFw, bwops, end_el-start_el, &(paramVec[start_el]), &(dataVec[start_el]));
	} // for gbatch
#ifndef FORCE_ALL_OMP
	); // for_each
#endif

}


/**
 * Get tloc, bloc ready for the next step.  If the new range is under
 * the ceiling.
 */
inline void nextLocsBi(
	const EbwtParams& ep,         // index params
	const uint8_t* ebwt,          // index data
	const BwtTopBotFw &bwt,
	SideLocus& tloc,              // top locus
	SideLocus& bloc)              // bot locus
{
	assert_gt(botf, 0);
	assert(botb > 0);
	{
		if(bwt.botf - bwt.topf == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(bwt.topf, ep, ebwt);
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				bwt.topf, bwt.botf, ep, ebwt, tloc, bloc);
			assert(bloc.valid());
		}
	}
	assert(bwt.botf - bwt.topf == 1 ||  bloc.valid());
	assert(bwt.botf - bwt.topf > 1  || !bloc.valid());
}

// return true, if we are already done
inline bool startSearchSeedBi(
	const EbwtParams& ep,         // index params
	const uint8_t* ebwt,          // index data
	const TIndexOffU *ftab,
	const TIndexOffU *eftab,
	const TIndexOffU fchr[5],     // index fchr
	const SeedAlignerSearchParams &p,
	SeedAlignerSearchState &sstate,
	BwtTopBotFw &bwt)
{

	assert_eq(sstate.step, 0);
	assert_gt(p.cs.n_seed_steps, 0);
	{
		const char *seq = p.cs.seq;
		// Just starting
		assert(!sstate.tloc.valid());
		assert(!sstate.bloc.valid());
		const int seed_step_min = p.cs.seed_step_min();
		int off = abs(seed_step_min)-1;
		// Check whether/how far we can jump using ftab or fchr
		int ftabLen = ep.ftabChars();
		if (ftabLen > 1 && ftabLen <= p.cs.maxjump()) {
			TIndexOffU fwi0 = Ebwt::ftabSeqToInt(ftabLen, true, seq, off - ftabLen + 1, false);
			Ebwt::ftabLoHi(ftab, eftab, ep,
					fwi0,
					bwt.topf, bwt.botf);
			if(bwt.botf - bwt.topf == 0) return true;
			sstate.step += ftabLen;
		} else if(p.cs.maxjump() > 0) {
			// Use fchr
			const int c = seq[off];
			assert_range(0, 3, c);
			bwt.topf = fchr[c];
			bwt.botf = fchr[c+1];
			if(bwt.botf - bwt.topf == 0) return true;
			sstate.step++;
		} else {
			assert_eq(0, p.cs.maxjump());
			bwt.topf = 0;
			bwt.botf = fchr[4];
		}
		if(sstate.step == p.cs.n_seed_steps) {
			return true;
		}
		nextLocsBi(ep, ebwt, bwt, sstate.tloc, sstate.bloc);
		assert(sstate.tloc.valid());
	} 
	assert(sstate.tloc.valid());
	assert(bwt.botf - bwt.topf == 1 ||  sstate.bloc.valid());
	assert(bwt.botf - bwt.topf > 1  || !sstate.bloc.valid());
	assert_geq(sstate.step, 0);

	return false;
}

/**
 * Given a seed, search.  Assumes zone 0 = no backtracking.
 *
 * Return a list of Seed hits.
 * 1. Edits
 * 2. Bidirectional BWT range(s) on either end
 */
template<uint8_t SS_SIZE>
void
SeedAligner::searchSeedBi(
                        const Ebwt* ebwt,       // forward index (BWT)
                        uint64_t& bwops_,         // Burrows-Wheeler operations
                        const uint8_t nparams,
			const SeedAlignerSearchParams paramVec[],
			SeedAlignerSearchData dataVec[])
{
	const EbwtParams& ep = ebwt->eh();
	const uint8_t* const ebwtPtr = ebwt->ebwt();
	const TIndexOffU * const ftab = ebwt->ftab();
	const TIndexOffU * const eftab = ebwt->eftab();
	const TIndexOffU * const fchr =  ebwt->fchr();

	uint8_t idxs[SS_SIZE]; // indexes into sstateVec and paramVec

	SeedAlignerSearchState sstateVec[SS_SIZE]; // work area
	assert(nparams<=SS_SIZE);

	uint8_t nleft = nparams; // will keep track of how many are not done yet

	{
	   uint8_t n=0;
           uint8_t iparam = 0; // iparam and n may diverge, if some are done at init stage
	   while (n<nleft) {
		const SeedAlignerSearchParams& p = paramVec[iparam];
		SeedAlignerSearchData&   sdata   = dataVec[iparam];
		SeedAlignerSearchState&  sstate  = sstateVec[iparam];
		sdata.resetData(p.cs.seq, p.cs.seq_len);
		sstate.reset();
		idxs[n] = iparam;
		iparam+=1;
		const bool done = startSearchSeedBi(
					ep, ebwtPtr, ftab, eftab, fchr,
					p, sstate, sdata.bwt);
		if(done) {
		        if(sstate.step == (int)p.cs.n_seed_steps) {
                		// Finished aligning seed
				sdata.set_reporting();
			}
			// done with this, swap with last and reduce nleft
			nleft-=1;
		} else {
			n+=1;
		}
	    }
	}


	while (nleft>0) {
	   // Note: We can do the params in any order we want
	   // but we must do the steps inside the same param in order
	   // We still want to do them sequentially, not in parallel, 
	   // to give time for the prefetch to do its job.
	   uint32_t bwops = 0;
	   uint8_t n=0;
           // logically a for (uint32_t n=0; n<nleft; n++) but with nleft potentially changing
	   while (n<nleft) {
		const uint8_t iparam = idxs[n];

		const SeedAlignerSearchParams& p = paramVec[iparam];
		SeedAlignerSearchData&   sdata   = dataVec[iparam];
		SeedAlignerSearchState&  sstate  = sstateVec[iparam];

		const int n_seed_steps = p.cs.n_seed_steps;
		if (sstate.step >= (int) n_seed_steps) {
			// done with this, swap with last and reduce nleft
			nleft-=1;
			if (n<nleft) idxs[n] = idxs[nleft];
			continue;
		}

		assert_gt(sdata.bwt.botf, sdata.bwt.topf);
		assert(sdata.bwt.botf - sdata.bwt.topf == 1 ||  sstate.bloc.valid());
		assert(sdata.bwt.botf - sdata.bwt.topf > 1  || !sstate.bloc.valid());
		assert(sstate.tloc.valid());

		const int seed_step_min = p.cs.seed_step_min();
		const char *seq = p.cs.seq;

		SeedAlignerSearchWorkState wstate(seed_step_min+sstate.step);
		sstate.step++; // get ready for the next iteration

		if(sstate.bloc.valid()) {
			// Range delimited by tloc/bloc has size >1.  If size == 1,
			// we use a simpler query (see if(!bloc.valid()) blocks below)
			bwops++;
			ebwt->mapBiLFEx(sstate.tloc, sstate.bloc, wstate.t, wstate.b);
		}
		int c = seq[wstate.off];  assert_range(0, 4, c);
		//
		if(c == 4) { // couldn't handle the N
			// done with this, swap with last and reduce nleft
			nleft-=1;
			if (n<nleft) idxs[n] = idxs[nleft];
			continue;
		}
		if(!sstate.bloc.valid()) {
			assert(wstate.bp[c] == wstate.tp[c]+1);
			// Range delimited by tloc/bloc has size 1
			bwops++;
			const TIndexOffU ntop = sdata.bwt.topf;
			wstate.t[c] = ebwt->mapLF1(ntop, sstate.tloc, c);
			if(wstate.t[c] == OFF_MASK) {
				// done with this, swap with last and reduce nleft
				nleft-=1;
				if (n<nleft) idxs[n] = idxs[nleft];
				continue;
			}
			assert_geq(wstate.t[c], fchr[c]);
			assert_lt(wstate.t[c],  fchr[c+1]);
			wstate.b[c] = wstate.t[c]+1;
			assert_gt(wstate.b[c], 0);
		}
		if(wstate.b[c] == wstate.t[c]) {
			// done with this, swap with last and reduce nleft
			nleft-=1;
			if (n<nleft) idxs[n] = idxs[nleft];
			continue;
		}
		sdata.bwt.set(wstate.t[c], wstate.b[c]);
		if(sstate.step == n_seed_steps) {
			sdata.set_reporting();
			// done with this, swap with last and reduce nleft
			nleft-=1;
			if (n<nleft) idxs[n] = idxs[nleft];
			continue;
		}
		nextLocsBi(ep, ebwtPtr, sdata.bwt, sstate.tloc, sstate.bloc);
		// not done, move to the next element
		n+=1;
	   } // while n
	   bwops_ += bwops;
	} // while nleft

	return;
}

#ifdef ALIGNER_SEED_MAIN

#include <getopt.h>
#include <string>

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(const char *errmsg, const char *arg) {
	long l;
	char *endPtr = NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	throw 1;
	return -1;
}

enum {
	ARG_NOFW = 256,
	ARG_NORC,
	ARG_MM,
	ARG_SHMEM,
	ARG_TESTS,
	ARG_RANDOM_TESTS,
	ARG_SEED
};

static const char *short_opts = "vCt";
static struct option long_opts[] = {
	{(char*)"verbose",  no_argument,       0, 'v'},
	{(char*)"timing",   no_argument,       0, 't'},
	{(char*)"nofw",     no_argument,       0, ARG_NOFW},
	{(char*)"norc",     no_argument,       0, ARG_NORC},
	{(char*)"mm",       no_argument,       0, ARG_MM},
	{(char*)"shmem",    no_argument,       0, ARG_SHMEM},
	{(char*)"tests",    no_argument,       0, ARG_TESTS},
	{(char*)"random",   required_argument, 0, ARG_RANDOM_TESTS},
	{(char*)"seed",     required_argument, 0, ARG_SEED},
};

static void printUsage(ostream& os) {
	os << "Usage: ac [options]* <index> <patterns>" << endl;
	os << "Options:" << endl;
	os << "  --mm                memory-mapped mode" << endl;
	os << "  --shmem             shared memory mode" << endl;
	os << "  --nofw              don't align forward-oriented read" << endl;
	os << "  --norc              don't align reverse-complemented read" << endl;
	os << "  -t/--timing         show timing information" << endl;
	os << "  -v/--verbose        talkative mode" << endl;
}

bool gNorc = false;
bool gNofw = false;
int gVerbose = 0;
int gGapBarrier = 1;
int gSnpPhred = 30;
bool gReportOverhangs = true;

extern void aligner_seed_tests();
extern void aligner_random_seed_tests(
	int num_tests,
	TIndexOffU qslo,
	TIndexOffU qshi,
	uint32_t seed);

/**
 * A way of feeding simply tests to the seed alignment infrastructure.
 */
int main(int argc, char **argv) {
	bool useMm = false;
	bool useShmem = false;
	bool mmSweep = false;
	bool noRefNames = false;
	bool sanity = false;
	bool timing = false;
	int option_index = 0;
	int seed = 777;
	int next_option;
	do {
		next_option = getopt_long(
			argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 'v':       gVerbose = true; break;
			case 't':       timing   = true; break;
			case ARG_NOFW:  gNofw    = true; break;
			case ARG_NORC:  gNorc    = true; break;
			case ARG_MM:    useMm    = true; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_SEED:  seed = parseInt("", optarg); break;
			case ARG_TESTS: {
				aligner_seed_tests();
				aligner_random_seed_tests(
					100,     // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					18);   // pseudo-random seed
				return 0;
			}
			case ARG_RANDOM_TESTS: {
				seed = parseInt("", optarg);
				aligner_random_seed_tests(
					100,   // num references
					100,   // queries per reference lo
					400,   // queries per reference hi
					seed); // pseudo-random seed
				return 0;
			}
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
			}
		}
	} while(next_option != -1);
	char *reffn;
	if(optind >= argc) {
		cerr << "No reference; quitting..." << endl;
		return 1;
	}
	reffn = argv[optind++];
	if(optind >= argc) {
		cerr << "No reads; quitting..." << endl;
		return 1;
	}
	string ebwtBase(reffn);
	BitPairReference ref(
		ebwtBase,    // base path
		false,       // whether we expect it to be colorspace
		sanity,      // whether to sanity-check reference as it's loaded
		NULL,        // fasta files to sanity check reference against
		NULL,        // another way of specifying original sequences
		false,       // true -> infiles (2 args ago) contains raw seqs
		useMm,       // use memory mapping to load index?
		useShmem,    // use shared memory (not memory mapping)
		mmSweep,     // touch all the pages after memory-mapping the index
		gVerbose,    // verbose
		gVerbose);   // verbose but just for startup messages
	Timer *t = new Timer(cerr, "Time loading fw index: ", timing);
	Ebwt ebwtFw(
		ebwtBase,
		false,       // index is colorspace
		0,           // don't need entireReverse for fw index
		true,        // index is for the forward direction
		-1,          // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	t = new Timer(cerr, "Time loading bw index: ", timing);
	Ebwt ebwtBw(
		ebwtBase + ".rev",
		false,       // index is colorspace
		1,           // need entireReverse
		false,       // index is for the backward direction
		-1,          // offrate (irrelevant)
		useMm,       // whether to use memory-mapped files
		useShmem,    // whether to use shared memory
		mmSweep,     // sweep memory-mapped files
		!noRefNames, // load names?
		false,       // load SA sample?
		true,        // load ftab?
		false,       // load rstarts?
		NULL,        // reference map, or NULL if none is needed
		gVerbose,    // whether to be talkative
		gVerbose,    // talkative during initialization
		false,       // handle memory exceptions, don't pass them up
		sanity);
	delete t;
	for(int i = optind; i < argc; i++) {
	}
}
#endif
