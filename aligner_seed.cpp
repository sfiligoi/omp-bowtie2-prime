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
class SeedAligner::SeedAlignerSearchParams {
public:
	class CacheAndSeed {
	public:
		CacheAndSeed()
		: seq(NULL), n_seed_steps(0)
		, hasi0(false), fwi0(0), bwi0(0), pcache(NULL) {}

		CacheAndSeed(
			SeedSearchCache &_cache,         // local seed alignment cache
			const InstantiatedSeed& _seed,   // current instantiated seed
		        const Ebwt* ebwtFw, 	         // forward index (BWT)
			const Ebwt* ebwtBw               // backward/mirror index (BWT')

		) : seq(NULL), n_seed_steps(0)
		, hasi0(false), fwi0(0), bwi0(0), pcache(NULL) // just set a default
		{ reset(_cache,_seed,ebwtFw,ebwtBw); }

		void reset(
			SeedSearchCache &_cache,         // local seed alignment cache
			const InstantiatedSeed& _seed,   // current instantiated seed
		        const Ebwt* ebwtFw, 	         // forward index (BWT)
			const Ebwt* ebwtBw               // backward/mirror index (BWT')
		)
		{
			seq = _cache.getSeq();
			n_seed_steps = _seed.n_steps;

	                constexpr bool ltr = false; // seed_step_min > 0 i.e. n_seed_steps<0
        	        int off = abs(seed_step_min())-1;
			int ftabLen = ebwtFw->eh().ftabChars();
			hasi0 = (ftabLen > 1 && ftabLen <= maxjump());
			if(hasi0) {
				if(!ltr) {
					assert_geq(off+1, ftabLen-1);
					off = off - ftabLen + 1;
				}
				// startSearchSeedBi will need them, start prefetching now
				fwi0 = ebwtFw->ftabSeqToInt( seq, off, false);
				ebwtFw->ftabLoHiPrefetch(fwi0);
				if(ebwtBw!=NULL) {
					bwi0 = ebwtBw->ftabSeqToInt( seq, off, false);
					ebwtBw->ftabLoHiPrefetch(bwi0);
				}
			}
			pcache = &_cache; //  keep track of the big object
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
		const char *seq;                // sequence os the local seed alignment cache
		int n_seed_steps;               // steps in the current instantiated seed
		bool hasi0;
		TIndexOffU fwi0;                // Idx of fw ftab
		TIndexOffU bwi0;                // Idx of bw ftab
		SeedSearchCache *pcache; //  keep track of the big object
	};

	// create an empty bwt, tloc and bloc, with step=0
	// and constratins from seed, for initial searchSeedBi invocation
	SeedAlignerSearchParams(
		SeedSearchCache &cache,         // local seed alignment cache
		const InstantiatedSeed& seed,   // current instantiated seed
	        const Ebwt* ebwtFw, 	        // forward index (BWT)
		const Ebwt* ebwtBw)             // backward/mirror index (BWT')
	: cs(cache, seed, ebwtFw, ebwtBw)
	, step(0)
	, depth(0)
	, bwt()
	, tloc()
	, bloc()
	{}

	// create an empty bwt, tloc and bloc, with step=0
	// and constratins from seed, for initial searchSeedBi invocation
	SeedAlignerSearchParams()
	: cs()
	, step(0)
	, depth(0)
	, bwt()
	, tloc()
	, bloc()
	{}

	SeedAlignerSearchParams& operator=(const SeedAlignerSearchParams& other) = default;

	void reset(
		SeedSearchCache &cache,         // local seed alignment cache
		const InstantiatedSeed& seed,   // current instantiated seed
	        const Ebwt* ebwtFw, 	        // forward index (BWT)
		const Ebwt* ebwtBw)             // backward/mirror index (BWT')
	{
	  cs.reset(cache, seed, ebwtFw, ebwtBw);
	  step = 0;
	  depth = 0;
	  bwt.set(0,0,0,0);
	  tloc.invalidate();
	  bloc.invalidate();
	}

	void checkCV() const {
	}

	CacheAndSeed cs;      // local seed alignment cache and associated instatiated seed
	int step;             // depth into steps[] array, i.e. step_min+step
	int depth;            // recursion depth
	BwtTopBot bwt;        // The 4 BWT idxs
	SideLocus tloc;       // locus for top (perhaps unititialized)
	SideLocus bloc;       // locus for bot (perhaps unititialized)
};


//
// Some static methods for constructing some standard SeedPolicies
//

/**
 * Given a read, depth and orientation, extract a seed data structure
 * from the read and fill in the steps & zones arrays.  The Seed
 * contains the sequence and quality values.
 */

#if 0
bool
Seed::instantiate(
	const Read& read,
	const char *seq, // seed read sequence
	const Scoring& pens,
	int depth,
	int seedoffidx,
	bool fw,
	InstantiatedSeed& is) const
{
	assert(overall != NULL);
	int seedlen = len;
	if((int)read.length() < seedlen) {
		// Shrink seed length to fit read if necessary
		seedlen = (int)read.length();
	}
	assert_gt(seedlen, 0);
	is.n_steps = seedlen;
	// Fill in 'steps' and 'zones'
	//
	// The 'steps' list indicates which read character should be
	// incorporated at each step of the search process.  Often we will
	// simply proceed from one end to the other, in which case the
	// 'steps' list is ascending or descending.  In some cases (e.g.
	// the 2mm case), we might want to switch directions at least once
	// during the search, in which case 'steps' will jump in the
	// middle.  When an element of the 'steps' list is negative, this
	// indicates that the next
	//
	// The 'zones' list indicates which zone constraint is active at
	// each step.  Each element of the 'zones' list is a pair; the
	// first pair element indicates the applicable zone when
	// considering either mismatch or delete (ref gap) events, while
	// the second pair element indicates the applicable zone when
	// considering insertion (read gap) events.  When either pair
	// element is a negative number, that indicates that we are about
	// to leave the zone for good, at which point we may need to
	// evaluate whether we have reached the zone's budget.
	//
	switch(type) {
		case SEED_TYPE_EXACT: {
			is.step_min = -seedlen;
			break;
		}
		case SEED_TYPE_LEFT_TO_RIGHT: {
			is.step_min = 1;
			break;
		}
		case SEED_TYPE_RIGHT_TO_LEFT: {
			is.step_min = -seedlen;
			break;
		}
		case SEED_TYPE_INSIDE_OUT: {
			fprintf(stderr,"Unsupported SEED_TYPE_INSIDE_OUT\n");
			throw 1;
		}
		default:
			throw 1;
	}
	// Take a sweep through the seed sequence.  Consider where the Ns
	// occur and how zones are laid out.  Calculate the maximum number
	// of positions we can jump over initially (e.g. with the ftab) and
	// perhaps set this function's return value to false, indicating
	// that the arrangements of Ns prevents the seed from aligning.
	bool streak = true;
	is.maxjump = 0;
	bool ret = true;
	bool ltr = (is.step_min > 0); // true -> left-to-right
	for(size_t i = 0; i < seedlen; i++) {
		const int stepi = is.step_min+i;
		assert_neq(0, stepi);
		int off = abs(stepi)-1;
		int c = seq[off];  assert_range(0, 4, c);
		if(ltr != (stepi > 0)) // changed direction
		{
			streak = false;
		}
		if(c == 4) {
			// cons.canN always false since mms and edit both ==0
			// Seed disqualified due to arrangement of Ns
			return false;
		}
		if(streak) is.maxjump++;
	}
	is.seedoff = depth;
	is.seedoffidx = seedoffidx;
	is.fw = fw;
	is.s = *this;
	return ret;
}
#endif

bool
InstantiatedSeed::instantiateExact(
	const int seed_len,
	const Read& read,
	const char *seq) // seed read sequence
{
	int seedlen = seed_len;
	if((int)read.length() < seedlen) {
		// Shrink seed length to fit read if necessary
		seedlen = (int)read.length();
	}
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

#if 0
/**
 * Return a set of 2 seeds encapsulating a half-and-half 1mm strategy.
 */
void
Seed::oneMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0;
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
}

/**
 * Return a set of 3 seeds encapsulating search roots for:
 *
 * 1. Starting from the left-hand side and searching toward the
 *    right-hand side allowing 2 mismatches in the right half.
 * 2. Starting from the right-hand side and searching toward the
 *    left-hand side allowing 2 mismatches in the left half.
 * 3. Starting (effectively) from the center and searching out toward
 *    both the left and right-hand sides, allowing one mismatch on
 *    either side.
 *
 * This is not exhaustive.  There are 2 mismatch cases mised; if you
 * imagine the seed as divided into four successive quarters A, B, C
 * and D, the cases we miss are when mismatches occur in A and C or B
 * and D.
 */
void
Seed::twoMmSeeds(int ln, EList<Seed>& pols, Constraint& oall) {
	oall.init();
	// Seed policy 1: left-to-right search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 2: right-to-left search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(2);
	pols.back().zones[1].mmsCeil = 1; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::exact(); // not used
	pols.back().overall = &oall;
	// Seed policy 3: inside-out search
	pols.expand();
	pols.back().len = ln;
	pols.back().type = SEED_TYPE_INSIDE_OUT;
	pols.back().zones[0] = Constraint::exact();
	pols.back().zones[1] = Constraint::mmBased(1);
	pols.back().zones[1].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().zones[2] = Constraint::mmBased(1);
	pols.back().zones[2].mmsCeil = 0; // Must have used at least 1 mismatch
	pols.back().overall = &oall;
}
#endif

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
void SeedAligner::instantiateSeeds(
	const int seed_len,        // search seed length
	size_t off,                // offset into read to start extracting
	int per,                   // interval between seeds
	const Read& read,          // read to align
	bool nofw,                 // don't align forward read
	bool norc,                 // don't align revcomp read
	SeedResults& sr,           // holds all the seed hits
	int insts[3])              // counters
{
	assert_gt(read.length(), 0);
	// Check whether read has too many Ns
	const int len = seed_len;
	// Calc # seeds within read interval
	int nseeds = 1;
	if((int)read.length() - (int)off > len) {
		nseeds += ((int)read.length() - (int)off - len) / per;
	}
	insts[0] = insts[1] = insts[2] = 0;

	const int min_len = std::min<int>(len, (int)read.length());
	sr.reset(off, per, nseeds, min_len);

	// For each seed position
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			// Skip this orientation b/c user specified --nofw or --norc
			continue;
		}
		// For each seed position
		for(int i = 0; i < nseeds; i++) {
			int depth = i * per + (int)off;
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
					seed_len,
					read,
					sr.seqs(fw,i)))
				{
					// Can we fill this seed hit in from the cache?
					insts[0]++;
					insts[fw ? 1 : 2]++;
				}
			}
		}
	}
}

/**
 * We assume that all seeds are the same length.
 *
 * For each seed:
 *
 * 1. Instantiate all seeds, retracting them if necessary.
 * 2. Calculate zone boundaries for each seed
 */
void SeedAligner::searchAllSeeds(
	const Ebwt* ebwtFw,          // BWT index
	const Ebwt* ebwtBw,          // BWT' index
	const Scoring& pens,         // scoring scheme
	AlignmentCacheIface& cache,  // local cache for seed alignments
	SeedResults& sr,             // holds all the seed hits
	PerReadMetrics& prm)         // per-read metrics
{
	assert(ebwtFw != NULL);
	assert(ebwtFw->isInMemory());
	assert(sr.repOk(&cache.current()));
	ebwtFw_ = ebwtFw;
	ebwtBw_ = ebwtBw;
	sc_ = &pens;
	bwops_ = bwedits_ = 0;
	uint64_t possearches = 0, seedsearches = 0, ooms = 0;

	SeedSearchMultiCache& mcache = mcache_;
	auto& paramVec = paramVec_;

	for(int fwi = 0; fwi < 2; fwi++) {
		const bool fw = (fwi == 0);
                size_t i =0;
		// For each instantiated seed, but batched
		while (i < sr.numOffs()) {
		   const size_t ibatch_max = std::min(i+ibatch_size,sr.numOffs());
		   mcache.clear();
		   paramVec.clear();
		   // start aligning and find list of seeds to search
		   for(; i < ibatch_max; i++) {
			assert(sr.repOk(&cache.current()));
			InstantiatedSeed& is = sr.instantiatedSeed(fw, i);
			if(!is.isValid()) {
				// Cache hit in an across-read cache
				continue;
			}
			mcache.emplace_back_noresize(sr.seqs(fw,i), i, fw);
			const size_t mnr = mcache.size()-1;
			SeedSearchCache &srcache = mcache[mnr];
			{
				possearches++;
				{
					// Set seq and qual appropriately, using the seed sequences
					// and qualities already installed in SeedResults
					assert_eq(fw, is.fw);
					assert_eq(i, (int)is.seedoffidx);
					paramVec.expand_noresize();
					paramVec.back().reset(srcache, is, ebwtFw_, ebwtBw_);
					seedsearches++;
				}
			}
		   } // internal i (batch) loop

		   // do the searches
		   if (!paramVec.empty()) {
			const size_t nparams = paramVec.size();
			assert(ebwtBw_==NULL);
			searchSeedBi(ebwtFw_, sstateVec_.ptr(), bwops_, nparams, &(paramVec[0]));

			for (size_t n=0; n<nparams; n++) {
				SeedAlignerSearchState& sstate = sstateVec_[n];
				if (sstate.need_reporting) {
					SeedAlignerSearchParams& p= paramVec[n];
					assert(p.prevEdit==NULL);
					// Finished aligning seed
					const char *seq = p.cs.seq;
					auto& bwt = p.bwt;
					p.cs.pcache->addOnTheFly(seq, sr.seqs_len(), bwt.topf, bwt.botf, bwt.topb, bwt.botb);
					
				}
			}
		   }

		   // finish aligning and add to SeedResult
		   for (size_t mnr=0; mnr<mcache.size(); mnr++) {
			SeedSearchCache &srcache = mcache[mnr];
			// Tell the cache that we've started aligning, so the cache can
			// expect a series of on-the-fly updates
			int ret = srcache.beginAlign(cache, sr.seqs_len());
			if(ret == -1) {
				// Out of memory when we tried to add key to map
				ooms++;
				continue;
			}
			assert(srcache.aligning());
			if(!srcache.addAllCached()){
				// Memory exhausted during copy
				ooms++;
				continue;
			}
			srcache.finishAlign();
			assert(!srcache.aligning());
			if(srcache.qvValid()) {
				sr.add(
					srcache.getQv(),   // range of ranges in cache
					cache.current(), // cache
					mcache.getSeedOffIdx(mnr),     // seed index (from 5' end)
					mcache.getFw(mnr));   // whether seed is from forward read
			}
		   } // mnr loop
		} // external i while
	} // for fwi
}

bool SeedAligner::sanityPartial(
	const Ebwt*        ebwtFw, // BWT index
	const Ebwt*        ebwtBw, // BWT' index
	const BTDnaString& seq,
	size_t dep,
	size_t len,
	bool do1mm,
	TIndexOffU topfw,
	TIndexOffU botfw,
	TIndexOffU topbw,
	TIndexOffU botbw)
{
	tmpdnastr_.clear();
	for(size_t i = dep; i < len; i++) {
		tmpdnastr_.append(seq[i]);
	}
	TIndexOffU top_fw = 0, bot_fw = 0;
	ebwtFw->contains(tmpdnastr_, &top_fw, &bot_fw);
	assert_eq(top_fw, topfw);
	assert_eq(bot_fw, botfw);
	if(do1mm && ebwtBw != NULL) {
		tmpdnastr_.reverse();
		TIndexOffU top_bw = 0, bot_bw = 0;
		ebwtBw->contains(tmpdnastr_, &top_bw, &bot_bw);
		assert_eq(top_bw, topbw);
		assert_eq(bot_bw, botbw);
	}
	return true;
}

inline void exactSweepInit(
	const Ebwt&        ebwt,
	const BTDnaString& seq,
	const int          ftabLen,
	const size_t       len,
	size_t            &dep,
	TIndexOffU        &top, 
	TIndexOffU        &bot
	)
{
	top = bot = 0;

	const size_t left = len - dep;
	assert_gt(left, 0);
	bool doFtab = ftabLen > 1 && left >= (size_t)ftabLen;
	if(doFtab) {
		const size_t endi = len-dep-1;
		// Does N interfere with use of Ftab?
		for(size_t i = 0; i < (size_t)ftabLen; i++) {
			int c = seq[endi-i];
			if(c > 3) {
				doFtab = false;
				break;
			}
		}
	}
	if(doFtab) {
		// Use ftab
		ebwt.ftabLoHi(seq.buf(), left - ftabLen, false, top, bot);
		dep += (size_t)ftabLen;
	} else {
		// Use fchr
		int c = seq[len-dep-1];
		if(c < 4) {
			top = ebwt.fchr()[c];
			bot = ebwt.fchr()[c+1];
		}
		dep++;
	}
}

inline void exactSweepMapLF(
	const Ebwt&        ebwt,
	const BTDnaString& seq,
	const size_t       len,
	const size_t       dep,
	const SideLocus   &tloc, 
	const SideLocus  &bloc,
	TIndexOffU        &top, 
	TIndexOffU        &bot,
	uint64_t          &bwops           // Burrows-Wheeler operations
)
{
	int c = seq[len-dep-1];
	if(c > 3) {
		top = bot = 0;
	} else {
		if(bloc.valid()) {
			bwops += 2;
			top = ebwt.mapLF(tloc, c);
			bot = ebwt.mapLF(bloc, c);
		} else {
			bwops++;
			top = ebwt.mapLF1(top, tloc, c);
			if(top == OFF_MASK) {
				top = bot = 0;
			} else {
				bot = top+1;
			}
		}
	}
}


inline bool exactSweepStep(
	const Ebwt&        ebwt,    // BWT index
	const TIndexOffU   top, 
	const TIndexOffU   bot,
	const size_t       mineMax, // don't care about edit bounds > this
	SideLocus         &tloc, 
	SideLocus         &bloc,
	size_t            &mineCnt, // minimum # edits
	size_t            &nedit,
	bool              &done
	)
{
	if(bot <= top) {
		nedit++;
		if(nedit >= mineMax) {
			mineCnt = nedit;
			done = true;
		}
		return true;
	}
	INIT_LOCS(top, bot, tloc, bloc, ebwt);
	return false;
}

/**
 * Sweep right-to-left and left-to-right using exact matching.  Remember all
 * the SA ranges encountered along the way.  Report exact matches if there are
 * any.  Calculate a lower bound on the number of edits in an end-to-end
 * alignment.
 */
size_t SeedAligner::exactSweep(
	const Ebwt&        ebwt,    // BWT index
	const Read&        read,    // read to align
	const Scoring&     sc,      // scoring scheme
	bool               nofw,    // don't align forward read
	bool               norc,    // don't align revcomp read
	size_t             mineMax, // don't care about edit bounds > this
	size_t&            mineFw,  // minimum # edits for forward read
	size_t&            mineRc,  // minimum # edits for revcomp read
	bool               repex,   // report 0mm hits?
	SeedResults&       hits)    // holds all the seed hits (and exact hit)
{
	assert_gt(mineMax, 0);
	const size_t len = read.length();
	const int ftabLen = ebwt.eh().ftabChars();

	size_t nelt = 0;

	std::array<SideLocus,2> tloc;
	std::array<SideLocus,2> bloc;
	TIndexOffU top[2] = {0, 0};
	TIndexOffU bot[2] = {0, 0};

	size_t dep[2] = {0, 0};
	size_t nedit[2] = {0, 0};
	bool doInit[2] = {true, true};

	size_t prefetch_count = 0;
	bool done[2] = {nofw, norc};

	for(int fwi = 0; fwi < 2; fwi++) {
		if (!done[fwi]) {
			bool fw = (fwi == 0);
			const BTDnaString& seq = fw ? read.patFw : read.patRc;
			assert(!seq.empty());
			__builtin_prefetch(&(seq[len-1]));
			if (len>48) __builtin_prefetch(&(seq[len-49])); // HW prefetch prediction assumes forward, help it
		}
	}

	while( (dep[0] < len && !done[0]) || (dep[1] < len && !done[1]) ) {
		prefetch_count++;
		if (prefetch_count>=48) { // cache line is 64 bytes, but we may skip some deps
			for(int fwi = 0; fwi < 2; fwi++) {
				if (dep[fwi] < len && !done[fwi]) {
					bool fw = (fwi == 0);
					const BTDnaString& seq = fw ? read.patFw : read.patRc;
					const size_t left = len-dep[fwi];
					if (left>48) {
						__builtin_prefetch(&(seq[left-49])); // HW prefetch prediction assumes forward, help it
					}
				}
			}
			prefetch_count=0;
		}
		// by doing both fw in the internal loop, I give the prefetch in exactSweepStep to be effective
		for(int fwi = 0; fwi < 2; fwi++) {
			if (dep[fwi] < len && !done[fwi]) {
				bool fw = (fwi == 0);
				const BTDnaString& seq = fw ? read.patFw : read.patRc;

				if (doInit[fwi]) {
					exactSweepInit(ebwt, seq, ftabLen, len,            // in
							dep[fwi], top[fwi], bot[fwi]);          // out
					if ( exactSweepStep(ebwt, top[fwi], bot[fwi], mineMax,
							tloc[fwi], bloc[fwi],
							fw ? mineFw : mineRc,
							nedit[fwi], done[fwi]) ) {
						continue;
					}
					doInit[fwi]=false;
				}

				if (dep[fwi]< len) {
					exactSweepMapLF(ebwt, seq, len, dep[fwi], tloc[fwi], bloc[fwi],
							top[fwi], bot[fwi], bwops_);

					if ( exactSweepStep(ebwt, top[fwi], bot[fwi], mineMax,
								tloc[fwi], bloc[fwi],
								fw ? mineFw : mineRc,
								nedit[fwi], done[fwi]) ) {
						doInit[fwi]=true;
					}
					dep[fwi]++;
				}
			}
		}
	}

	for(int fwi = 0; fwi < 2; fwi++) {
		if( (!done[fwi]) && (dep[fwi] >= len) ) {
			const bool fw = (fwi == 0);

			// Set the minimum # edits
			if(fw) { mineFw = nedit[fwi]; } else { mineRc = nedit[fwi]; }
			// Done
			if(nedit[fwi] == 0 && bot[fwi] > top[fwi]) {
				if(repex) {
					// This is an exact hit
					int64_t score = len * sc.match();
					if(fw) {
						hits.addExactEeFw(top[fwi], bot[fwi], NULL, NULL, fw, score);
						assert(ebwt.contains(fw ? read.patFw : read.patRc, NULL, NULL));
					} else {
						hits.addExactEeRc(top[fwi], bot[fwi], NULL, NULL, fw, score);
						assert(ebwt.contains(fw ? read.patFw : read.patRc, NULL, NULL));
					}
				}
				nelt += (bot[fwi] - top[fwi]);
			}
		}
	}
	return nelt;
}

#if 0
/**
 * Search for end-to-end exact hit for read.  Return true iff one is found.
 */
bool SeedAligner::oneMmSearch(
	const Ebwt*        ebwtFw, // BWT index
	const Ebwt*        ebwtBw, // BWT' index
	const Read&        read,   // read to align
	const Scoring&     sc,     // scoring
	int64_t            minsc,  // minimum score
	bool               nofw,   // don't align forward read
	bool               norc,   // don't align revcomp read
	bool               repex,  // report 0mm hits?
	bool               rep1mm, // report 1mm hits?
	SeedResults&       hits)   // holds all the seed hits (and exact hit)
{
	assert(!rep1mm || ebwtBw != NULL);
	const size_t len = read.length();
	int nceil = sc.nCeil.f<int>((double)len);
	size_t ns = read.ns();
	if(ns > 1) {
		// Can't align this with <= 1 mismatches
		return false;
	} else if(ns == 1 && !rep1mm) {
		// Can't align this with 0 mismatches
		return false;
	}
	assert_geq(len, 2);
	assert(!rep1mm || ebwtBw->eh().ftabChars() == ebwtFw->eh().ftabChars());
#ifndef NDEBUG
	if(ebwtBw != NULL) {
		for(int i = 0; i < 4; i++) {
			assert_eq(ebwtBw->fchr()[i], ebwtFw->fchr()[i]);
		}
	}
#endif
	size_t halfFw = len >> 1;
	size_t halfBw = len >> 1;
	if((len & 1) != 0) {
		halfBw++;
	}
	assert_geq(halfFw, 1);
	assert_geq(halfBw, 1);
	SideLocus tloc, bloc;
	TIndexOffU t[4], b[4];   // dest BW ranges for BWT
	t[0] = t[1] = t[2] = t[3] = 0;
	b[0] = b[1] = b[2] = b[3] = 0;
	TIndexOffU tp[4], bp[4]; // dest BW ranges for BWT'
	tp[0] = tp[1] = tp[2] = tp[3] = 0;
	bp[0] = bp[1] = bp[2] = bp[3] = 0;
	TIndexOffU top = 0, bot = 0, topp = 0, botp = 0;
	// Align fw read / rc read
	bool results = false;
	for(int fwi = 0; fwi < 2; fwi++) {
		bool fw = (fwi == 0);
		if( fw && nofw) continue;
		if(!fw && norc) continue;
		// Align going right-to-left, left-to-right
		int lim = rep1mm ? 2 : 1;
		for(int ebwtfwi = 0; ebwtfwi < lim; ebwtfwi++) {
			bool ebwtfw = (ebwtfwi == 0);
			const Ebwt* ebwt  = (ebwtfw ? ebwtFw : ebwtBw);
			const Ebwt* ebwtp = (ebwtfw ? ebwtBw : ebwtFw);
			assert(rep1mm || ebwt->fw());
			const BTDnaString& seq =
				(fw ? (ebwtfw ? read.patFw : read.patFwRev) :
				      (ebwtfw ? read.patRc : read.patRcRev));
			assert(!seq.empty());
			const BTString& qual =
				(fw ? (ebwtfw ? read.qual    : read.qualRev) :
				      (ebwtfw ? read.qualRev : read.qual));
			int ftabLen = ebwt->eh().ftabChars();
			size_t nea = ebwtfw ? halfFw : halfBw;
			// Check if there's an N in the near portion
			bool skip = false;
			for(size_t dep = 0; dep < nea; dep++) {
				if(seq[len-dep-1] > 3) {
					skip = true;
					break;
				}
			}
			if(skip) {
				continue;
			}
			size_t dep = 0;
			// Align near half
			if(ftabLen > 1 && (size_t)ftabLen <= nea) {
				// Use ftab to jump partway into near half
				bool rev = !ebwtfw;
				ebwt->ftabLoHi(seq.buf(), len - ftabLen, rev, top, bot);
				if(rep1mm) {
					ebwtp->ftabLoHi(seq.buf(), len - ftabLen, rev, topp, botp);
					assert_eq(bot - top, botp - topp);
				}
				if(bot - top == 0) {
					continue;
				}
				int c = seq[len - ftabLen];
				t[c] = top; b[c] = bot;
				tp[c] = topp; bp[c] = botp;
				dep = ftabLen;
				// initialize tloc, bloc??
			} else {
				// Use fchr to jump in by 1 pos
				int c = seq[len-1];
				assert_range(0, 3, c);
				top = topp = tp[c] = ebwt->fchr()[c];
				bot = botp = bp[c] = ebwt->fchr()[c+1];
				if(bot - top == 0) {
					continue;
				}
				dep = 1;
				// initialize tloc, bloc??
			}
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
			assert(sanityPartial(ebwt, ebwtp, seq, len-dep, len, rep1mm, top, bot, topp, botp));
			bool do_continue = false;
			for(; dep < nea; dep++) {
				assert_lt(dep, len);
				int rdc = seq[len - dep - 1];
				tp[0] = tp[1] = tp[2] = tp[3] = topp;
				bp[0] = bp[1] = bp[2] = bp[3] = botp;
				if(bloc.valid()) {
					bwops_++;
					t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
					ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					top = t[rdc]; bot = b[rdc];
					if(bot <= top) {
						do_continue = true;
						break;
					}
					topp = tp[rdc]; botp = bp[rdc];
					assert(!rep1mm || bot - top == botp - topp);
				} else {
					assert_eq(bot, top+1);
					assert(!rep1mm || botp == topp+1);
					bwops_++;
					top = ebwt->mapLF1(top, tloc, rdc);
					if(top == OFF_MASK) {
						do_continue = true;
						break;
					}
					bot = top + 1;
					t[rdc] = top; b[rdc] = bot;
					tp[rdc] = topp; bp[rdc] = botp;
					assert(!rep1mm || b[rdc] - t[rdc] == bp[rdc] - tp[rdc]);
					// topp/botp stay the same
				}
				INIT_LOCS(top, bot, tloc, bloc, *ebwt);
				assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
			}
			if(do_continue) {
				continue;
			}
			// Align far half
			for(; dep < len; dep++) {
				int rdc = seq[len-dep-1];
				int quc = qual[len-dep-1];
				if(rdc > 3 && nceil == 0) {
					break;
				}
				tp[0] = tp[1] = tp[2] = tp[3] = topp;
				bp[0] = bp[1] = bp[2] = bp[3] = botp;
				int clo = 0, chi = 3;
				bool match = true;
				if(bloc.valid()) {
					bwops_++;
					t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
					ebwt->mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					match = rdc < 4;
					top = t[rdc]; bot = b[rdc];
					topp = tp[rdc]; botp = bp[rdc];
				} else {
					assert_eq(bot, top+1);
					assert(!rep1mm || botp == topp+1);
					bwops_++;
					clo = ebwt->mapLF1(top, tloc);
					match = (clo == rdc);
					assert_range(-1, 3, clo);
					if(clo < 0) {
						break; // Hit the $
					} else {
						t[clo] = top;
						b[clo] = bot = top + 1;
					}
					bp[clo] = botp;
					tp[clo] = topp;
					assert(!rep1mm || bot - top == botp - topp);
					assert(!rep1mm || b[clo] - t[clo] == bp[clo] - tp[clo]);
					chi = clo;
				}
				//assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
				if(rep1mm && (ns == 0 || rdc > 3)) {
					for(int j = clo; j <= chi; j++) {
						if(j == rdc || b[j] == t[j]) {
							// Either matches read or isn't a possibility
							continue;
						}
						// Potential mismatch - next, try
						size_t depm = dep + 1;
						TIndexOffU topm = t[j], botm = b[j];
						TIndexOffU topmp = tp[j], botmp = bp[j];
						assert_eq(botm - topm, botmp - topmp);
						TIndexOffU tm[4], bm[4];   // dest BW ranges for BWT
						tm[0] = t[0]; tm[1] = t[1];
						tm[2] = t[2]; tm[3] = t[3];
						bm[0] = b[0]; bm[1] = t[1];
						bm[2] = b[2]; bm[3] = t[3];
						TIndexOffU tmp[4], bmp[4]; // dest BW ranges for BWT'
						tmp[0] = tp[0]; tmp[1] = tp[1];
						tmp[2] = tp[2]; tmp[3] = tp[3];
						bmp[0] = bp[0]; bmp[1] = tp[1];
						bmp[2] = bp[2]; bmp[3] = tp[3];
						SideLocus tlocm, blocm;
						INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
						for(; depm < len; depm++) {
							int rdcm = seq[len - depm - 1];
							tmp[0] = tmp[1] = tmp[2] = tmp[3] = topmp;
							bmp[0] = bmp[1] = bmp[2] = bmp[3] = botmp;
							if(blocm.valid()) {
								bwops_++;
								tm[0] = tm[1] = tm[2] = tm[3] =
								bm[0] = bm[1] = bm[2] = bm[3] = 0;
								ebwt->mapBiLFEx(tlocm, blocm, tm, bm, tmp, bmp);
								SANITY_CHECK_4TUP(tm, bm, tmp, bmp);
								topm = tm[rdcm]; botm = bm[rdcm];
								topmp = tmp[rdcm]; botmp = bmp[rdcm];
								if(botm <= topm) {
									break;
								}
							} else {
								assert_eq(botm, topm+1);
								assert_eq(botmp, topmp+1);
								bwops_++;
								topm = ebwt->mapLF1(topm, tlocm, rdcm);
								if(topm == OFF_MASK) {
									break;
								}
								botm = topm + 1;
								// topp/botp stay the same
							}
							INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
						}
						if(depm == len) {
							// Success; this is a 1MM hit
							size_t off5p = dep;  // offset from 5' end of read
							size_t offstr = dep; // offset into patFw/patRc
							if(fw == ebwtfw) {
								off5p = len - off5p - 1;
							}
							if(!ebwtfw) {
								offstr = len - offstr - 1;
							}
							Edit e((uint32_t)off5p, j, rdc, EDIT_TYPE_MM, false);
							results = true;
							int64_t score = (len - 1) * sc.match();
							int pen = sc.score(rdc, (int)(1 << j), quc - 33);
							score += pen;
							bool valid = score >= minsc;
							if(valid) { // Note: Should be always valid, but needs to prove it
								TIndexOffU toprep = ebwtfw ? topm : topmp;
								TIndexOffU botrep = ebwtfw ? botm : botmp;
								assert_eq(toprep, toptmp);
								assert_eq(botrep, bottmp);
								hits.add1mmEe(toprep, botrep, &e, NULL, fw, score);
							}
						}
					}
				}
				if(bot > top && match) {
					assert_lt(rdc, 4);
					if(dep == len-1) {
						// Success; this is an exact hit
						if(ebwtfw && repex) {
							if(fw) {
								results = true;
								int64_t score = len * sc.match();
								hits.addExactEeFw(
									ebwtfw ? top : topp,
									ebwtfw ? bot : botp,
									NULL, NULL, fw, score);
								assert(ebwtFw->contains(seq, NULL, NULL));
							} else {
								results = true;
								int64_t score = len * sc.match();
								hits.addExactEeRc(
									ebwtfw ? top : topp,
									ebwtfw ? bot : botp,
									NULL, NULL, fw, score);
								assert(ebwtFw->contains(seq, NULL, NULL));
							}
						}
						break; // End of far loop
					} else {
						INIT_LOCS(top, bot, tloc, bloc, *ebwt);
						assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
					}
				} else {
					break; // End of far loop
				}
			} // for(; dep < len; dep++)
		} // for(int ebwtfw = 0; ebwtfw < 2; ebwtfw++)
	} // for(int fw = 0; fw < 2; fw++)
	return results;
}
#endif

/**
 * Get tloc, bloc ready for the next step.  If the new range is under
 * the ceiling.
 */
inline void
SeedAligner::nextLocsBi(
        const Ebwt* ebwt,             // forward index (BWT)
	const int seed_step,          // current instantiated seed step
	SideLocus& tloc,              // top locus
	SideLocus& bloc,              // bot locus
	TIndexOffU topf,              // top in BWT
	TIndexOffU botf,              // bot in BWT
	TIndexOffU topb,              // top in BWT'
	TIndexOffU botb               // bot in BWT'
	)
{
	assert_gt(botf, 0);
	assert(botb > 0);
	assert(botf-topf == botb-topb);
	{
		if(botf - topf == 1) {
			// Already down to 1 row; just init top locus
			tloc.initFromRow(topf, ebwt->eh(), ebwt->ebwt());
			bloc.invalidate();
		} else {
			SideLocus::initFromTopBot(
				topf, botf, ebwt->eh(), ebwt->ebwt(), tloc, bloc);
			assert(bloc.valid());
		}
	}
	assert(botf - topf == 1 ||  bloc.valid());
	assert(botf - topf > 1  || !bloc.valid());
}

// return true, if we are already done
bool
SeedAligner::startSearchSeedBi(
                        const Ebwt* ebwt,       // forward index (BWT)
			SeedAligner::SeedAlignerSearchParams &p)
{
	const char *seq = p.cs.seq;
	const int n_seed_steps = p.cs.n_seed_steps;

	assert_gt(n_seed_steps, 0);
	if(p.step == n_seed_steps) {
		return true;
	}
#ifndef NDEBUG
	if(p.depth > 0) {
		assert(p.bwt.botf - p.bwt.topf == 1 ||  p.bloc.valid());
		assert(p.bwt.botf - p.bwt.topf > 1  || !p.bloc.valid());
	}
#endif
	if(p.step == 0) {
		// Just starting
		assert(p.prevEdit == NULL);
		assert(!p.tloc.valid());
		assert(!p.bloc.valid());
		const int seed_step_min = p.cs.seed_step_min();
		int off = abs(seed_step_min)-1;
		// Check whether/how far we can jump using ftab or fchr
		int ftabLen = ebwt->eh().ftabChars();
		if(p.cs.hasi0) { //if(ftabLen > 1 && ftabLen <= p.cs.maxjump)
			ebwt->ftabLoHi(p.cs.fwi0, p.bwt.topf, p.bwt.botf);
			if(p.bwt.botf - p.bwt.topf == 0) return true;
			p.step += ftabLen;
		} else if(p.cs.maxjump() > 0) {
			// Use fchr
			int c = seq[off];
			assert_range(0, 3, c);
			p.bwt.topf = p.bwt.topb = ebwt->fchr()[c];
			p.bwt.botf = p.bwt.botb = ebwt->fchr()[c+1];
			if(p.bwt.botf - p.bwt.topf == 0) return true;
			p.step++;
		} else {
			assert_eq(0, p.cs.maxjump());
			p.bwt.topf = p.bwt.topb = 0;
			p.bwt.botf = p.bwt.botb = ebwt->fchr()[4];
		}
		if(p.step == n_seed_steps) {
			return true;
		}
		nextLocsBi(ebwt, seed_step_min+p.step, p.tloc, p.bloc, p.bwt);
		assert(p.tloc.valid());
	} else assert(p.prevEdit != NULL);
	assert(p.tloc.valid());
	assert(p.bwt.botf - p.bwt.topf == 1 ||  p.bloc.valid());
	assert(p.bwt.botf - p.bwt.topf > 1  || !p.bloc.valid());
	assert_geq(p.step, 0);

	return false;
}

class SeedAlignerSearchSave {
public:
	SeedAlignerSearchSave(
		Constraint &cons, Constraint &ovCons,
		SideLocus &tloc, SideLocus &bloc)
	: orgCons(cons), orgOvCons(ovCons)
	, orgTloc(tloc), orgBloc(bloc)
	, oldCons(cons), oldOvCons(ovCons)
	, oldTloc(tloc), oldBloc(bloc)
	{}

	~SeedAlignerSearchSave() {
		orgCons = oldCons; orgOvCons = oldOvCons;
		orgTloc = oldTloc; orgBloc = oldBloc;
	}

private:
	// reference to the original variables
	Constraint &orgCons;
	Constraint &orgOvCons;
	SideLocus &orgTloc;
	SideLocus &orgBloc;

	// copy of the originals
	const Constraint oldCons;
	const Constraint oldOvCons;
	const SideLocus oldTloc;
	const SideLocus oldBloc;
};

#if 0
// State used for recursion
class SeedAlignerSearchRecState {
public:
	SeedAlignerSearchRecState(
			int j,
			int c,
			const SeedAlignerSearchState &sstate,
			DoublyLinkedList<Edit> *prevEdit)
	: bwt(sstate.tf[j], sstate.bf[j], sstate.tb[j], sstate.bb[j])
	, edit(sstate.off, j, c, EDIT_TYPE_MM, false)
	, editl()
	, _prevEdit(prevEdit)
	{
		assert(_prevEdit == NULL || _prevEdit->next == NULL);
		editl.payload = edit;
		if(_prevEdit != NULL) {
			_prevEdit->next = &editl;
			editl.prev = _prevEdit;
		}
		assert(editl.next == NULL);
	}

	~SeedAlignerSearchRecState() {
		if(_prevEdit != NULL) _prevEdit->next = NULL;
	}

	BwtTopBot bwt;
	Edit edit;
	DoublyLinkedList<Edit> editl;
private:
	DoublyLinkedList<Edit> *_prevEdit;
};
#endif

/**
 * Given a seed, search.  Assumes zone 0 = no backtracking.
 *
 * Return a list of Seed hits.
 * 1. Edits
 * 2. Bidirectional BWT range(s) on either end
 */
void
SeedAligner::searchSeedBi(
                        const Ebwt* ebwt,       // forward index (BWT)
                        SeedAlignerSearchState* sstateVec,
                        uint64_t& bwops_,         // Burrows-Wheeler operations
                        const uint32_t nparams, SeedAlignerSearchParams paramVec[])
{
	uint32_t nleft = nparams; // will keep track of how many are not done yet

	{
           uint32_t ncompleted = 0;
	   for (uint32_t n=0; n<nparams; n++) {
		SeedAlignerSearchParams& p= paramVec[n];
		SeedAlignerSearchState& sstate = sstateVec[n];
		const bool done = startSearchSeedBi(ebwt, p);
		sstate.done = done;
		sstate.need_reporting = false;
		if(done) {
		        if(p.step == (int)p.cs.n_seed_steps) {
                		// Finished aligning seed
				p.checkCV();
				sstate.need_reporting = true;
			}
			ncompleted++;
		}
	    }
	    nleft -= ncompleted;
	}


	while (nleft>0) {
	   // Note: We can do the params in any order we want
	   // but we must do the steps inside the same param in order
	   // We still want to do them sequentially, not in parallel, 
	   // to give time for the prefetch to do its job.
	   // Will loop over all of them, and just check which ones are invalid
	   uint32_t bwops = 0;
           uint32_t ncompleted = 0;
           for (uint32_t n=0; n<nparams; n++) {
                SeedAlignerSearchState& sstate = sstateVec[n];
		if (sstate.done) continue;

		SeedAlignerSearchParams& p= paramVec[n];
		const int n_seed_steps = p.cs.n_seed_steps;
		if (p.step >= (int) n_seed_steps) {
			sstate.done = true;
			ncompleted++;
			continue;
		}
		const int seed_step_min = p.cs.seed_step_min();
		size_t i = p.step; // call the stepIdx i for historical reasons
		p.step++; // get ready for the next iteration

		const char *seq = p.cs.seq;

		assert_gt(p.bwt.botf, p.bwt.topf);
		assert(p.bwt.botf - p.bwt.topf == 1 ||  p.bloc.valid());
		assert(p.bwt.botf - p.bwt.topf > 1  || !p.bloc.valid());
		assert(p.bwt.botf-p.bwt.topf == p.bwt.botb-p.bwt.topb);
		assert(p.tloc.valid());
		SeedAlignerSearchWorkState wstate(seed_step_min+i, p.bwt);
		__builtin_prefetch(&(seq[wstate.off]));
		if(p.bloc.valid()) {
			// Range delimited by tloc/bloc has size >1.  If size == 1,
			// we use a simpler query (see if(!bloc.valid()) blocks below)
			bwops++;
			ebwt->mapBiLFEx(p.tloc, p.bloc, wstate.t, wstate.b, wstate.tp, wstate.bp);
		}
		int c = seq[wstate.off];  assert_range(0, 4, c);
		//
		if(c == 4) { // couldn't handle the N
			sstate.done = true;
			ncompleted++;
			continue;
		}
		if(!p.bloc.valid()) {
			assert(wstate.bp[c] == wstate.tp[c]+1);
			// Range delimited by tloc/bloc has size 1
			bwops++;
			wstate.t[c] = ebwt->mapLF1(wstate.ntop, p.tloc, c);
			if(wstate.t[c] == OFF_MASK) {
				sstate.done = true;
				ncompleted++;
				continue;
			}
			assert_geq(wstate.t[c], ebwt->fchr()[c]);
			assert_lt(wstate.t[c],  ebwt->fchr()[c+1]);
			wstate.b[c] = wstate.t[c]+1;
			assert_gt(wstate.b[c], 0);
		}
		assert(wstate.bf[c]-wstate.tf[c] == wstate.bb[c]-wstate.tb[c]);
		if(wstate.b[c] == wstate.t[c]) {
			sstate.done = true;
			ncompleted++;
			continue;
		}
		p.bwt.set(wstate.tf[c], wstate.bf[c], wstate.tb[c], wstate.bb[c]);
		if(i+1 == n_seed_steps) {
			p.checkCV();
			sstate.need_reporting = true;
			sstate.done = true;
			ncompleted++;
			continue;
		}
		nextLocsBi(ebwt, seed_step_min+i+1, p.tloc, p.bloc, p.bwt);
	   } // for n
	   nleft -= ncompleted;
	   bwops_ += bwops;
	} // while

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
