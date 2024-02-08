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

#include <cassert>
#include <vector>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <thread>
#include <time.h>
#include <utility>
#include <execution>

#ifndef _WIN32
#include <signal.h>
#endif

#include "alphabet.h"
#include "assert_helpers.h"
#include "endian_swap.h"
#include "bt2_idx.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "aln_sink.h"
#include "pat.h"
#include "threading.h"
#include "ds.h"
#include "aligner_metrics.h"
#include "sam.h"
#include "aligner_seed.h"
#include "aligner_seed_policy.h"
#include "aligner_driver.h"
#include "aligner_sw.h"
#include "aligner_sw_driver.h"
#include "aligner_cache.h"
#include "util.h"
#include "pe.h"
#include "simple_func.h"
#include "presets.h"
#include "opts.h"
#include "outq.h"
#include "aligner_seed2.h"
#include "bt2_search.h"

using namespace std;

static int FNAME_SIZE;
static std::atomic<int> thread_counter;
static EList<string> mates1;  // mated reads (first mate)
static EList<string> mates2;  // mated reads (second mate)
static EList<string> mates12; // mated reads (1st/2nd interleaved in 1 file)
static string adjIdxBase;
int gVerbose;             // be talkative
static bool startVerbose; // be talkative at startup
int gQuiet;               // print nothing but the alignments
static int sanityCheck;   // enable expensive sanity checks
static int format;        // default read format is FASTQ
static bool interleaved;  // reads are interleaved
static string origString; // reference text, or filename(s)
static int seed;          // srandom() seed
static int timing;        // whether to report basic timing data
static bool allHits;      // for multihits, report just one
static bool showVersion;  // just print version and quit?
static int ipause;        // pause before maching?
static int gTrim5;        // amount to trim from 5' end
static int gTrim3;        // amount to trim from 3' end
static pair<short, size_t> trimTo; // trim reads exceeding given length from either 3' or 5'-end
static int offRate;       // keep default offRate
static bool solexaQuals;  // quality strings are solexa quals, not phred, and subtract 64 (not 33)
static bool phred64Quals; // quality chars are phred, but must subtract 64 (not 33)
static bool integerQuals; // quality strings are space-separated strings of integers, not ASCII
static int nthreads;      // number of pthreads operating concurrently
static int thread_ceiling;// maximum number of threads user wants bowtie to use
static int outType;       // style of output
static bool noRefNames;   // true -> print reference indexes; not names
static uint32_t khits;    // number of hits per read; >1 is much slower
static uint32_t mhits;    // don't report any hits if there are > mhits
static int partitionSz;   // output a partitioning key in first field
static int readsPerBatch; // # reads to read from input file at once
static bool fileParallel; // separate threads read separate input files in parallel
static bool useShmem;     // use shared memory to hold the index
static bool useMm;        // use memory-mapped files to hold the index
static bool mmSweep;      // sweep through memory-mapped files immediately after mapping
int gMinInsert;           // minimum insert size
int gMaxInsert;           // maximum insert size
bool gMate1fw;            // -1 mate aligns in fw orientation on fw strand
bool gMate2fw;            // -2 mate aligns in rc orientation on fw strand
bool gFlippedMatesOK;     // allow mates to be in wrong order
bool gDovetailMatesOK;    // allow one mate to extend off the end of the other
bool gContainMatesOK;     // allow one mate to contain the other in PE alignment
bool gOlapMatesOK;        // allow mates to overlap in PE alignment
bool gExpandToFrag;       // incr max frag length to =larger mate len if necessary
bool gReportDiscordant;   // find and report discordant paired-end alignments
bool gReportMixed;        // find and report unpaired alignments for paired reads
static uint32_t cacheLimit;      // ranges w/ size > limit will be cached
static uint32_t cacheSize;       // # words per range cache
bool gNofw; // don't align fw orientation of read
bool gNorc; // don't align rc orientation of read
static uint32_t fastaContLen;
static uint32_t fastaContFreq;
static bool hadoopOut; // print Hadoop status and summary messages
static bool fullRef;
static bool samTruncQname; // whether to truncate QNAME to 255 chars
static bool samAppendComment; // append FASTA/FASTQ comment to SAM record
static bool samOmitSecSeqQual; // omit SEQ/QUAL for 2ndary alignments?
static bool samNoUnal; // don't print records for unaligned reads
static bool samNoHead; // don't print any header lines in SAM output
static bool samNoSQ;   // don't print @SQ header lines
static bool sam_print_as;
static bool sam_print_xs;  // XS:i
static bool sam_print_xss; // Xs:i and Ys:i
static bool sam_print_yn;  // YN:i and Yn:i
static bool sam_print_xn;
static bool sam_print_x0;
static bool sam_print_x1;
static bool sam_print_xm;
static bool sam_print_xo;
static bool sam_print_xg;
static bool sam_print_nm;
static bool sam_print_md;
static bool sam_print_yf;
static bool sam_print_yi;
static bool sam_print_ym;
static bool sam_print_yp;
static bool sam_print_yt;
static bool sam_print_ys;
static bool sam_print_zs;
static bool sam_print_xr;
static bool sam_print_xt;
static bool sam_print_xd;
static bool sam_print_xu;
static bool sam_print_yl;
static bool sam_print_ye;
static bool sam_print_yu;
static bool sam_print_xp;
static bool sam_print_yr;
static bool sam_print_zb;
static bool sam_print_zr;
static bool sam_print_zf;
static bool sam_print_zm;
static bool sam_print_zi;
static bool sam_print_zp;
static bool sam_print_zu;
static bool sam_print_zt;
static bool preserve_tags;     // Only applies when aligning BAM files
static bool align_paired_reads; // Process only the paired reads in BAM file
static bool gSeedLenIsSet;
static bool qcFilter;
bool gReportOverhangs;        // false -> filter out alignments that fall off the end of a reference sequence
static string rgid;           // ID: setting for @RG header line
static string rgs;            // SAM outputs for @RG header line
static string rgs_optflag;    // SAM optional flag to add corresponding to @RG ID
static bool msample;          // whether to report a random alignment when maxed-out via -m/-M
int      gGapBarrier;         // # diags on top/bot only to be entered diagonally
int gDefaultSeedLen;
static EList<string> qualities;
static EList<string> qualities1;
static EList<string> qualities2;
static string polstr;         // temporary holder for policy string
static bool  msNoCache;       // true -> disable local cache
static int   bonusMatchType;  // how to reward matches
static int   bonusMatch;      // constant reward if bonusMatchType=constant
static int   penMmcType;      // how to penalize mismatches
static int   penMmcMax;       // max mm penalty
static int   penMmcMin;       // min mm penalty
static int   penNType;        // how to penalize Ns in the read
static int   penN;            // constant if N pelanty is a constant
static bool  penNCatPair;     // concatenate mates before N filtering?
static bool  localAlign;      // do local alignment in DP steps
static bool  noisyHpolymer;   // set to true if gap penalties should be reduced to be consistent with a sequencer that under- and overcalls homopolymers
static int   penRdGapConst;   // constant cost of extending a gap in the read
static int   penRfGapConst;   // constant cost of extending a gap in the reference
static int   penRdGapLinear;  // coeff of linear term for cost of gap extension in read
static int   penRfGapLinear;  // coeff of linear term for cost of gap extension in ref
static SimpleFunc scoreMin;   // minimum valid score as function of read len
static SimpleFunc nCeil;      // max # Ns allowed as function of read len
static SimpleFunc msIval;     // interval between seeds as function of read len
static double descConsExp;    // how to adjust score minimum as we descent further into index-assisted alignment
static bool descPrioritizeRoots; // whether to prioritize search roots with scores
static size_t descLanding;    // don't place a search root if it's within this many positions of end
static SimpleFunc descentTotSz;    // maximum space a DescentDriver can use in bytes
static SimpleFunc descentTotFmops; // maximum # FM ops a DescentDriver can perform
static int    multiseedMms;   // mismatches permitted in a multiseed seed
static int    multiseedLen;   // length of multiseed seeds
static size_t multiseedOff;   // offset to begin extracting seeds
static uint32_t seedCacheLocalMB;   // # MB to use for non-shared seed alignment cacheing
static uint32_t seedCacheCurrentMB; // # MB to use for current-read seed hit cacheing
static uint32_t exactCacheCurrentMB; // # MB to use for current-read seed hit cacheing
static size_t maxhalf;        // max width on one side of DP table
static bool scUnMapped;       // consider soft-clipped bases unmapped when calculating TLEN
static bool xeq;              // use X/= instead of M in CIGAR string
static size_t maxIters;       // stop after this many extend loop iterations
static size_t maxUg;          // stop after this many ungap extends
static size_t maxDp;          // stop after this many DPs
static size_t maxItersIncr;   // amt to add to maxIters for each -k > 1
static size_t maxEeStreak;    // stop after this many end-to-end fails in a row
static size_t maxUgStreak;    // stop after this many ungap fails in a row
static size_t maxDpStreak;    // stop after this many dp fails in a row
static size_t maxStreakIncr;  // amt to add to streak for each -k > 1
static size_t maxMateStreak;  // stop seed range after this many mate-find fails
static bool doExtend;         // extend seed hits
static bool enable8;          // use 8-bit SSE where possible?
static size_t cminlen;        // longer reads use checkpointing
static size_t cpow2;          // checkpoint interval log2
static bool doTri;            // do triangular mini-fills?
static string defaultPreset;  // default preset; applied immediately
static bool ignoreQuals;      // all mms incur same penalty, regardless of qual
static string wrapper;        // type of wrapper script, so we can print correct usage
static EList<string> queries; // list of query files
static string outfile;        // write SAM output to this file
static int mapqv;             // MAPQ calculation version
static int tighten;           // -M tighten mode (0=none, 1=best, 2=secbest+1)
static size_t do1mmMinLen;    // length below which we disable 1mm e2e search
static size_t seedBoostThresh;   // if average non-zero position has more than this many elements
static size_t nSeedRounds;    // # seed rounds
static bool reorder;          // true -> reorder SAM recs in -p mode
static float sampleFrac;      // only align random fraction of input reads
static bool arbitraryRandom;  // pseudo-randoms no longer a function of read properties
static bool bowtie2p5;
static string logDps;         // log seed-extend dynamic programming problems
static string logDpsOpp;      // log mate-search dynamic programming problems

static string bt2index;      // read Bowtie 2 index from files with this prefix
static EList<pair<int, string> > extra_opts;
static size_t extra_opts_cur;

#ifdef USE_SRA
static EList<string> sra_accs;
#endif

#define DMAX std::numeric_limits<double>::max()

static void set_format(int &current_format, file_format format) {
	if (current_format == UNKNOWN)
		current_format = format;
	else {
		std::cerr << file_format_names[current_format] << " and "
			  << file_format_names[format] << " formats are "
			  << "mutually exclusive." << std::endl;
		exit(1);
	}
}

static void resetOptions() {
	mates1.clear();
	mates2.clear();
	mates12.clear();
	adjIdxBase	    = "";
	gVerbose            = 0;
	startVerbose	    = 0;
	gQuiet		    = false;
	sanityCheck	    = 0;	// enable expensive sanity checks
	format		    = UNKNOWN;	// default read format is FASTQ
	interleaved	    = false;	// reads are not interleaved by default
	origString	    = "";	// reference text, or filename(s)
	seed		    = 0;	// srandom() seed
	timing		    = 0;	// whether to report basic timing data
	allHits		    = false;	// for multihits, report just one
	showVersion	    = false;	// just print version and quit?
	ipause		    = 0;	// pause before maching?
	gTrim5		    = 0;	// amount to trim from 5' end
	gTrim3		    = 0;	// amount to trim from 3' end
	trimTo		    = pair<short, size_t>(5, 0);	// default: don't do any trimming
	offRate		    = -1;	// keep default offRate
	solexaQuals	    = false;	// quality strings are solexa quals, not phred, and subtract 64 (not 33)
	phred64Quals	    = false;	// quality chars are phred, but must subtract 64 (not 33)
	integerQuals	    = false;	// quality strings are space-separated strings of integers, not ASCII
	nthreads	    = 1;	// number of pthreads operating concurrently
	thread_ceiling	    = 0;	// max # threads user asked for
	FNAME_SIZE	    = 4096;
	outType		    = OUTPUT_SAM;	// style of output
	noRefNames	    = false;	// true -> print reference indexes; not names
	khits		    = 1;	// number of hits per read; >1 is much slower
	mhits		    = 50;	// stop after finding this many alignments+1
	partitionSz	    = 0;	// output a partitioning key in first field
	readsPerBatch	    = 16;	// # reads to read from input file at once
	fileParallel	    = false;	// separate threads read separate input files in parallel
	useShmem	    = false;	// use shared memory to hold the index
	useMm		    = false;	// use memory-mapped files to hold the index
	mmSweep		    = false;	// sweep through memory-mapped files immediately after mapping
	gMinInsert	    = 0;	// minimum insert size
	gMaxInsert	    = 500;	// maximum insert size
	gMate1fw	    = true;	// -1 mate aligns in fw orientation on fw strand
	gMate2fw	    = false;	// -2 mate aligns in rc orientation on fw strand
	gFlippedMatesOK     = false;	// allow mates to be in wrong order
	gDovetailMatesOK    = false;	// allow one mate to extend off the end of the other
	gContainMatesOK     = true;	// allow one mate to contain the other in PE alignment
	gOlapMatesOK        = true;	// allow mates to overlap in PE alignment
	gExpandToFrag       = true;	// incr max frag length to =larger mate len if necessary
	gReportDiscordant   = true;	// find and report discordant paired-end alignments
	gReportMixed        = true;	// find and report unpaired alignments for paired reads

	cacheLimit	    = 5;	// ranges w/ size > limit will be cached
	cacheSize	    = 0;	// # words per range cache
	gNofw		    = false;	// don't align fw orientation of read
	gNorc		    = false;	// don't align rc orientation of read
	fastaContLen	    = 0;
	fastaContFreq	    = 0;
	hadoopOut	    = false;	// print Hadoop status and summary messages
	fullRef		    = false;	// print entire reference name instead of just up to 1st space
	samTruncQname       = true;	// whether to truncate QNAME to 255 chars
	samAppendComment    = false;	// append FASTA/Q comment to SAM record
	samOmitSecSeqQual   = false;	// omit SEQ/QUAL for 2ndary alignments?
	samNoUnal           = false;	// omit SAM records for unaligned reads
	samNoHead	    = false;	// don't print any header lines in SAM output
	samNoSQ		    = false;	// don't print @SQ header lines
	sam_print_as        = true;
	sam_print_xs        = true;
	sam_print_xss       = false;	// Xs:i and Ys:i
	sam_print_yn        = false;	// YN:i and Yn:i
	sam_print_xn        = true;
	sam_print_x0        = true;
	sam_print_x1        = true;
	sam_print_xm        = true;
	sam_print_xo        = true;
	sam_print_xg        = true;
	sam_print_nm        = true;
	sam_print_md        = true;
	sam_print_yf        = true;
	sam_print_yi        = false;
	sam_print_ym        = false;
	sam_print_yp        = false;
	sam_print_yt        = true;
	sam_print_ys        = true;
	sam_print_zs        = false;
	sam_print_xr        = false;
	sam_print_xt        = false;
	sam_print_xd        = false;
	sam_print_xu        = false;
	sam_print_yl        = false;
	sam_print_ye        = false;
	sam_print_yu        = false;
	sam_print_xp        = false;
	sam_print_yr        = false;
	sam_print_zb        = false;
	sam_print_zr        = false;
	sam_print_zf        = false;
	sam_print_zm        = false;
	sam_print_zi        = false;
	sam_print_zp        = false;
	sam_print_zu        = false;
	sam_print_zt        = false;
	preserve_tags       = false;
	align_paired_reads  = false;
	gSeedLenIsSet	    = false;
	gDefaultSeedLen	    = DEFAULT_SEEDLEN;
	qcFilter            = false;	// don't believe upstream qc by default
	rgid		    = "";	// SAM outputs for @RG header line
	rgs		    = "";	// SAM outputs for @RG header line
	rgs_optflag	    = "";	// SAM optional flag to add corresponding to @RG ID
	msample		    = true;
	gGapBarrier	    = 4;	// disallow gaps within this many chars of either end of alignment
	qualities.clear();
	qualities1.clear();
	qualities2.clear();
	polstr.clear();
	msNoCache	    = true;	// true -> disable local cache
	bonusMatchType	    = DEFAULT_MATCH_BONUS_TYPE;
	bonusMatch	    = DEFAULT_MATCH_BONUS;
	penMmcType	    = DEFAULT_MM_PENALTY_TYPE;
	penMmcMax	    = DEFAULT_MM_PENALTY_MAX;
	penMmcMin	    = DEFAULT_MM_PENALTY_MIN;
	penNType	    = DEFAULT_N_PENALTY_TYPE;
	penN		    = DEFAULT_N_PENALTY;
	penNCatPair	    = DEFAULT_N_CAT_PAIR;	// concatenate mates before N filtering?
	localAlign	    = false;	// do local alignment in DP steps
	noisyHpolymer	    = false;
	penRdGapConst	    = DEFAULT_READ_GAP_CONST;
	penRfGapConst	    = DEFAULT_REF_GAP_CONST;
	penRdGapLinear	    = DEFAULT_READ_GAP_LINEAR;
	penRfGapLinear	    = DEFAULT_REF_GAP_LINEAR;
	scoreMin.init  (SIMPLE_FUNC_LINEAR, DEFAULT_MIN_CONST,   DEFAULT_MIN_LINEAR);
	nCeil.init     (SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
	msIval.init    (SIMPLE_FUNC_LINEAR, 1.0f, DMAX, DEFAULT_IVAL_B, DEFAULT_IVAL_A);
	descConsExp	    = 2.0;
	descPrioritizeRoots = false;
	descLanding	    = 20;
	descentTotSz.init(SIMPLE_FUNC_LINEAR, 1024.0, DMAX, 0.0, 1024.0);
	descentTotFmops.init(SIMPLE_FUNC_LINEAR, 100.0, DMAX, 0.0, 10.0);
	multiseedMms	    = DEFAULT_SEEDMMS;
	multiseedLen	    = gDefaultSeedLen;
	multiseedOff	    = 0;
	seedCacheLocalMB    = 32;	// # MB to use for non-shared seed alignment cacheing
	seedCacheCurrentMB  = 20;	// # MB to use for current-read seed hit cacheing
	exactCacheCurrentMB = 20;	// # MB to use for current-read seed hit cacheing
	maxhalf		    = 15;	// max width on one side of DP table
	scUnMapped	    = false;	// consider soft clipped bases unmapped when calculating TLEN
	xeq		    = false;	// use =/X instead of M in CIGAR string
	maxIters	    = 400;	// max iterations of extend loop
	maxUg		    = 300;	// stop after this many ungap extends
	maxDp		    = 300;	// stop after this many dp extends
	maxItersIncr	    = 20;	// amt to add to maxIters for each -k > 1
	maxEeStreak	    = 15;	// stop after this many end-to-end fails in a row
	maxUgStreak	    = 15;	// stop after this many ungap fails in a row
	maxDpStreak	    = 15;	// stop after this many dp fails in a row
	maxStreakIncr	    = 10;	// amt to add to streak for each -k > 1
	maxMateStreak	    = 10;	// in PE: abort seed range after N mate-find fails
	doExtend	    = true;	// do seed extensions
	enable8		    = true;	// use 8-bit SSE where possible?
	cminlen		    = 2000;	// longer reads use checkpointing
	cpow2		    = 4;	// checkpoint interval log2
	doTri		    = false;	// do triangular mini-fills?
	defaultPreset	    = "sensitive%LOCAL%";	// default preset; applied immediately
	extra_opts.clear();
	extra_opts_cur	    = 0;
	bt2index.clear();       // read Bowtie 2 index from files with this prefix
	ignoreQuals	    = false;	// all mms incur same penalty, regardless of qual
	wrapper.clear();        // type of wrapper script, so we can print correct usage
	queries.clear();        // list of query files
	outfile.clear();        // write SAM output to this file
	mapqv		    = 2;	// MAPQ calculation version
	tighten		    = 3;	// -M tightening mode
	seedBoostThresh	    = 300;	// if average non-zero position has more than this many elements
	nSeedRounds	    = 2;	// # rounds of seed searches to do for repetitive reads
	do1mmMinLen	    = 60;	// length below which we disable 1mm search
	reorder		    = false;	// reorder SAM records with -p > 1
	sampleFrac	    = 1.1f;	// align all reads
	arbitraryRandom	    = false;	// let pseudo-random seeds be a function of read properties
	bowtie2p5	    = false;
	logDps.clear();         // log seed-extend dynamic programming problems
	logDpsOpp.clear();      // log mate-search dynamic programming problems
#ifdef USE_SRA
	sra_accs.clear();
#endif
}

static const char *short_options = "bfF:qbzhcu:rv:s:aP:t3:5:w:p:k:M:1:2:I:X:CQ:N:i:L:U:x:S:g:O:D:R:";

static struct option long_options[] = {
	{(char*)"verbose",                     no_argument,        0,                   ARG_VERBOSE},
	{(char*)"startverbose",                no_argument,        0,                   ARG_STARTVERBOSE},
	{(char*)"quiet",                       no_argument,        0,                   ARG_QUIET},
	{(char*)"sanity",                      no_argument,        0,                   ARG_SANITY},
	{(char*)"pause",                       no_argument,        &ipause,             1},
	{(char*)"orig",                        required_argument,  0,                   ARG_ORIG},
	{(char*)"all",                         no_argument,        0,                   'a'},
	{(char*)"solexa-quals",                no_argument,        0,                   ARG_SOLEXA_QUALS},
	{(char*)"integer-quals",               no_argument,        0,                   ARG_INTEGER_QUALS},
	{(char*)"int-quals",                   no_argument,        0,                   ARG_INTEGER_QUALS},
	{(char*)"metrics",                     required_argument,  0,                   ARG_METRIC_IVAL},
	{(char*)"metrics-file",                required_argument,  0,                   ARG_METRIC_FILE},
	{(char*)"metrics-stderr",              no_argument,        0,                   ARG_METRIC_STDERR},
	{(char*)"metrics-per-read",            no_argument,        0,                   ARG_METRIC_PER_READ},
	{(char*)"met-read",                    no_argument,        0,                   ARG_METRIC_PER_READ},
	{(char*)"met",                         required_argument,  0,                   ARG_METRIC_IVAL},
	{(char*)"met-file",                    required_argument,  0,                   ARG_METRIC_FILE},
	{(char*)"met-stderr",                  no_argument,        0,                   ARG_METRIC_STDERR},
	{(char*)"time",                        no_argument,        0,                   't'},
	{(char*)"trim3",                       required_argument,  0,                   '3'},
	{(char*)"trim5",                       required_argument,  0,                   '5'},
	{(char*)"seed",                        required_argument,  0,                   ARG_SEED},
	{(char*)"qupto",                       required_argument,  0,                   'u'},
	{(char*)"upto",                        required_argument,  0,                   'u'},
	{(char*)"version",                     no_argument,        0,                   ARG_VERSION},
	{(char*)"reads-per-batch",             required_argument,  0,                   ARG_READS_PER_BATCH},
	{(char*)"filepar",                     no_argument,        0,                   ARG_FILEPAR},
	{(char*)"help",                        no_argument,        0,                   'h'},
	{(char*)"threads",                     required_argument,  0,                   'p'},
	{(char*)"khits",                       required_argument,  0,                   'k'},
	{(char*)"minins",                      required_argument,  0,                   'I'},
	{(char*)"maxins",                      required_argument,  0,                   'X'},
	{(char*)"quals",                       required_argument,  0,                   'Q'},
	{(char*)"Q1",                          required_argument,  0,                   ARG_QUALS1},
	{(char*)"Q2",                          required_argument,  0,                   ARG_QUALS2},
	{(char*)"refidx",                      no_argument,        0,                   ARG_REFIDX},
	{(char*)"partition",                   required_argument,  0,                   ARG_PARTITION},
	{(char*)"ff",                          no_argument,        0,                   ARG_FF},
	{(char*)"fr",                          no_argument,        0,                   ARG_FR},
	{(char*)"rf",                          no_argument,        0,                   ARG_RF},
	{(char*)"cachelim",                    required_argument,  0,                   ARG_CACHE_LIM},
	{(char*)"cachesz",                     required_argument,  0,                   ARG_CACHE_SZ},
	{(char*)"nofw",                        no_argument,        0,                   ARG_NO_FW},
	{(char*)"norc",                        no_argument,        0,                   ARG_NO_RC},
	{(char*)"skip",                        required_argument,  0,                   's'},
	{(char*)"12",                          required_argument,  0,                   ARG_ONETWO},
	{(char*)"tab5",                        required_argument,  0,                   ARG_TAB5},
	{(char*)"tab6",                        required_argument,  0,                   ARG_TAB6},
	{(char*)"interleaved",                 required_argument,  0,                   ARG_INTERLEAVED},
	{(char*)"phred33-quals",               no_argument,        0,                   ARG_PHRED33},
	{(char*)"phred64-quals",               no_argument,        0,                   ARG_PHRED64},
	{(char*)"phred33",                     no_argument,        0,                   ARG_PHRED33},
	{(char*)"phred64",                     no_argument,        0,                   ARG_PHRED64},
	{(char*)"solexa1.3-quals",             no_argument,        0,                   ARG_PHRED64},
	{(char*)"mm",                          no_argument,        0,                   ARG_MM},
	{(char*)"shmem",                       no_argument,        0,                   ARG_SHMEM},
	{(char*)"mmsweep",                     no_argument,        0,                   ARG_MMSWEEP},
	{(char*)"hadoopout",                   no_argument,        0,                   ARG_HADOOPOUT},
	{(char*)"fullref",                     no_argument,        0,                   ARG_FULLREF},
	{(char*)"usage",                       no_argument,        0,                   ARG_USAGE},
	{(char*)"sam-no-qname-trunc",          no_argument,        0,                   ARG_SAM_NO_QNAME_TRUNC},
	{(char*)"sam-omit-sec-seq",            no_argument,        0,                   ARG_SAM_OMIT_SEC_SEQ},
	{(char*)"omit-sec-seq",                no_argument,        0,                   ARG_SAM_OMIT_SEC_SEQ},
	{(char*)"sam-no-head",                 no_argument,        0,                   ARG_SAM_NOHEAD},
	{(char*)"sam-nohead",                  no_argument,        0,                   ARG_SAM_NOHEAD},
	{(char*)"sam-noHD",                    no_argument,        0,                   ARG_SAM_NOHEAD},
	{(char*)"sam-no-hd",                   no_argument,        0,                   ARG_SAM_NOHEAD},
	{(char*)"sam-nosq",                    no_argument,        0,                   ARG_SAM_NOSQ},
	{(char*)"sam-no-sq",                   no_argument,        0,                   ARG_SAM_NOSQ},
	{(char*)"sam-noSQ",                    no_argument,        0,                   ARG_SAM_NOSQ},
	{(char*)"no-head",                     no_argument,        0,                   ARG_SAM_NOHEAD},
	{(char*)"no-hd",                       no_argument,        0,                   ARG_SAM_NOHEAD},
	{(char*)"no-sq",                       no_argument,        0,                   ARG_SAM_NOSQ},
	{(char*)"no-HD",                       no_argument,        0,                   ARG_SAM_NOHEAD},
	{(char*)"no-SQ",                       no_argument,        0,                   ARG_SAM_NOSQ},
	{(char*)"no-unal",                     no_argument,        0,                   ARG_SAM_NO_UNAL},
	{(char*)"sam-RG",                      required_argument,  0,                   ARG_SAM_RG},
	{(char*)"sam-rg",                      required_argument,  0,                   ARG_SAM_RG},
	{(char*)"sam-rg-id",                   required_argument,  0,                   ARG_SAM_RGID},
	{(char*)"RG",                          required_argument,  0,                   ARG_SAM_RG},
	{(char*)"rg",                          required_argument,  0,                   ARG_SAM_RG},
	{(char*)"rg-id",                       required_argument,  0,                   ARG_SAM_RGID},
	{(char*)"snpphred",                    required_argument,  0,                   ARG_SNPPHRED},
	{(char*)"snpfrac",                     required_argument,  0,                   ARG_SNPFRAC},
	{(char*)"gbar",                        required_argument,  0,                   ARG_GAP_BAR},
	{(char*)"qseq",                        no_argument,        0,                   ARG_QSEQ},
	{(char*)"policy",                      required_argument,  0,                   ARG_ALIGN_POLICY},
	{(char*)"preset",                      required_argument,  0,                   'P'},
	{(char*)"seed-summ",                   no_argument,        0,                   ARG_SEED_SUMM},
	{(char*)"seed-summary",                no_argument,        0,                   ARG_SEED_SUMM},
	{(char*)"overhang",                    no_argument,        0,                   ARG_OVERHANG},
	{(char*)"no-cache",                    no_argument,        0,                   ARG_NO_CACHE},
	{(char*)"cache",                       no_argument,        0,                   ARG_USE_CACHE},
	{(char*)"454",                         no_argument,        0,                   ARG_NOISY_HPOLY},
	{(char*)"ion-torrent",                 no_argument,        0,                   ARG_NOISY_HPOLY},
	{(char*)"no-mixed",                    no_argument,        0,                   ARG_NO_MIXED},
	{(char*)"no-discordant",               no_argument,        0,                   ARG_NO_DISCORDANT},
	{(char*)"local",                       no_argument,        0,                   ARG_LOCAL},
	{(char*)"end-to-end",                  no_argument,        0,                   ARG_END_TO_END},
	{(char*)"ungapped",                    no_argument,        0,                   ARG_UNGAPPED},
	{(char*)"no-ungapped",                 no_argument,        0,                   ARG_UNGAPPED_NO},
	{(char*)"sse8",                        no_argument,        0,                   ARG_SSE8},
	{(char*)"no-sse8",                     no_argument,        0,                   ARG_SSE8_NO},
	{(char*)"scan-narrowed",               no_argument,        0,                   ARG_SCAN_NARROWED},
	{(char*)"qc-filter",                   no_argument,        0,                   ARG_QC_FILTER},
	{(char*)"bwa-sw-like",                 no_argument,        0,                   ARG_BWA_SW_LIKE},
	{(char*)"multiseed",                   required_argument,  0,                   ARG_MULTISEED_IVAL},
	{(char*)"ma",                          required_argument,  0,                   ARG_SCORE_MA},
	{(char*)"mp",                          required_argument,  0,                   ARG_SCORE_MMP},
	{(char*)"np",                          required_argument,  0,                   ARG_SCORE_NP},
	{(char*)"rdg",                         required_argument,  0,                   ARG_SCORE_RDG},
	{(char*)"rfg",                         required_argument,  0,                   ARG_SCORE_RFG},
	{(char*)"score-min",                   required_argument,  0,                   ARG_SCORE_MIN},
	{(char*)"min-score",                   required_argument,  0,                   ARG_SCORE_MIN},
	{(char*)"n-ceil",                      required_argument,  0,                   ARG_N_CEIL},
	{(char*)"dpad",                        required_argument,  0,                   ARG_DPAD},
	{(char*)"mapq-print-inputs",           no_argument,        0,                   ARG_SAM_PRINT_YI},
	{(char*)"very-fast",                   no_argument,        0,                   ARG_PRESET_VERY_FAST},
	{(char*)"fast",                        no_argument,        0,                   ARG_PRESET_FAST},
	{(char*)"sensitive",                   no_argument,        0,                   ARG_PRESET_SENSITIVE},
	{(char*)"very-sensitive",              no_argument,        0,                   ARG_PRESET_VERY_SENSITIVE},
	{(char*)"very-fast-local",             no_argument,        0,                   ARG_PRESET_VERY_FAST_LOCAL},
	{(char*)"fast-local",                  no_argument,        0,                   ARG_PRESET_FAST_LOCAL},
	{(char*)"sensitive-local",             no_argument,        0,                   ARG_PRESET_SENSITIVE_LOCAL},
	{(char*)"very-sensitive-local",        no_argument,        0,                   ARG_PRESET_VERY_SENSITIVE_LOCAL},
	{(char*)"seedlen",                     required_argument,  0,                   'L'},
	{(char*)"seedmms",                     required_argument,  0,                   'N'},
	{(char*)"seedival",                    required_argument,  0,                   'i'},
	{(char*)"ignore-quals",                no_argument,        0,                   ARG_IGNORE_QUALS},
	{(char*)"index",                       required_argument,  0,                   'x'},
	{(char*)"arg-desc",                    no_argument,        0,                   ARG_DESC},
	{(char*)"wrapper",                     required_argument,  0,                   ARG_WRAPPER},
	{(char*)"unpaired",                    required_argument,  0,                   'U'},
	{(char*)"output",                      required_argument,  0,                   'S'},
	{(char*)"mapq-v",                      required_argument,  0,                   ARG_MAPQ_V},
	{(char*)"dovetail",                    no_argument,        0,                   ARG_DOVETAIL},
	{(char*)"no-dovetail",                 no_argument,        0,                   ARG_NO_DOVETAIL},
	{(char*)"contain",                     no_argument,        0,                   ARG_CONTAIN},
	{(char*)"no-contain",                  no_argument,        0,                   ARG_NO_CONTAIN},
	{(char*)"overlap",                     no_argument,        0,                   ARG_OVERLAP},
	{(char*)"no-overlap",                  no_argument,        0,                   ARG_NO_OVERLAP},
	{(char*)"tighten",                     required_argument,  0,                   ARG_TIGHTEN},
	{(char*)"exact-upfront",               no_argument,        0,                   ARG_EXACT_UPFRONT},
	{(char*)"1mm-upfront",                 no_argument,        0,                   ARG_1MM_UPFRONT},
	{(char*)"no-exact-upfront",            no_argument,        0,                   ARG_EXACT_UPFRONT_NO},
	{(char*)"no-1mm-upfront",              no_argument,        0,                   ARG_1MM_UPFRONT_NO},
	{(char*)"1mm-minlen",                  required_argument,  0,                   ARG_1MM_MINLEN},
	{(char*)"seed-off",                    required_argument,  0,                   'O'},
	{(char*)"seed-boost",                  required_argument,  0,                   ARG_SEED_BOOST_THRESH},
	{(char*)"read-times",                  no_argument,        0,                   ARG_READ_TIMES},
	{(char*)"show-rand-seed",              no_argument,        0,                   ARG_SHOW_RAND_SEED},
	{(char*)"dp-fail-streak",              required_argument,  0,                   ARG_DP_FAIL_STREAK_THRESH},
	{(char*)"ee-fail-streak",              required_argument,  0,                   ARG_EE_FAIL_STREAK_THRESH},
	{(char*)"ug-fail-streak",              required_argument,  0,                   ARG_UG_FAIL_STREAK_THRESH},
	{(char*)"fail-streak",                 required_argument,  0,                   'D'},
	{(char*)"dp-fails",                    required_argument,  0,                   ARG_DP_FAIL_THRESH},
	{(char*)"ug-fails",                    required_argument,  0,                   ARG_UG_FAIL_THRESH},
	{(char*)"extends",                     required_argument,  0,                   ARG_EXTEND_ITERS},
	{(char*)"no-extend",                   no_argument,        0,                   ARG_NO_EXTEND},
	{(char*)"mapq-extra",                  no_argument,        0,                   ARG_MAPQ_EX},
	{(char*)"seed-rounds",                 required_argument,  0,                   'R'},
	{(char*)"reorder",                     no_argument,        0,                   ARG_REORDER},
	{(char*)"passthrough",                 no_argument,        0,                   ARG_READ_PASSTHRU},
	{(char*)"sample",                      required_argument,  0,                   ARG_SAMPLE},
	{(char*)"cp-min",                      required_argument,  0,                   ARG_CP_MIN},
	{(char*)"cp-ival",                     required_argument,  0,                   ARG_CP_IVAL},
	{(char*)"tri",                         no_argument,        0,                   ARG_TRI},
	{(char*)"nondeterministic",            no_argument,        0,                   ARG_NON_DETERMINISTIC},
	{(char*)"non-deterministic",           no_argument,        0,                   ARG_NON_DETERMINISTIC},
	{(char*)"local-seed-cache-sz",         required_argument,  0,                   ARG_LOCAL_SEED_CACHE_SZ},
	{(char*)"seed-cache-sz",               required_argument,  0,                   ARG_CURRENT_SEED_CACHE_SZ},
	{(char*)"no-unal",                     no_argument,        0,                   ARG_SAM_NO_UNAL},
	{(char*)"test-25",                     no_argument,        0,                   ARG_TEST_25},
// TODO: following should be a function of read length?
	{(char*)"desc-kb",                     required_argument,  0,                   ARG_DESC_KB},
	{(char*)"desc-landing",                required_argument,  0,                   ARG_DESC_LANDING},
	{(char*)"desc-exp",                    required_argument,  0,                   ARG_DESC_EXP},
	{(char*)"desc-prioritize",             no_argument,        0,                   ARG_DESC_PRIORITIZE},
	{(char*)"desc-fmops",                  required_argument,  0,                   ARG_DESC_FMOPS},
	{(char*)"log-dp",                      required_argument,  0,                   ARG_LOG_DP},
	{(char*)"log-dp-opp",                  required_argument,  0,                   ARG_LOG_DP_OPP},
	{(char*)"soft-clipped-unmapped-tlen",  no_argument,        0,                   ARG_SC_UNMAPPED},
	{(char*)"xeq",                         no_argument,        0,                   ARG_XEQ},
	{(char*)"thread-ceiling",              required_argument,  0,                   ARG_THREAD_CEILING},
	{(char*)"thread-piddir",               required_argument,  0,                   ARG_THREAD_PIDDIR},
	{(char*)"trim-to",                     required_argument,  0,                   ARG_TRIM_TO},
	{(char*)"preserve-tags",               no_argument,        0,                   ARG_PRESERVE_TAGS},
	{(char*)"align-paired-reads",          no_argument,        0,                   ARG_ALIGN_PAIRED_READS},
#ifdef USE_SRA
	{(char*)"sra-acc",                     required_argument,  0,                   ARG_SRA_ACC},
#endif
	{(char*)"sam-append-comment",          no_argument,        0,                   ARG_SAM_APPEND_COMMENT},
	{(char*)0,                             0,                  0,                   0} //  terminator
};

/**
 * Print out a concise description of what options are taken and whether they
 * take an argument.
 */
static void printArgDesc(ostream& out) {
	// struct option {
	//   const char *name;
	//   int has_arg;
	//   int *flag;
	//   int val;
	// };
	size_t i = 0;
	while(long_options[i].name != 0) {
		out << long_options[i].name << "\t"
		    << (long_options[i].has_arg == no_argument ? 0 : 1)
		    << endl;
		i++;
	}
	size_t solen = strlen(short_options);
	for(i = 0; i < solen; i++) {
		// Has an option?  Does if next char is :
		if(i == solen-1) {
			assert_neq(':', short_options[i]);
			cout << (char)short_options[i] << "\t" << 0 << endl;
		} else {
			if(short_options[i+1] == ':') {
				// Option with argument
				cout << (char)short_options[i] << "\t" << 1 << endl;
				i++; // skip the ':'
			} else {
				// Option with no argument
				cout << (char)short_options[i] << "\t" << 0 << endl;
			}
		}
	}
}

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Bowtie 2 version " << string(BOWTIE2_VERSION).c_str() << " by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)" << endl;
	string tool_name = "bowtie2-align";
	if(wrapper == "basic-0") {
		tool_name = "bowtie2";
	}
	out << "Usage: " << endl
#ifdef USE_SRA
	    << "  " << tool_name.c_str() << " [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc> | -b <bam>} [-S <sam>]" << endl
#else
	    << "  " << tool_name.c_str() << " [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]" << endl
#endif
	    << endl
	    <<     "  <bt2-idx>  Index filename prefix (minus trailing .X." + gEbwt_ext + ")." << endl
	    <<     "             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible." << endl
	    <<     "  <m1>       Files with #1 mates, paired with files in <m2>." << endl;
	if(wrapper == "basic-0") {
		out << "             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." << endl;
	}
	out <<     "  <m2>       Files with #2 mates, paired with files in <m1>." << endl;
	if(wrapper == "basic-0") {
		out << "             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." << endl;
	}
	out <<     "  <r>        Files with unpaired reads." << endl;
	if(wrapper == "basic-0") {
		out << "             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." << endl;
	}
	out <<     "  <i>        Files with interleaved paired-end FASTQ/FASTA reads" << endl;
	if(wrapper == "basic-0") {
		out << "             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." << endl;
	}
#ifdef USE_SRA
	out <<     "  <acc>      Files are SRA accessions. Accessions not found in local storage will\n"
	    <<     "             be fetched from NCBI." << endl;
#endif
	out <<     "  <bam>      Files are unaligned BAM sorted by read name." << endl;
	out <<     "  <sam>      File for SAM output (default: stdout)" << endl
	    << endl
	    << "  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be" << endl
	    << "  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'." << endl
		// Wrapper script should write <bam> line next
	    << endl
	    << "Options (defaults in parentheses):" << endl
	    << endl
	    << " Input:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  --tab5             query input files are TAB5 .tab5" << endl
	    << "  --tab6             query input files are TAB6 .tab6" << endl
	    << "  --qseq             query input files are in Illumina's qseq format" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    << "  -r                 query input files are raw one-sequence-per-line" << endl
	    << "  -F k:<int>,i:<int> query input files are continuous FASTA where reads" << endl
	    << "                     are substrings (k-mers) extracted from a FASTA file <s>" << endl
	    << "                     and aligned at offsets 1, 1+i, 1+2i ... end of reference" << endl
	    << "  -c                 <m1>, <m2>, <r> are sequences themselves, not files" << endl
	    << "  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)" << endl
	    << "  -u/--upto <int>    stop after first <int> reads/pairs (no limit)" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)" << endl
	    << "  --trim-to [3:|5:]<int> trim reads exceeding <int> bases from either 3' or 5' end" << endl
	    << "                     If the read end is not specified then it defaults to 3 (0)" << endl
	    << "  --phred33          qualities are Phred+33 (default)" << endl
	    << "  --phred64          qualities are Phred+64" << endl
	    << "  --int-quals        qualities encoded as space-delimited integers" << endl
	    << endl
	    << " Presets:                 Same as:" << endl
	    << "  For --end-to-end:" << endl
	    << "   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50" << endl
	    << "   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50" << endl
	    << "   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)" << endl
	    << "   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50" << endl
	    << endl
	    << "  For --local:" << endl
	    << "   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00" << endl
	    << "   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75" << endl
	    << "   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)" << endl
	    << "   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50" << endl
	    << endl
	    << " Alignment:" << endl
	    << "  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)" << endl
	    << "  -L <int>           length of seed substrings; must be >3, <32 (22)" << endl
	    << "  -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)" << endl
	    << "  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)" << endl
	    << "  --dpad <int>       include <int> extra ref chars on sides of DP table (15)" << endl
	    << "  --gbar <int>       disallow gaps within <int> nucs of read extremes (4)" << endl
	    << "  --ignore-quals     treat all quality values as 30 on Phred scale (off)" << endl
	    << "  --nofw             do not align forward (original) version of read (off)" << endl
	    << "  --norc             do not align reverse-complement version of read (off)" << endl
	    << "  --no-1mm-upfront   do not allow 1 mismatch alignments before attempting to" << endl
	    << "                     scan for the optimal seeded alignments"
	    << endl
	    << "  --end-to-end       entire read must align; no clipping (on)" << endl
	    << "   OR" << endl
	    << "  --local            local alignment; ends might be soft clipped (off)" << endl
	    << endl
	    << " Scoring:" << endl
	    << "  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) " << endl
	    << "  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)" << endl
	    << "  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)" << endl
	    << "  --rdg <int>,<int>  read gap open, extend penalties (5,3)" << endl
	    << "  --rfg <int>,<int>  reference gap open, extend penalties (5,3)" << endl
	    << "  --score-min <func> min acceptable alignment score w/r/t read length" << endl
	    << "                     (G,20,8 for local, L,-0.6,-0.6 for end-to-end)" << endl
	    << endl
	    << " Reporting:" << endl
	    << "  (default)          look for multiple alignments, report best, with MAPQ" << endl
	    << "   OR" << endl
	    << "  -k <int>           report up to <int> alns per read; MAPQ not meaningful" << endl
	    << "   OR" << endl
	    << "  -a/--all           report all alignments; very slow, MAPQ not meaningful" << endl
	    << endl
	    << " Effort:" << endl
	    << "  -D <int>           give up extending after <int> failed extends in a row (15)" << endl
	    << "  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)" << endl
	    << endl
	    << " Paired-end:" << endl
	    << "  -I/--minins <int>  minimum fragment length (0)" << endl
	    << "  -X/--maxins <int>  maximum fragment length (500)" << endl
	    << "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)" << endl
	    << "  --no-mixed         suppress unpaired alignments for paired reads" << endl
	    << "  --no-discordant    suppress discordant alignments for paired reads" << endl
	    << "  --dovetail         concordant when mates extend past each other" << endl
	    << "  --no-contain       not concordant when one mate alignment contains other" << endl
	    << "  --no-overlap       not concordant when mates overlap at all" << endl
	    << endl
	    << " BAM:" << endl
	    << "  --align-paired-reads" << endl
	    << "                     Bowtie2 will, by default, attempt to align unpaired BAM reads." << endl
	    << "                     Use this option to align paired-end reads instead." << endl
	    << "  --preserve-tags    Preserve tags from the original BAM record by" << endl
	    << "                     appending them to the end of the corresponding SAM output." << endl
	    << endl
	    << " Output:" << endl;
	//if(wrapper == "basic-0") {
	//	out << "  --bam              output directly to BAM (by piping through 'samtools view')" << endl;
	//}
	out << "  -t/--time          print wall-clock time taken by search phases" << endl;
	if(wrapper == "basic-0") {
		out << "  --un <path>        write unpaired reads that didn't align to <path>" << endl
		    << "  --al <path>        write unpaired reads that aligned at least once to <path>" << endl
		    << "  --un-conc <path>   write pairs that didn't align concordantly to <path>" << endl
		    << "  --al-conc <path>   write pairs that aligned concordantly at least once to <path>" << endl
		    << "    (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g." << endl
		    << "    --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)" << endl;
	}
	out << "  --quiet            print nothing to stderr except serious errors" << endl
		//  << "  --refidx           refer to ref. seqs by 0-based index rather than name" << endl
	    << "  --met-file <path>  send metrics to file at <path> (off)" << endl
	    << "  --met-stderr       send metrics to stderr (off)" << endl
	    << "  --met <int>        report internal counters & metrics every <int> secs (1)" << endl
		// Following is supported in the wrapper instead
	    << "  --no-unal          suppress SAM records for unaligned reads" << endl
	    << "  --no-head          suppress header lines, i.e. lines starting with @" << endl
	    << "  --no-sq            suppress @SQ header lines" << endl
	    << "  --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field" << endl
	    << "  --rg <text>        add <text> (\"lab:value\") to @RG line of SAM header." << endl
	    << "                     Note: @RG line only printed when --rg-id is set." << endl
	    << "  --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments." << endl
	    << "  --sam-no-qname-trunc" << endl
	    << "                     Suppress standard behavior of truncating readname at first whitespace " << endl
	    << "                     at the expense of generating non-standard SAM." << endl
	    << "  --xeq              Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record." << endl
	    << "  --soft-clipped-unmapped-tlen" << endl
	    << "                     Exclude soft-clipped bases when reporting TLEN" << endl
	    << "  --sam-append-comment" << endl
	    << "                     Append FASTA/FASTQ comment to SAM record" << endl
	    << endl
	    << " Performance:" << endl
		//    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
	    << "  -p/--threads <int> number of alignment threads to launch (1)" << endl
	    << "  --reorder          force SAM output order to match order of input reads" << endl
#ifdef BOWTIE_MM
	    << "  --mm               use memory-mapped I/O for index; many 'bowtie's can share" << endl
#endif
#ifdef BOWTIE_SHARED_MEM
		//<< "  --shmem            use shared mem for index; many 'bowtie's can share" << endl
#endif
	    << endl
	    << " Other:" << endl
	    << "  --qc-filter        filter out reads that are bad according to QSEQ filter" << endl
	    << "  --seed <int>       seed for random number generator (0)" << endl
	    << "  --non-deterministic" << endl
	    << "                     seed rand. gen. arbitrarily instead of using read attributes" << endl
		//  << "  --verbose          verbose output for debugging" << endl
	    << "  --version          print version information and quit" << endl
	    << "  -h/--help          print this usage message" << endl;
	if(wrapper.empty()) {
		cerr << endl
		     << "*** Warning ***" << endl
		     << "'bowtie2-align' was run directly.  It is recommended that you run the wrapper script 'bowtie2' instead." << endl
		     << endl;
	}
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, int upper, const char *errmsg, const char *arg) {
	long l;
	char *endPtr= NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower || l > upper) {
			cerr << errmsg << endl;
			printUsage(cerr);
			throw 1;
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

/**
 * Upper is maximum int by default.
 */
static int parseInt(int lower, const char *errmsg, const char *arg) {
	return parseInt(lower, std::numeric_limits<int>::max(), errmsg, arg);
}

/**
 * Parse a T string 'str'.
 */
template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

/**
 * Parse a pair of Ts from a string, 'str', delimited with 'delim'.
 */
template<typename T>
pair<T, T> parsePair(const char *str, char delim) {
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	pair<T, T> ret;
	ret.first = parse<T>(ss[0].c_str());
	ret.second = parse<T>(ss[1].c_str());
	return ret;
}

/**
 * Parse a pair of Ts from a string, 'str', delimited with 'delim'.
 */
template<typename T>
void parseTuple(const char *str, char delim, EList<T>& ret) {
	string s(str);
	EList<string> ss;
	tokenize(s, delim, ss);
	for(size_t i = 0; i < ss.size(); i++) {
		ret.push_back(parse<T>(ss[i].c_str()));
	}
}

static string applyPreset(const string& sorig, Presets& presets) {
	string s = sorig;
	size_t found = s.find("%LOCAL%");
	if(found != string::npos) {
		s.replace(found, strlen("%LOCAL%"), localAlign ? "-local" : "");
	}
	if(gVerbose) {
		cerr << "Applying preset: '" << s.c_str() << "' using preset menu '"
		     << presets.name() << "'" << endl;
	}
	string pol;
	presets.apply(s, pol, extra_opts);
	return pol;
}

static bool saw_M;
static bool saw_a;
static bool saw_k;
static bool saw_trim3;
static bool saw_trim5;
static bool saw_trim_to;
static bool saw_bam;
static bool saw_preserve_tags;
static bool saw_align_paired_reads;
static EList<string> presetList;

/**
 * TODO: Argument parsing is very, very flawed.  The biggest problem is that
 * there are two separate worlds of arguments, the ones set via polstr, and
 * the ones set directly in variables.  This makes for nasty interactions,
 * e.g., with the -M option being resolved at an awkward time relative to
 * the -k and -a options.
 */
static void parseOption(int next_option, const char *arg) {
	switch (next_option) {
	case ARG_TEST_25: bowtie2p5 = true; break;
	case ARG_DESC_KB: descentTotSz = SimpleFunc::parse(arg, 0.0, 1024.0, 1024.0, DMAX); break;
	case ARG_DESC_FMOPS: descentTotFmops = SimpleFunc::parse(arg, 0.0, 10.0, 100.0, DMAX); break;
	case ARG_LOG_DP: logDps = arg; break;  // NOTE: Deprecated, noop
	case ARG_LOG_DP_OPP: logDpsOpp = arg; break; // NOTE: Deprecated, noop
	case ARG_DESC_LANDING: {
		descLanding = parse<int>(arg);
		if(descLanding < 1) {
			cerr << "Error: --desc-landing must be greater than or equal to 1" << endl;
			throw 1;
		}
		break;
	}
	case ARG_DESC_EXP: {
		descConsExp = parse<double>(arg);
		if(descConsExp < 0.0) {
			cerr << "Error: --desc-exp must be greater than or equal to 0" << endl;
			throw 1;
		}
		break;
	}
	case ARG_DESC_PRIORITIZE: descPrioritizeRoots = true; break;
	case '1': tokenize(arg, ",", mates1); break;
	case '2': tokenize(arg, ",", mates2); break;
	case ARG_ONETWO: tokenize(arg, ",", mates12); set_format(format, TAB_MATE5); break;
	case ARG_TAB5:   tokenize(arg, ",", mates12); set_format(format, TAB_MATE5); break;
	case ARG_TAB6:   tokenize(arg, ",", mates12); set_format(format, TAB_MATE6); break;
	case ARG_INTERLEAVED: {
		tokenize(arg, ",", mates12);
		interleaved = true;
		break;
	}
	case 'b': {
		set_format(format, BAM);
		saw_bam = true;
		break;
	}
	case 'f': set_format(format, FASTA); break;
	case 'F': {
		set_format(format, FASTA_CONT);
		pair<uint32_t, uint32_t> p = parsePair<uint32_t>(arg, ',');
		fastaContLen = p.first;
		fastaContFreq = p.second;
		break;
	}
	case ARG_BWA_SW_LIKE: {
		cerr << "WARNING: BWA_SW_LIKE not supported" << endl; 
		break;
	}
	case 'q': set_format(format, FASTQ); break;
	case 'r': set_format(format, RAW); break;
	case 'c': set_format(format, CMDLINE); break;
	case ARG_QSEQ: set_format(format, QSEQ); break;
	case 'I':
		gMinInsert = parseInt(0, "-I arg must be positive", arg);
		break;
	case 'X':
		gMaxInsert = parseInt(1, "-X arg must be at least 1", arg);
		break;
	case ARG_NO_DISCORDANT: gReportDiscordant = false; break;
	case ARG_NO_MIXED: gReportMixed = false; break;
	case 's':
		cerr << "WARNING: skipReads not supported" << endl; 
		break;
	case ARG_FF: gMate1fw = true;  gMate2fw = true;  break;
	case ARG_RF: gMate1fw = false; gMate2fw = true;  break;
	case ARG_FR: gMate1fw = true;  gMate2fw = false; break;
	case ARG_SHMEM: useShmem = true; break;
	case ARG_SEED_SUMM:
		cerr << "WARNING: seedSumm not supported" << endl; 
		break;
	case ARG_SC_UNMAPPED: scUnMapped = true; break;
	case ARG_XEQ: xeq = true; break;
	case ARG_PRESERVE_TAGS: {
		preserve_tags = true;
		saw_preserve_tags = true;
		break;
	}
	case ARG_ALIGN_PAIRED_READS: {
		align_paired_reads = true;
		saw_align_paired_reads = true;
		break;
	}
	case ARG_MM: {
#ifdef BOWTIE_MM
		useMm = true;
		break;
#else
		cerr << "Memory-mapped I/O mode is disabled because bowtie was not compiled with" << endl
		     << "BOWTIE_MM defined.  Memory-mapped I/O is not supported under Windows.  If you" << endl
		     << "would like to use memory-mapped I/O on a platform that supports it, please" << endl
		     << "refrain from specifying BOWTIE_MM=0 when compiling Bowtie." << endl;
		throw 1;
#endif
	}
	case ARG_MMSWEEP: mmSweep = true; break;
	case ARG_HADOOPOUT: hadoopOut = true; break;
	case ARG_SOLEXA_QUALS: solexaQuals = true; break;
	case ARG_INTEGER_QUALS: integerQuals = true; break;
	case ARG_PHRED64: phred64Quals = true; break;
	case ARG_PHRED33: solexaQuals = false; phred64Quals = false; break;
	case ARG_OVERHANG: gReportOverhangs = true; break;
	case ARG_NO_CACHE: msNoCache = true; break;  // default
	case ARG_USE_CACHE: 
		cerr << "WARNING: USE_CACHE not supported" << endl; 
		break;
	case ARG_LOCAL_SEED_CACHE_SZ:
		seedCacheLocalMB = (uint32_t)parseInt(1, "--local-seed-cache-sz arg must be at least 1", arg);
		break;
	case ARG_CURRENT_SEED_CACHE_SZ:
		seedCacheCurrentMB = (uint32_t)parseInt(1, "--seed-cache-sz arg must be at least 1", arg);
		break;
	case ARG_REFIDX: noRefNames = true; break;
	case ARG_FULLREF: fullRef = true; break;
	case ARG_GAP_BAR:
		gGapBarrier = parseInt(1, "--gbar must be no less than 1", arg);
		break;
	case ARG_SEED:
		seed = parseInt(0, "--seed arg must be at least 0", arg);
		break;
	case ARG_NON_DETERMINISTIC:
		cerr << "WARNING: arbitraryRandom not supported" << endl; 
		break;
	case 'u':
		cerr << "WARNING: qupto not supported" << endl; 
		break;
	case 'Q':
		tokenize(arg, ",", qualities);
		integerQuals = true;
		break;
	case ARG_QUALS1:
		tokenize(arg, ",", qualities1);
		integerQuals = true;
		break;
	case ARG_QUALS2:
		tokenize(arg, ",", qualities2);
		integerQuals = true;
		break;
	case ARG_CACHE_LIM:
		cacheLimit = (uint32_t)parseInt(1, "--cachelim arg must be at least 1", arg);
		break;
	case ARG_CACHE_SZ:
		cacheSize = (uint32_t)parseInt(1, "--cachesz arg must be at least 1", arg);
		cacheSize *= (1024 * 1024); // convert from MB to B
		break;
	case ARG_WRAPPER: wrapper = arg; break;
	case 'p':
		nthreads = parseInt(1, "-p/--threads arg must be at least 1", arg);
		break;
	case ARG_THREAD_CEILING:
		thread_ceiling = parseInt(0, "--thread-ceiling must be at least 0", arg);
		break;
	case ARG_THREAD_PIDDIR:
		cerr << "WARNING: THREAD_PIDDIR not supported" << endl; 
		break;
	case ARG_FILEPAR:
		fileParallel = true;
		break;
	case '3': gTrim3 = parseInt(0, "-3/--trim3 arg must be at least 0", arg); break;
	case '5': gTrim5 = parseInt(0, "-5/--trim5 arg must be at least 0", arg); break;
	case ARG_TRIM_TO: {
		if (strlen(arg) > 1 && arg[1] != ':') {
			trimTo.first = 3;
			trimTo.second = parseInt(0, "--trim-to: the number of bases to trim must be at least 0", arg);
			break;
		}
		pair<int, int> res = parsePair<int>(arg, ':');
		if (res.first != 3 && res.first != 5) {
			cerr << "--trim-to: trim position must be either 3 or 5" << endl;
			printUsage(cerr);
			throw 1;
		}
		if(res.second < 0) {
			cerr << "--trim-to: the number bases to trim must be at least 0" << endl;
			printUsage(cerr);
			throw 1;
		}
		trimTo = static_cast<pair<short, size_t> >(res);
		break;
	}
	case 'h': printUsage(cout); throw 0; break;
	case ARG_USAGE: printUsage(cout); throw 0; break;
		//
		// NOTE that unlike in Bowtie 1, -M, -a and -k are mutually
		// exclusive here.
		//
	case 'M': {
		msample = true;
		mhits = parse<uint32_t>(arg);
		if(saw_a || saw_k) {
			cerr << "Warning: -M, -k and -a are mutually exclusive. "
			     << "-M will override" << endl;
			khits = 1;
		}
		assert_eq(1, khits);
		saw_M = true;
		cerr << "Warning: -M is deprecated.  Use -D and -R to adjust " <<
			"effort instead." << endl;
		break;
	}
	case ARG_EXTEND_ITERS: {
		maxIters = parse<size_t>(arg);
		break;
	}
	case ARG_NO_EXTEND: {
		doExtend = false;
		break;
	}
	case 'R': { polstr += ";ROUNDS="; polstr += arg; break; }
	case 'D': { polstr += ";DPS=";    polstr += arg; break; }
	case ARG_DP_MATE_STREAK_THRESH: {
		maxMateStreak = parse<size_t>(arg);
		break;
	}
	case ARG_DP_FAIL_STREAK_THRESH: {
		maxDpStreak = parse<size_t>(arg);
		break;
	}
	case ARG_EE_FAIL_STREAK_THRESH: {
		maxEeStreak = parse<size_t>(arg);
		break;
	}
	case ARG_UG_FAIL_STREAK_THRESH: {
		maxUgStreak = parse<size_t>(arg);
		break;
	}
	case ARG_DP_FAIL_THRESH: {
		maxDp = parse<size_t>(arg);
		break;
	}
	case ARG_UG_FAIL_THRESH: {
		maxUg = parse<size_t>(arg);
		break;
	}
	case ARG_SEED_BOOST_THRESH: {
		seedBoostThresh = parse<size_t>(arg);
		break;
	}
	case 'a': {
		cerr << "WARNING: allHits not supported" << endl;
		//msample = false;
		//allHits = true;
		//mhits = 0; // disable -M
		//if(saw_M || saw_k) {
		//	cerr << "Warning: -M, -k and -a are mutually exclusive. "
		//	     << "-a will override" << endl;
		//}
		//saw_a = true;
		break;
	}
	case 'k': {
		msample = false;
		khits = (uint32_t)parseInt(1, "-k arg must be at least 1", arg);
		mhits = 0; // disable -M
		if(saw_M || saw_a) {
			cerr << "Warning: -M, -k and -a are mutually exclusive. "
			     << "-k will override" << endl;
		}
		saw_k = true;
		break;
	}
	case ARG_VERBOSE: gVerbose = 1; break;
	case ARG_STARTVERBOSE: startVerbose = true; break;
	case ARG_QUIET: gQuiet = true; break;
	case ARG_SANITY: sanityCheck = true; break;
	case 't': timing = true; break;
	case ARG_METRIC_IVAL: {
		cerr << "WARNING: metricsIval not supported" << endl;
		break;
	}
	case ARG_METRIC_FILE: {
		cerr << "WARNING: metricsFile not supported" << endl;
		break;
        }
	case ARG_METRIC_STDERR: {
		cerr << "WARNING: metricsStderr not supported" << endl;
		break;
        }
	case ARG_METRIC_PER_READ: {
		cerr << "WARNING: metricsPerRead not supported" << endl; 
		break;
        }
	case ARG_NO_FW: gNofw = true; break;
	case ARG_NO_RC: gNorc = true; break;
	case ARG_SAM_NO_QNAME_TRUNC: samTruncQname = false; break;
	case ARG_SAM_APPEND_COMMENT: samAppendComment = true; break;
	case ARG_SAM_OMIT_SEC_SEQ: samOmitSecSeqQual = true; break;
	case ARG_SAM_NO_UNAL: samNoUnal = true; break;
	case ARG_SAM_NOHEAD: samNoHead = true; break;
	case ARG_SAM_NOSQ: samNoSQ = true; break;
	case ARG_SAM_PRINT_YI: sam_print_yi = true; break;
	case ARG_REORDER: reorder = true; break;
	case ARG_MAPQ_EX: {
		sam_print_zt = true;
		break;
	}
	case ARG_SHOW_RAND_SEED: {
		sam_print_zs = true;
		break;
	}
	case ARG_SAMPLE:
		cerr << "WARNING: sampleFrac not supported" << endl; 
		break;
	case ARG_CP_MIN:
		cminlen = parse<size_t>(arg);
		break;
	case ARG_CP_IVAL:
		cpow2 = parse<size_t>(arg);
		break;
	case ARG_TRI:
		doTri = true;
		break;
	case ARG_READ_PASSTHRU: {
		sam_print_xr = true;
		break;
	}
	case ARG_READ_TIMES: {
		cerr << "WARNING: Read_Times not supported" << endl;
		break;
	}
	case ARG_SAM_RG: {
		string argstr = arg;
		if(argstr.substr(0, 3) == "ID:") {
			rgid = "\t";
			rgid += argstr;
			rgs_optflag = "RG:Z:" + argstr.substr(3);
		} else {
			rgs += '\t';
			rgs += argstr;
		}
		break;
	}
	case ARG_SAM_RGID: {
		string argstr = arg;
		rgid = "\t";
		rgid = "\tID:" + argstr;
		rgs_optflag = "RG:Z:" + argstr;
		break;
	}
	case ARG_PARTITION: partitionSz = parse<int>(arg); break;
	case ARG_READS_PER_BATCH:
		readsPerBatch = parseInt(1, "--reads-per-batch arg must be at least 1", arg);
		break;
	case ARG_DPAD:
		maxhalf = parseInt(0, "--dpad must be no less than 0", arg);
		break;
	case ARG_ORIG:
		if(arg == NULL || strlen(arg) == 0) {
			cerr << "--orig arg must be followed by a string" << endl;
			printUsage(cerr);
			throw 1;
		}
		origString = arg;
		break;
	case ARG_LOCAL: {
		cerr << "WARNING: localAlign not supported" << endl; 
		break;
	}
	case ARG_END_TO_END: localAlign = false; break;
	case ARG_SSE8: enable8 = true; break;
	case ARG_SSE8_NO: enable8 = false; break;
	case ARG_UNGAPPED:  // silently ignore
                break;
	case ARG_UNGAPPED_NO: // silently ignore
                break;
	case ARG_NO_DOVETAIL: gDovetailMatesOK = false; break;
	case ARG_NO_CONTAIN:  gContainMatesOK  = false; break;
	case ARG_NO_OVERLAP:  gOlapMatesOK     = false; break;
	case ARG_DOVETAIL:    gDovetailMatesOK = true;  break;
	case ARG_CONTAIN:     gContainMatesOK  = true;  break;
	case ARG_OVERLAP:     gOlapMatesOK     = true;  break;
	case ARG_QC_FILTER: qcFilter = true; break;
	case ARG_IGNORE_QUALS: ignoreQuals = true; break;
	case ARG_MAPQ_V: mapqv = parse<int>(arg); break;
	case ARG_TIGHTEN: tighten = parse<int>(arg); break;
	case ARG_EXACT_UPFRONT:    /* noop */ break;
	case ARG_1MM_UPFRONT:      /* noop */ break;
	case ARG_EXACT_UPFRONT_NO:
		cerr << "WARNING: No doExactUpFront not supported" << endl; 
		break;
	case ARG_1MM_UPFRONT_NO:
		cerr << "WARNING: do1mmUpFront not supported" << endl; 
		break;
	case ARG_1MM_MINLEN:       do1mmMinLen = parse<size_t>(arg); break;
	case ARG_NOISY_HPOLY: noisyHpolymer = true; break;
	case 'x': bt2index = arg; break;
	case ARG_PRESET_VERY_FAST_LOCAL:
		cerr << "WARNING: localAlign not supported" << endl; 
		break;
	case ARG_PRESET_VERY_FAST: {
		presetList.push_back("very-fast%LOCAL%"); break;
	}
	case ARG_PRESET_FAST_LOCAL:
		cerr << "WARNING: localAlign not supported" << endl; 
		break;
	case ARG_PRESET_FAST: {
		presetList.push_back("fast%LOCAL%"); break;
	}
	case ARG_PRESET_SENSITIVE_LOCAL:
		cerr << "WARNING: localAlign not supported" << endl; 
		break;
	case ARG_PRESET_SENSITIVE: {
		presetList.push_back("sensitive%LOCAL%"); break;
	}
	case ARG_PRESET_VERY_SENSITIVE_LOCAL:
		cerr << "WARNING: localAlign not supported" << endl; 
		break;
	case ARG_PRESET_VERY_SENSITIVE: {
		presetList.push_back("very-sensitive%LOCAL%"); break;
	}
	case 'P': { presetList.push_back(arg); break; }
	case ARG_ALIGN_POLICY: {
		if(strlen(arg) > 0) {
			polstr += ";"; polstr += arg;
		}
		break;
	}
	case 'N': {
		int64_t len = parse<size_t>(arg);
		if (len < 0 || len > 1) {
			cerr << "Error: -N argument must be within the interval [0,1]; was " << arg << endl;
			throw 1;
		}
		polstr += ";SEED=";
		polstr += arg;
		break;
	}
	case 'L': {
		int64_t len = parse<size_t>(arg);
		if(len < 1 || len > 32) {
			cerr << "Error: -L argument must be within the interval [1,32]; was " << arg << endl;
			throw 1;
		}
		polstr += ";SEEDLEN=";
		polstr += arg;
		break;
	}
	case 'O':
		multiseedOff = parse<size_t>(arg);
		break;
	case 'i': {
		EList<string> args;
		tokenize(arg, ",", args);
		if(args.size() > 3 || args.size() == 0) {
			cerr << "Error: expected 3 or fewer comma-separated "
			     << "arguments to -i option, got "
			     << args.size() << endl;
			throw 1;
		}
		// Interval-settings arguments
		polstr += (";IVAL=" + args[0]); // Function type
		if(args.size() > 1) {
			polstr += ("," + args[1]);  // Constant term
		}
		if(args.size() > 2) {
			polstr += ("," + args[2]);  // Coefficient
		}
		break;
	}
	case ARG_MULTISEED_IVAL: {
		polstr += ";";
		// Split argument by comma
		EList<string> args;
		tokenize(arg, ",", args);
		if(args.size() > 5 || args.size() == 0) {
			cerr << "Error: expected 5 or fewer comma-separated "
			     << "arguments to --multiseed option, got "
			     << args.size() << endl;
			throw 1;
		}
		// Seed mm and length arguments
		polstr += "SEED=";
		polstr += (args[0]); // # mismatches
		if(args.size() >  1) polstr += (";SEEDLEN=" + args[1]); // length
		if(args.size() >  2) polstr += (";IVAL=" + args[2]); // Func type
		if(args.size() >  3) polstr += ("," + args[ 3]); // Constant term
		if(args.size() >  4) polstr += ("," + args[ 4]); // Coefficient
		break;
	}
	case ARG_N_CEIL: {
		// Split argument by comma
		EList<string> args;
		tokenize(arg, ",", args);
		if(args.size() > 3) {
			cerr << "Error: expected 3 or fewer comma-separated "
			     << "arguments to --n-ceil option, got "
			     << args.size() << endl;
			throw 1;
		}
		if(args.size() == 0) {
			cerr << "Error: expected at least one argument to --n-ceil option" << endl;
			throw 1;
		}
		polstr += ";NCEIL=";
		if(args.size() == 3) {
			polstr += (args[0] + "," + args[1] + "," + args[2]);
		} else {
			polstr += ("L," + args[0]);
			if(args.size() > 1) {
				polstr += ("," + (args[1]));
			}
		}
		break;
	}
	case ARG_SCORE_MA:  polstr += ";MA=";    polstr += arg; break;
	case ARG_SCORE_MMP: {
		EList<string> args;
		tokenize(arg, ",", args);
		if(args.size() > 2 || args.size() == 0) {
			cerr << "Error: expected 1 or 2 comma-separated "
			     << "arguments to --mmp option, got " << args.size() << endl;
			throw 1;
		}
		if(args.size() >= 1) {
			polstr += ";MMP=Q,";
			polstr += args[0];
			if(args.size() >= 2) {
				polstr += ",";
				polstr += args[1];
			}
		}
		break;
	}
	case ARG_SCORE_NP:  polstr += ";NP=C";   polstr += arg; break;
	case ARG_SCORE_RDG: polstr += ";RDG=";   polstr += arg; break;
	case ARG_SCORE_RFG: polstr += ";RFG=";   polstr += arg; break;
	case ARG_SCORE_MIN: {
		polstr += ";";
		EList<string> args;
		tokenize(arg, ",", args);
		if(args.size() > 3 || args.size() == 0) {
			cerr << "Error: expected 3 or fewer comma-separated "
			     << "arguments to --n-ceil option, got "
			     << args.size() << endl;
			throw 1;
		}
		polstr += ("MIN=" + args[0]);
		if(args.size() > 1) {
			polstr += ("," + args[1]);
		}
		if(args.size() > 2) {
			polstr += ("," + args[2]);
		}
		break;
	}
	case ARG_DESC: printArgDesc(cout); throw 0;
	case 'S': outfile = arg; break;
	case 'U': {
		EList<string> args;
		tokenize(arg, ",", args);
		for(size_t i = 0; i < args.size(); i++) {
			queries.push_back(args[i]);
		}
		break;
	}
#ifdef USE_SRA
        case ARG_SRA_ACC: {
		tokenize(arg, ",", sra_accs);
		set_format(format, SRA_FASTA);
		break;
        }
#endif
        case ARG_VERSION: showVersion = 1; break;
	default:
		printUsage(cerr);
		throw 1;
	}
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	saw_M = false;
	saw_a = false;
	saw_k = false;
	saw_trim3 = false;
	saw_trim5 = false;
	saw_trim_to = false;
	saw_bam = false;
	saw_preserve_tags = false;
	saw_align_paired_reads = false;
	presetList.clear();
	if(startVerbose) { cerr << "Parsing options: "; logTime(cerr, true); }
	while(true) {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		const char * arg = optarg;
		if(next_option == EOF) {
			if(extra_opts_cur < extra_opts.size()) {
				next_option = extra_opts[extra_opts_cur].first;
				arg = extra_opts[extra_opts_cur].second.c_str();
				extra_opts_cur++;
			} else {
				break;
			}
		}
		parseOption(next_option, arg);
	}

	if (!localAlign && scUnMapped) {
		cerr << "ERROR: --soft-clipped-unmapped-tlen can only be set for local alignments." << endl;
		exit(1);
	}

	if ((saw_trim3 || saw_trim5) && saw_trim_to) {
		cerr << "ERROR: --trim5/--trim3 and --trim-to are mutually exclusive "
		     << "options." << endl;
		exit(1);
	}

	if (!saw_bam && saw_preserve_tags) {
		cerr << "--preserve_tags can only be used when aligning BAM reads." << endl;
		exit(1);
	}

	if (!saw_bam && saw_align_paired_reads) {
		cerr << "--align-paired-reads can only be used when aligning BAM reads." << endl;
		exit(1);
	}
	// Now parse all the presets.  Might want to pick which presets version to
	// use according to other parameters.
	unique_ptr<Presets> presets(new PresetsV0());
	// Apply default preset
	if (presetList.empty())
		polstr = applyPreset(defaultPreset, *presets.get()) + polstr;
	else {
		for (size_t i = presetList.size(); i != 0; i--)
			polstr = applyPreset(presetList[i-1], *presets.get()) + polstr;
	}
	for(size_t i = 0; i < extra_opts.size(); i++) {
		next_option = extra_opts[extra_opts_cur].first;
		const char *arg = extra_opts[extra_opts_cur].second.c_str();
		parseOption(next_option, arg);
	}
	// Remove initial semicolons
	while(!polstr.empty() && polstr[0] == ';') {
		polstr = polstr.substr(1);
	}
	if(gVerbose) {
		cerr << "Final policy string: '" << polstr.c_str() << "'" << endl;
	}
	size_t failStreakTmp = 0;
	SeedAlignmentPolicy::parseString(
		polstr,
		localAlign,
		noisyHpolymer,
		ignoreQuals,
		bonusMatchType,
		bonusMatch,
		penMmcType,
		penMmcMax,
		penMmcMin,
		penNType,
		penN,
		penRdGapConst,
		penRfGapConst,
		penRdGapLinear,
		penRfGapLinear,
		scoreMin,
		nCeil,
		penNCatPair,
		multiseedMms,
		multiseedLen,
		msIval,
		failStreakTmp,
		nSeedRounds);
	if(failStreakTmp > 0) {
		maxEeStreak = failStreakTmp;
		maxUgStreak = failStreakTmp;
		maxDpStreak = failStreakTmp;
	}
	if(saw_a || saw_k) {
		msample = false;
		mhits = 0;
	} else {
		assert_gt(mhits, 0);
		msample = true;
	}
	if (format == UNKNOWN)
		set_format(format, FASTQ);
	if(mates1.size() != mates2.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
		     << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
		     << "sequences must be specified with -1 and -2." << endl;
		throw 1;
	}
	if(interleaved && (format != FASTA && format != FASTQ)) {
		cerr << "Error: --interleaved only works in combination with FASTA (-f) and FASTQ (-q) formats." << endl;
		throw 1;
	}
        if (samAppendComment && (format != FASTA && format != FASTQ)) {
		cerr << "Error --sam-append-comment only works with FASTA (-f) and FASTQ (-q) formats. " << endl;
		throw 1;
        }
        if(qualities.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with -Q but -f was not" << endl
		     << "enabled.  -Q works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q1 but -f was not" << endl
		     << "enabled.  --Q1 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities2.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q2 but -f was not" << endl
		     << "enabled.  --Q2 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() > 0 && mates1.size() != qualities1.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << qualities1.size() << endl
		     << "quality files were specified with --Q1.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -1 and --Q1." << endl;
		throw 1;
	}
	if(qualities2.size() > 0 && mates2.size() != qualities2.size()) {
		cerr << "Error: " << mates2.size() << " mate files/sequences were specified with -2, but " << qualities2.size() << endl
		     << "quality files were specified with --Q2.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -2 and --Q2." << endl;
		throw 1;
	}
	if(!rgs.empty() && rgid.empty()) {
		cerr << "Warning: --rg was specified without --rg-id also "
		     << "being specified.  @RG line is not printed unless --rg-id "
		     << "is specified." << endl;
	}
	// Check for duplicate mate input files
	if(format != CMDLINE) {
		for(size_t i = 0; i < mates1.size(); i++) {
			for(size_t j = 0; j < mates2.size(); j++) {
				if(mates1[i] == mates2[j] && !gQuiet) {
					cerr << "Warning: Same mate file \"" << mates1[i].c_str() << "\" appears as argument to both -1 and -2" << endl;
				}
			}
		}
	}
	if(useShmem && useMm && !gQuiet) {
		cerr << "Warning: --shmem overrides --mm..." << endl;
		useMm = false;
	}
	if(gGapBarrier < 1) {
		cerr << "Warning: --gbar was set less than 1 (=" << gGapBarrier
		     << "); setting to 1 instead" << endl;
		gGapBarrier = 1;
	}
	if(bonusMatch > 0 && !scoreMin.alwaysPositive()) {
		cerr << "Error: the match penalty is greater than 0 (" << bonusMatch
		     << ") but the --score-min function can be less than or equal to "
		     << "zero.  Either let the match penalty be 0 or make --score-min "
		     << "always positive." << endl;
		throw 1;
	}
	if(multiseedMms >= multiseedLen) {
		assert_gt(multiseedLen, 0);
		cerr << "Warning: seed mismatches (" << multiseedMms
		     << ") is less than seed length (" << multiseedLen
		     << "); setting mismatches to " << (multiseedMms-1)
		     << " instead" << endl;
		multiseedMms = multiseedLen-1;
	}
	sam_print_zm = sam_print_zm && bowtie2p5;
#ifndef NDEBUG
	if(!gQuiet) {
		cerr << "Warning: Running in debug mode.  Please use debug mode only "
		     << "for diagnosing errors, and not for typical use of Bowtie 2."
		     << endl;
	}
#endif
}

static const char *argv0 = NULL;

#if 0
/// Create a PatternSourcePerThread for the current thread according
/// to the global params and return a pointer to it
static PatternSourcePerThreadFactory*
createPatsrcFactory(
	PatternComposer& patcomp,
	const PatternParams& pp,
	int tid)
{
	PatternSourcePerThreadFactory *patsrcFact;
	patsrcFact = new PatternSourcePerThreadFactory(patcomp, pp, tid);
	assert(patsrcFact != NULL);
	return patsrcFact;
}
#endif

#define PTHREAD_ATTRS (PTHREAD_CREATE_JOINABLE | PTHREAD_CREATE_DETACHED)

static PatternSourceReadAheadFactory* multiseed_readahead_factory;
static Ebwt*                    multiseed_ebwtFw;
static Ebwt*                    multiseed_ebwtBw;
static Scoring*                 multiseed_sc;
static BitPairReference*        multiseed_refs;
static AlnSink*                 multiseed_msink;

// Cyclic rotations
#define ROTL(n, x) (((x) << (n)) | ((x) >> (32-n)))
#define ROTR(n, x) (((x) >> (n)) | ((x) << (32-n)))

static inline void printMmsSkipMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1,
	int seedmms)
{
	ostringstream os;
	if(paired) {
		os << "Warning: skipping mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because length (" << (mate1 ? ps.read_a().patFw.length() : ps.read_b().patFw.length())
		   << ") <= # seed mismatches (" << seedmms << ")" << endl;
	} else {
		os << "Warning: skipping read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because length (" << (mate1 ? ps.read_a().patFw.length() : ps.read_b().patFw.length())
		   << ") <= # seed mismatches (" << seedmms << ")" << endl;
	}
	cerr << os.str().c_str();
}

static inline void printLenSkipMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if(paired) {
		os << "Warning: skipping mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because it was < 2 characters long" << endl;
	} else {
		os << "Warning: skipping read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "' because it was < 2 characters long" << endl;
	}
	cerr << os.str().c_str();
}

static inline void printLocalScoreMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if(paired) {
		os << "Warning: minimum score function gave negative number in "
		   << "--local mode for mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	} else {
		os << "Warning: minimum score function gave negative number in "
		   << "--local mode for read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	}
	cerr << os.str().c_str();
}

static inline void printEEScoreMsg(
	const PatternSourcePerThread& ps,
	bool paired,
	bool mate1)
{
	ostringstream os;
	if(paired) {
		os << "Warning: minimum score function gave positive number in "
		   << "--end-to-end mode for mate #" << (mate1 ? '1' : '2')
		   << " of read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	} else {
		os << "Warning: minimum score function gave positive number in "
		   << "--end-to-end mode for read '" << (mate1 ? ps.read_a().name : ps.read_b().name)
		   << "; setting to 0 instead" << endl;
	}
	cerr << os.str().c_str();
}

/**
 * Initialize the minsc and maxpen arrays given information about the reads,
 * the alignment policy and the scoring scheme.
 */
static void setupMinScores(
	const PatternSourcePerThread& ps,
	bool paired,
	const Scoring& sc,
	const size_t *rdlens,
	TAlScore *minsc,
	TAlScore *maxpen)
{
	{
		minsc[0] = scoreMin.f<TAlScore>(rdlens[0]);
		if(paired) minsc[1] = scoreMin.f<TAlScore>(rdlens[1]);
		{
			if(minsc[0] > 0) {
				if(!gQuiet) printEEScoreMsg(ps, paired, true);
				minsc[0] = 0;
			}
			if(paired && minsc[1] > 0) {
				if(!gQuiet) printEEScoreMsg(ps, paired, false);
				minsc[1] = 0;
			}
		}
	}
	// Given minsc, calculate maxpen
	{
		assert_leq(minsc[0], 0);
		maxpen[0] = -minsc[0];
		if(paired) {
			assert_leq(minsc[1], 0);
			maxpen[1] = -minsc[1];
		} else {
			maxpen[1] = std::numeric_limits<TAlScore>::min();
		}
	}
}

#ifdef PER_THREAD_TIMING
/// Based on http://stackoverflow.com/questions/16862620/numa-get-current-node-core
void get_cpu_and_node(int& cpu, int& node) {
	unsigned long a,d,c;
	__asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
	node = (c & 0xFFF000)>>12;
	cpu = c & 0xFFF;
}
#endif

static inline bool have_next_read(std::unique_ptr<PatternSourceReadAhead> &g_psrah) {
  bool have_read;
  if (g_psrah) { // default case, we have already a ps
	  PatternSourcePerThread* ps = g_psrah.get()->ptr();
	  have_read = ps->nextReadPairReady();
	  if (!have_read) {
	     // finished current batch, but there may be more
	     g_psrah.reset(); // destroy the old object before allocating the new one
	     g_psrah.reset(new PatternSourceReadAhead(*multiseed_readahead_factory));
	     ps = g_psrah.get()->ptr();
	     have_read = ps->nextReadPairReady();
  	  }
  } else { // first read, just create ps
	     g_psrah.reset(new PatternSourceReadAhead(*multiseed_readahead_factory));
	     PatternSourcePerThread* ps = g_psrah.get()->ptr();
	     have_read = ps->nextReadPairReady();
  }
  return have_read;
}

// a couple helper functions to allow default constructors to be used
class AlignmentCacheIfaceBT2 : public AlignmentCacheIface {
public:
	AlignmentCacheIfaceBT2()
	: AlignmentCacheIface(new AlignmentCache(seedCacheCurrentMB * 1024 * 1024))
	{}

	virtual ~AlignmentCacheIfaceBT2() {delete current_;}

	AlignmentCacheIfaceBT2(const AlignmentCacheIfaceBT2& other) = delete;
	AlignmentCacheIfaceBT2& operator=(const AlignmentCacheIfaceBT2& other) = delete;
};

class SwDriverBT2 : public SwDriver {
public:
	SwDriverBT2()
	: SwDriver(exactCacheCurrentMB * 1024 * 1024) {}

	SwDriverBT2(SwDriverBT2&& other) = default;
	SwDriverBT2& operator=(SwDriverBT2&& other) noexcept = default;

	SwDriverBT2(const SwDriverBT2& other) = delete;
	SwDriverBT2& operator=(const SwDriverBT2& other) = delete;
};

class AlnSinkWrapOne : public AlnSinkWrap {
public:

	AlnSinkWrapOne(
		AlnSink& g,                // AlnSink being wrapped
		const ReportingParams& rp, // Parameters governing reporting
		Mapq& mapq,                // Mapq calculator
		size_t threadId,           // Thread ID
		BTPerThreadAllocators& allocs) :
	AlnSinkWrap(g, rp, mapq, threadId, allocs) ,
	prm() ,
	rpm() {}

	// Set filter subsets, reset counters and return overal filter
	bool setAndComputeFilter(const bool nfilt_, const bool scfilt_, const bool lenfilt_, const bool qcfilt_) {
		nfilt   = nfilt_;
		scfilt  = scfilt_;
		lenfilt = lenfilt_;
		qcfilt  = qcfilt_;

		// also reset counters
		seedsTried = 0;
		nUniqueSeeds = 0;
		nRepeatSeeds = 0;
		seedHitTot = 0;
		for (size_t i=0; i<2; i++) {
			seedsTriedMS[i] = 0;
			nUniqueSeedsMS[i] = 0;
			nRepeatSeedsMS[i] = 0;
			seedHitTotMS[i] = 0;
		}

		return nfilt_ && scfilt_ && lenfilt_ && qcfilt_;
	}

	void updateSHSCounters(const SeedResults& shs) {
		// shs contain what we need to know to update our seed
		// summaries for this seeding
		if(!shs.empty()) {
			nUniqueSeeds += shs.numUniqueSeeds();
			nUniqueSeedsMS[0] += shs.numUniqueSeedsStrand(true);
			nUniqueSeedsMS[1] += shs.numUniqueSeedsStrand(false);
			nRepeatSeeds += shs.numRepeatSeeds();
			nRepeatSeedsMS[0] += shs.numRepeatSeedsStrand(true);
			nRepeatSeedsMS[1] += shs.numRepeatSeedsStrand(false);
			seedHitTot += shs.numElts();
			seedHitTotMS[0] += shs.numEltsFw();
			seedHitTotMS[1] += shs.numEltsRc();
		}
	}

	/**
	 * Inform global, shared AlnSink object that we're finished with
	 * this read.  The global AlnSink is responsible for updating
	 * counters, creating the output record, and delivering the record
	 * to the appropriate output stream.
	 */
	void finishReadOne(
		const SeedResults *sr,          // seed alignment results
		bool               exhaust,     // mate  exhausted?
		RandomSource&      rnd,         // pseudo-random generator
		const Scoring& sc,              // scoring scheme
		bool suppressSeedSummary = true,
		bool suppressAlignments = false,
		bool scUnMapped = false,
		bool xeq = false) {
			finishRead(
					sr,                   // seed results for mate 1
					NULL,                 // seed results for mate 2 (NULL, since unpired)
					exhaust,              // exhausted seed hits for mate 1?
					false,                // exhausted seed hits for mate 2? (false, since upaired)
					nfilt,
					false,
					scfilt,
					false,
					lenfilt,
					true,
					qcfilt,
					true,
					rnd,                  // pseudo-random generator
					rpm,                  // reporting metrics
					prm,                  // per-read metrics
					sc,                   // scoring scheme
					true,                 // suppress seed summaries?
					false,                // suppress alignments?
					scUnMapped,           // Consider soft-clipped bases unmapped when calculating TLEN
					xeq);
	}

	void finalizePRM(size_t totnucs) {
		if(seedsTried > 0) {
			float invSeedsTried = (float) 1.0f/seedsTried;
			prm.seedPctUnique = nUniqueSeeds * invSeedsTried;
			prm.seedPctRep    = nRepeatSeeds * invSeedsTried;
			prm.seedHitAvg    = seedHitTot   * invSeedsTried;
		} else {
			prm.seedPctUnique = -1.0f;
			prm.seedPctRep = -1.0f;
			prm.seedHitAvg = -1.0f;
		}
		for(int i = 0; i < 2; i++) {
			if(seedsTriedMS[i] > 0) {
				float invSeedsTried = (float) 1.0f/seedsTriedMS[i];
				prm.seedPctUniqueMS[i] = nUniqueSeedsMS[i] * invSeedsTried;
				prm.seedPctRepMS[i]    = nRepeatSeedsMS[i] * invSeedsTried;
				prm.seedHitAvgMS[i]    = seedHitTotMS[i]   * invSeedsTried;
			} else {
				prm.seedPctUniqueMS[i] = -1.0f;
				prm.seedPctRepMS[i] = -1.0f;
				prm.seedHitAvgMS[i] = -1.0f;
			}
		}
		for(int i = 2; i < 4; i++) {
				prm.seedPctUniqueMS[i] = -1.0f;
				prm.seedPctRepMS[i] = -1.0f;
				prm.seedHitAvgMS[i] = -1.0f;
		}
		if (totnucs > 0) {
			float invtotnucs = (float) 1.0f/totnucs;
			prm.seedsPerNuc = seedsTried * invtotnucs;
			for(int i = 0; i < 2; i++) {
				prm.seedsPerNucMS[i] = seedsTriedMS[i] * invtotnucs;
			}
			for(int i = 2; i < 4; i++) {
				prm.seedsPerNucMS[i] = 0;
			}
		} else {
			prm.seedsPerNuc = -1;
			for(int i = 0; i < 4; i++) {
				prm.seedsPerNucMS[i] = -1;
			}
		}
	}

public:
	// per-read metrics
	PerReadMetrics prm;
	// global metrics, per-thread but will merge at the end
	ReportingMetrics rpm;

	// Keep track of whether mate was filtered out due Ns last time
	bool nfilt;
	// Keep track of whether mate was filtered out due to not having
	// enough characters to rise about the score threshold.
	bool scfilt;
	// Keep track of whether mate was filtered out due to not having
	// more characters than the number of mismatches permitted in a seed.
	bool lenfilt;
	// Keep track of whether mate was were filtered out by upstream qc
	bool qcfilt;

	// counters
	size_t seedsTried;
	size_t seedsTriedMS[2];
	size_t nUniqueSeeds;
	size_t nRepeatSeeds;
	size_t seedHitTot;
	size_t nUniqueSeedsMS[2];
	size_t nRepeatSeedsMS[2];
	size_t seedHitTotMS[2];

};


class msWorkerObjs {
public:
	msWorkerObjs() = default;

	msWorkerObjs(msWorkerObjs&& o) = default;
	msWorkerObjs(const msWorkerObjs& o) = delete;

	msWorkerObjs &operator=(msWorkerObjs&& o) noexcept = default;
	msWorkerObjs &operator=(const msWorkerObjs& o) noexcept = delete;

	void set_alloc(BTAllocator *alloc, bool propagate_alloc=true)
	{
		al.set_alloc(alloc,propagate_alloc);
		sd.set_alloc(alloc,propagate_alloc);
		sw.set_alloc(alloc,propagate_alloc);
		shs.set_alloc(alloc,propagate_alloc);
		seed.set_alloc(alloc,propagate_alloc);
	}

	// Interfaces for alignment and seed caches
	AlignmentCacheIfaceBT2 ca;

	SeedAligner al;
	SwDriverBT2 sd;
	SwAligner sw;
	SeedResults shs;
	RandomSource rnd;
	EList<Seed> seed;
};

class msWorkerConsts {
public:

	msWorkerConsts()
		: ref(*multiseed_refs)
		, sc(*multiseed_sc)
		, ebwtFw(*multiseed_ebwtFw)
		, ebwtBw(multiseed_ebwtBw)
		, nofw( gNofw )
		, norc( gNorc )
		, extend( doExtend )
		, doEnable8( enable8 )
		, tri( doTri )
		, doTighten( tighten)
		, maxHalf( maxhalf )
		, cMinLen( cminlen )
		, cPow2( cpow2 )
		, mxKHMul( (khits > 1) ? (khits-1) : 0 )
		, streak( maxDpStreak + mxKHMul * maxStreakIncr )
		, mxDp(   maxDp       + mxKHMul * maxItersIncr  )
		, mxUg(   maxUg       + mxKHMul * maxItersIncr  )
		, mxIter( maxIters    + mxKHMul * maxItersIncr  )
	{}

	const BitPairReference& ref;
	const Scoring&          sc;
	const Ebwt&             ebwtFw;
	const Ebwt* const       ebwtBw;

	const bool nofw;
	const bool norc;

	const bool extend;
	const bool doEnable8;
	const bool tri; 
	const int  doTighten; 
	const size_t maxHalf;
	const size_t cMinLen; 
	const size_t cPow2;

	// Calculate streak length
	const size_t mxKHMul;
	const size_t streak;
	const size_t mxDp;
	const size_t mxUg;
	const size_t mxIter;
};

/**
 * Called once per thread.  Sets up per-thread pointers to the shared global
 * data structures, creates per-thread structures, then enters the alignment
 * loop.  The general flow of the alignment loop is:
 *
 * - Get the next read/pair
 * - Check if this read/pair is identical to the previous
 *   + If identical, check whether we can skip any or all alignment stages.  If
 *     we can skip all stages, report the result immediately and move to next
 *     read/pair
 *   + If not identical, continue
 * -
 */
 /* Work in progress: There is only one such invocation, parallelization inside */
static void multiseedSearchWorker(const size_t num_parallel_tasks) {
	assert(multiseed_ebwtFw != NULL);
	assert(multiseedMms == 0 || multiseed_ebwtBw != NULL);
	AlnSink&                msink    = *multiseed_msink;

	constexpr bool paired = false;

	BTAllocator worker_alloc;
	BTPerThreadAllocators mate_allocs(num_parallel_tasks);

	// allocate to heap to make it GPU accessible
	const msWorkerConsts *msconsts = new msWorkerConsts;
 
	{
		// Sinks: these are so that we can print tables encoding counts for
		// events of interest on a per-read, per-seed, per-join, or per-SW
		// level.  These in turn can be used to diagnose performance
		// problems, or generally characterize performance.

		// Instantiate an object for holding reporting-related parameters.
		// Use new to make it GPU-accessible
		ReportingParams* rp= new ReportingParams(
			khits,             // -k
			mhits,             // -m/-M
			0,                 // penalty gap (not used now)
			msample,           // true -> -M was specified, otherwise assume -m
			gReportDiscordant, // report discordang paired-end alignments?
			gReportMixed);     // report unpaired alignments for paired reads?

		// Note: Cannot use std:vector due to GPU compute not having access to the CPU stack

		// Instantiate a mapping quality calculator
		std::vector< unique_ptr<Mapq> > bmapq;
		std::vector<AlnSinkWrapOne> v_msinkwrap;

		bmapq.reserve(num_parallel_tasks);
		v_msinkwrap.reserve(num_parallel_tasks);
		
		for (size_t mate=0; mate<num_parallel_tasks; mate++) {
			bmapq.emplace_back(new_mapq(mapqv, scoreMin, msconsts->sc));
			v_msinkwrap.emplace_back(
							msink,        // global sink
							*rp,           // reporting parameters
							*bmapq[mate], // MAPQ calculator
							mate,         // thread id
							mate_allocs); // memory pools
		}
		AlnSinkWrapOne *g_msinkwrap = v_msinkwrap.data();

		// Major per-thread objects
		msWorkerObjs* g_msobjs = new msWorkerObjs[num_parallel_tasks];
		for (size_t mate=0; mate<num_parallel_tasks; mate++) g_msobjs[mate].set_alloc(&(mate_allocs[mate]));

		// These arrays move between CPU and GPU often
		// Keep them together

		// Whether we're done with g_psrah mate
		bool *done_reading = new(worker_alloc.allocate(num_parallel_tasks,sizeof(bool))) bool[num_parallel_tasks];
		for (size_t mate=0; mate<num_parallel_tasks; mate++) {
			done_reading[mate] = false;
		}
		// Whether we're done with mate
		bool *done = new(worker_alloc.allocate(num_parallel_tasks,sizeof(bool))) bool[num_parallel_tasks];

		size_t* nelt = new(worker_alloc.allocate(num_parallel_tasks,sizeof(size_t))) size_t[num_parallel_tasks];

		// Used by thread with threadid == 1 to measure time elapsed
		time_t iTime = time(0);

		// Keep track of whether last search was exhaustive for mates 1 and 2
		bool *exhaustive = new bool[num_parallel_tasks];
		// Keep track of whether mates 1/2 were filtered out last time through
		bool *filt = new bool[num_parallel_tasks];


		// read object
		typedef Read* ReadPtr;
		Read* *rds = new ReadPtr[num_parallel_tasks];
		// Calculate the minimum valid score threshold for the read
		TAlScore* minsc = new TAlScore[num_parallel_tasks];

		int* nceil = new int[num_parallel_tasks];

                std::vector< std::unique_ptr<PatternSourceReadAhead> > g_psrah(num_parallel_tasks);
		std::vector<PatternSourcePerThread*> ps(num_parallel_tasks);

		// matemap vector is local to CPU
		std::vector<size_t> matemap(num_parallel_tasks);
		size_t current_num_parallel_tasks = num_parallel_tasks;

		while (true) { // will exit with a break
		   mate_allocs.ensure_spare();
		   {
			bool found_unread = false;
			// Note: Will use mate to distinguish between tread-specific elements
			// Note2: we could potentially run this as OMP, but the overhead is likely higher than the speedup
			for (size_t mate=0; mate<num_parallel_tasks; mate++) {
			    msWorkerObjs& msobj = g_msobjs[mate];
			    while (!done_reading[mate]) { // External loop, including filtering

				while (!done_reading[mate]) { // Internal loop, just buffer reads and retries
					done_reading[mate] = !have_next_read(g_psrah[mate]);
					if (!done_reading[mate]) {
						pair<bool, bool> ret = g_psrah[mate].get()->nextReadPair();
						bool success = ret.first;
						if(!success) {
							done_reading[mate] = ret.second;
							if (!done_reading[mate]) {
								// just retry the inner while
								continue;
							}
						}
					}
					// if we got here, do not retry inner loop anymore
					break;
				} // internal while loop

				if (!done_reading[mate]) {
					ps[mate] = g_psrah[mate].get()->ptr();
					rds[mate] = &ps[mate]->read_a();

					AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
					TReadId rdid = rds[mate]->rdid;

					// Align this read/pair
					//
					// Check if there is metrics reporting for us to do.
					//
					msinkwrap.prm.reset(); // per-read metrics
					msinkwrap.prm.doFmString = false;

					msobj.ca.nextRead(); // clear the cache
					assert(!msobj.ca.aligning());
					const size_t rdlen = rds[mate]->length();
					msinkwrap.nextRead(
						rds[mate],
						NULL,
						rdid,
						msconsts->sc.qualitiesMatter());
					assert(msinkwrap.inited());
					minsc[mate] = scoreMin.f<TAlScore>(rdlen);
					if(minsc[mate] > 0) {
						if(!gQuiet) printEEScoreMsg(*ps[mate], paired, true);
						minsc[mate] = 0;
					}
					{
						// N filter; does the read have too many Ns?
						bool nfilt = msconsts->sc.nFilter(rds[mate]->patFw);

						// Score filter; does the read enough character to rise above
						// the score threshold?
						bool scfilt =msconsts->sc.scoreFilter(minsc[mate], rdlen);
						bool lenfilt = true;
						if(rdlen <= (size_t)multiseedMms || rdlen < 2) {
							if(!gQuiet) printMmsSkipMsg(*ps[mate], paired, true, multiseedMms);
							lenfilt = false;
						}
						if(rdlen < 2) {
							if(!gQuiet) printLenSkipMsg(*ps[mate], paired, true);
							lenfilt = false;
						}
						bool qcfilt = true;
						if(qcFilter) {
							qcfilt = (rds[mate]->filter != '0');
						}
						filt[mate] = msinkwrap.setAndComputeFilter(nfilt, scfilt, lenfilt, qcfilt);
					}
					msinkwrap.prm.nFilt += (filt[mate] ? 0 : 1);
					// For each mate...
					assert(msinkwrap.empty());
					msobj.sd.nextRead(false, rdlen, 0); // SwDriver
					exhaustive[mate] = false;
					msobj.rnd.init(rds[mate]->seed);

					msinkwrap.prm.maxDPFails = msconsts->streak;
					assert_gt(msconsts->streak, 0);
					msobj.shs.clear();

					// Whether we're done with mate
					done[mate] = !filt[mate];
					// Increment counters according to what got filtered
					if(filt[mate]) { // done[mate] == false
						msobj.shs.nextRead(*rds[mate]);
						assert(msobj.shs.empty());

						nceil[mate] = min((int)nCeil.f<int>((double)rdlen), (int)rdlen);
						nelt[mate] = 0;
						found_unread = true;
						break; // we found a good read, can get out
					} else { // done[mate] == true
						// we are already done... finalize and read next

						msinkwrap.finalizePRM(0);

						msinkwrap.finishReadOne(
							&msobj.shs,           // seed results for mate 1
							exhaustive[mate],     // exhausted seed hits for mate 1?
							msobj.rnd,            // pseudo-random generator
							msconsts->sc,                   // scoring scheme
							true,                 // suppress seed summaries?
							false,                // suppress alignments?
							scUnMapped,           // Consider soft-clipped bases unmapped when calculating TLEN
							xeq);
					}

		    		} // if (!done_reading[mate])

			   } // external while loop

			} // for mate - found_unread
			if (!found_unread) break; // nothing else to do
		   }

		   // always call ensure_spare from main CPU thread
		   mate_allocs.ensure_spare();

		   // ASSERT: done[:] == false
		   matemap.clear();
		   for (size_t mate=0; mate<num_parallel_tasks; mate++) {
			if ( (!done_reading[mate]) && (!done[mate]) )  {
				matemap.push_back(mate);
			}
		   }
		   current_num_parallel_tasks = matemap.size();

		   // Expected fraction of mates: 50%
		   // we can do all of the "mates" in parallel
#if defined(USE_ACC_STDPAR)
		   std::for_each(std::execution::par_unseq, matemap.begin(), matemap.end(),
			[g_msobjs,msconsts,nelt,rds,minsc](size_t mate) mutable {
#else

#if defined(USE_ACC_CPU)
#pragma omp parallel for default(shared)
#else
#pragma acc parallel loop gang vector vector_length(32)
#endif
		   for (size_t fidx=0; fidx<current_num_parallel_tasks; fidx++) {
			const size_t mate=matemap[fidx];
#endif
					msWorkerObjs& msobj = g_msobjs[mate];
					bool not_yfw = false;
					bool not_yrc = false;
					// 1-mismatch
					// Clear out the exact hits
					msobj.al.oneMmSearch(
									&msconsts->ebwtFw,        // BWT index
									msconsts->ebwtBw,         // BWT' index
									*rds[mate],     // read
									msconsts->sc,             // scoring scheme
									minsc[mate],    // minimum score
									not_yfw,        // don't align forward read
									not_yrc,        // don't align revcomp read
									false,          // do exact match
									true,           // do 1mm
									msobj.shs);     // seed hits (hits installed here)
					nelt[mate] = msobj.shs.num1mmE2eHits();
		   } // for mate
#if defined(USE_ACC_STDPAR)
		   );
#endif
	
		   // always call ensure_spare from main CPU thread
		   mate_allocs.ensure_spare();

		   // nelt[:] has been modified by the previous loop, but we only care about the non-0 ones
		   matemap.clear();
		   for (size_t mate=0; mate<num_parallel_tasks; mate++) {
			if ((!done_reading[mate])  && (!done[mate]) && (nelt[mate] != 0)) {
				matemap.push_back(mate);
			}
		   }
		   current_num_parallel_tasks = matemap.size();

		   // Expected fraction of mates: 20%
		   // we can do all of the "mates" in parallel
#pragma omp parallel for default(shared)
		   for (size_t fidx=0; fidx<current_num_parallel_tasks; fidx++) {
			const size_t mate=matemap[fidx];
					msWorkerObjs& msobj = g_msobjs[mate];
					AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
					const size_t rdlen = rds[mate]->length();

					// 1-mismatch, part 2
							// Unpaired dynamic programming driver
							int	ret = msobj.sd.extendSeeds(
									*rds[mate],     // read
									true,           // mate #1?
									msobj.shs,      // seed hits
									msconsts->ebwtFw,         // bowtie index
									msconsts->ebwtBw,         // rev bowtie index
									msconsts->ref,            // packed reference strings
									msobj.sw,       // dynamic prog aligner
									msconsts->sc,             // scoring scheme
									-1,             // # mms allowed in a seed
									0,              // length of a seed
									0,              // interval between seeds
									minsc[mate],    // minimum score for valid
									nceil[mate],    // N ceil for anchor
									msconsts->maxHalf,        // max width on one DP side
									msconsts->mxIter,         // max extend loop iters
									msconsts->mxUg,           // max # ungapped extends
									msconsts->mxDp,           // max # DPs
									msconsts->streak,         // stop after streak of this many end-to-end fails
									msconsts->streak,         // stop after streak of this many ungap fails
									msconsts->extend,       // extend seed hits
									msconsts->doEnable8,        // use 8-bit SSE where possible
									msconsts->cMinLen,        // checkpoint if read is longer
									msconsts->cPow2,          // checkpointer interval, log2
									msconsts->tri,          // triangular mini-fills?
									msconsts->doTighten,        // -M score tightening mode
									msobj.ca,       // seed alignment cache
									msobj.rnd,      // pseudo-random source
									msinkwrap.prm,  // per-read metrics
									&msinkwrap,     // for organizing hits
									true,           // report hits once found
									exhaustive[mate]);

							assert_gt(ret, 0);
							// Clear out the 1mm hits so that we don't try to
							// extend them again later!
							msobj.shs.clear1mmE2eHits();
							if(ret == EXTEND_EXHAUSTED_CANDIDATES) {
								// Not done yet
							} else if(ret == EXTEND_POLICY_FULFILLED) {
								// Policy is satisfied for this mate at least
								if(msinkwrap.state().doneWithMate(true)) {
									done[mate] = true;
								}
							} else if(ret == EXTEND_PERFECT_SCORE) {
								// We exhausted this mode at least
								done[mate] = true;
							} else if(ret == EXTEND_EXCEEDED_HARD_LIMIT) {
								// We exceeded a per-read limit
								done[mate] = true;
							} else if(ret == EXTEND_EXCEEDED_SOFT_LIMIT) {
								// Not done yet
							} else {
								//
								cerr << "Bad return value: " << ret << endl;
								// We do not have clean exception handling, so just report to the user 
								done[mate] = true;
							}
							if(!done[mate]) {
								TAlScore perfectScore = msconsts->sc.perfectScore(rdlen);
								if(minsc[mate] == perfectScore) {
									done[mate] = true;
								}
							}
		   } // for mate

		   // Expected fraction of mates: 95%+
                   for(size_t roundi = 0; roundi < nSeedRounds; roundi++) {	
			// always call ensure_spare from main CPU thread
		 	mate_allocs.ensure_spare();

		 	matemap.clear();
		  	for (size_t mate=0; mate<num_parallel_tasks; mate++) {
				const size_t rdlen = rds[mate]->length();
				// Calculate interval length for both mates
				const int interval = max((int)msIval.f<int>((double)rdlen), 1);
				// Calculate # seed rounds for each mate
				const size_t nrounds = min<size_t>(nSeedRounds, interval);
				if ((!done_reading[mate])  && (!done[mate]) && (roundi<nrounds) ) {
					matemap.push_back(mate);
				}
		   	}
			current_num_parallel_tasks = matemap.size();

		   	// we can do all of the "mates" in parallel
#pragma omp parallel for default(shared)
		   	for (size_t fidx=0; fidx<current_num_parallel_tasks; fidx++) {
				const size_t mate=matemap[fidx];
					msWorkerObjs& msobj = g_msobjs[mate];
					AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
					const size_t rdlen = rds[mate]->length();
					if(msinkwrap.state().doneWithMate(true)) {
						// Should never get in here, but just in case
						// Done with this mate
						done[mate] = true;
						continue;
					}

					// Calculate interval length for both mates
					const int interval = max((int)msIval.f<int>((double)rdlen), 1);
					// Calculate # seed rounds for each mate
					const size_t nrounds = min<size_t>(nSeedRounds, interval);
					Constraint gc = Constraint::penaltyFuncBased(scoreMin);
						msobj.ca.nextRead(); // Clear cache in preparation for new search
						msobj.shs.clearSeeds();
						assert(msobj.shs.empty());
						assert(msobj.shs.repOk(&msobj.ca.current()));
							size_t offset = (interval * roundi) / nrounds;
							assert(roundi == 0 || offset > 0);
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							//rnd.init(ROTL(rds[mate]->seed, 10));
							assert(msobj.shs.repOk(&msobj.ca.current()));
							// Set up seeds
							msobj.seed.clear();
							Seed::mmSeeds(
								multiseedMms,    // max # mms per seed
								multiseedLen,    // length of a multiseed seed
								msobj.seed,      // seeds
								gc);             // global constraint
							// Check whether the offset would drive the first seed
							// off the end
							if(offset > 0 && multiseedLen + offset > rds[mate]->length()) {
								done[mate] = true;
							}
			} // for mate
			// always call ensure_spare from main CPU thread
		 	mate_allocs.ensure_spare();

		 	matemap.clear();
		  	for (size_t mate=0; mate<num_parallel_tasks; mate++) {
				const size_t rdlen = rds[mate]->length();
				// Calculate interval length for both mates
				const int interval = max((int)msIval.f<int>((double)rdlen), 1);
				const size_t nrounds = min<size_t>(nSeedRounds, interval);
				if ((!done_reading[mate])  && (!done[mate]) && (roundi<nrounds) ) {
					matemap.push_back(mate);
				}
		   	}
			current_num_parallel_tasks = matemap.size();

		   	// we can do all of the "mates" in parallel
#pragma omp parallel for default(shared)
		   	for (size_t fidx=0; fidx<current_num_parallel_tasks; fidx++) {
				const size_t mate=matemap[fidx];
					msWorkerObjs& msobj = g_msobjs[mate];
					AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
					const size_t rdlen = rds[mate]->length();
					// Calculate interval length for both mates
					const int interval = max((int)msIval.f<int>((double)rdlen), 1);
					const size_t nrounds = min<size_t>(nSeedRounds, interval);
					size_t offset = (interval * roundi) / nrounds;
							std::pair<int, int> instFw, instRc;
							// Instantiate the seeds
							std::pair<int, int> inst = msobj.al.instantiateSeeds(
								msobj.seed,     // search seeds
								offset,         // offset to begin extracting
								interval,       // interval between seeds
								*rds[mate],     // read to align
								msconsts->sc,             // scoring scheme
								msconsts->nofw,          // don't align forward read
								msconsts->norc,          // don't align revcomp read
								msobj.ca,       // holds some seed hits from previous reads
								msobj.shs,      // holds all the seed hits
								instFw,
								instRc);
							assert(msobj.shs.repOk(&msobj.ca.current()));
							if(inst.first + inst.second == 0) {
								// No seed hits!  Done with this mate.
								assert(msobj.shs.empty());
								done[mate] = true;
							} else {
								msinkwrap.seedsTried += (inst.first + inst.second);
								msinkwrap.seedsTriedMS[0] = instFw.first + instFw.second;
								msinkwrap.seedsTriedMS[1] = instRc.first + instRc.second;
							}
			} // for mate
			// always call ensure_spare from main CPU thread
		 	mate_allocs.ensure_spare();

		 	matemap.clear();
		  	for (size_t mate=0; mate<num_parallel_tasks; mate++) {
				const size_t rdlen = rds[mate]->length();
				// Calculate interval length for both mates
				const int interval = max((int)msIval.f<int>((double)rdlen), 1);
				const size_t nrounds = min<size_t>(nSeedRounds, interval);
				if ((!done_reading[mate])  && (!done[mate]) && (roundi<nrounds) ) {
					matemap.push_back(mate);
				}
		   	}
			current_num_parallel_tasks = matemap.size();

		   	// we can do all of the "mates" in parallel
#pragma omp parallel for default(shared)
		   	for (size_t fidx=0; fidx<current_num_parallel_tasks; fidx++) {
				const size_t mate=matemap[fidx];
					msWorkerObjs& msobj = g_msobjs[mate];
					AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
							// Align seeds
							msobj.al.searchAllSeeds(
								msobj.seed,       // search seeds
								&msconsts->ebwtFw,          // BWT index
								msconsts->ebwtBw,           // BWT' index
								*rds[mate],       // read
								msconsts->sc,               // scoring scheme
								msobj.ca,         // alignment cache
								msobj.shs,        // store seed hits here
								msinkwrap.prm);   // per-read metrics
							assert(msobj.shs.repOk(&msobj.ca.current()));
							if(msobj.shs.empty()) {
								// No seed alignments!  Done with this mate.
								done[mate] = true;
							} else {
								// shs contain what we need to know to update our seed
								// summaries for this seeding
								msinkwrap.updateSHSCounters(msobj.shs);
							}

			} // for mate
			// always call ensure_spare from main CPU thread
		 	mate_allocs.ensure_spare();

		 	matemap.clear();
		  	for (size_t mate=0; mate<num_parallel_tasks; mate++) {
				const size_t rdlen = rds[mate]->length();
				// Calculate interval length for both mates
				const int interval = max((int)msIval.f<int>((double)rdlen), 1);
				const size_t nrounds = min<size_t>(nSeedRounds, interval);
				if ((!done_reading[mate])  && (!done[mate]) && (roundi<nrounds) ) {
					matemap.push_back(mate);
				}
		   	}
			current_num_parallel_tasks = matemap.size();

		   	// we can do all of the "mates" in parallel
#pragma omp parallel for default(shared)
		   	for (size_t fidx=0; fidx<current_num_parallel_tasks; fidx++) {
				const size_t mate=matemap[fidx];
					msWorkerObjs& msobj = g_msobjs[mate];
					AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
					if(msinkwrap.state().doneWithMate(true)) {
						// Should never get in here, but just in case
						// Done with this mate
						done[mate] = true;
						continue;
					}
					const size_t rdlen = rds[mate]->length();
					// Calculate interval length for both mates
					const int interval = max((int)msIval.f<int>((double)rdlen), 1);
								// Sort seed hits into ranks
								msobj.shs.rankSeedHits(msobj.rnd, msinkwrap.allHits());
								int ret = 0;
                                                                {
									// Unpaired dynamic programming driver
									ret = msobj.sd.extendSeeds(
										*rds[mate],     // read
										true,           // mate #1?
										msobj.shs,      // seed hits
										msconsts->ebwtFw,         // bowtie index
										msconsts->ebwtBw,         // rev bowtie index
										msconsts->ref,            // packed reference strings
										msobj.sw,       // dynamic prog aligner
										msconsts->sc,             // scoring scheme
										multiseedMms,   // # mms allowed in a seed
										multiseedLen,   // length of a seed
										interval,       // interval between seeds
										minsc[mate],    // minimum score for valid
										nceil[mate],    // N ceil for anchor
										msconsts->maxHalf,        // max width on one DP side
										msconsts->mxIter,         // max extend loop iters
										msconsts->mxUg,           // max # ungapped extends
										msconsts->mxDp,           // max # DPs
										msconsts->streak,         // stop after streak of this many end-to-end fails
										msconsts->streak,         // stop after streak of this many ungap fails
										msconsts->extend,       // extend seed hits
										msconsts->doEnable8,        // use 8-bit SSE where possible
										msconsts->cMinLen,        // checkpoint if read is longer
										msconsts->cPow2,          // checkpointer interval, log2
										msconsts->tri,          // triangular mini-fills?
										msconsts->doTighten,        // -M score tightening mode
										msobj.ca,       // seed alignment cache
										msobj.rnd,      // pseudo-random source
										msinkwrap.prm,  // per-read metrics
										&msinkwrap,     // for organizing hits
										true,           // report hits once found
										exhaustive[mate]);
								}
								assert_gt(ret, 0);
								if ((ret == EXTEND_EXHAUSTED_CANDIDATES) ||
								     (ret == EXTEND_EXCEEDED_SOFT_LIMIT) ||
                                                                     (ret == EXTEND_POLICY_FULFILLED)) {
									// Not done yet
								} else {
									done[mate] = true;
								}

			} // for mate

			// no point doing in parallel
		  	for (size_t mate=0; mate<num_parallel_tasks; mate++) {
				const size_t rdlen = rds[mate]->length();
				// Calculate interval length for both mates
				const int interval = max((int)msIval.f<int>((double)rdlen), 1);
				const size_t nrounds = min<size_t>(nSeedRounds, interval);
				if ((!done_reading[mate])  && (!done[mate]) && (roundi<nrounds) ) {
					msWorkerObjs& msobj = g_msobjs[mate];
					AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
							// We don't necessarily have to continue investigating both
							// mates.  We continue on a mate only if its average
							// interval length is high (> 1000)
							if(msobj.shs.averageHitsPerSeed() < seedBoostThresh) {
								done[mate] = true;
							}
							if(msinkwrap.state().doneWithMate(true)) {
								done[mate] = true;
							}
				}
			} // for mate
		   } // end loop over reseeding rounds

		   // always call ensure_spare from main CPU thread
		   mate_allocs.ensure_spare();

		   // we could do all of the "mates" in parallel, but the overhead is likely higher than the speedup
//pragma omp parallel for default(shared)
		   for (size_t mate=0; mate<num_parallel_tasks; mate++) {
			if (!done_reading[mate]) { // only do it for valid ones, to handle end tails
				msWorkerObjs& msobj = g_msobjs[mate];
				AlnSinkWrapOne& msinkwrap = g_msinkwrap[mate]; 
		
				const size_t rdlen = rds[mate]->length();
				size_t totnucs = 0;
				{
					size_t len = rdlen;
					if(!gNofw && !gNorc) {
						len *= 2;
					}
					totnucs += len;
				}
				msinkwrap.finalizePRM(totnucs);

				// Commit and report paired-end/unpaired alignments
				//uint32_t sd = rds[0]->seed ^ rds[1]->seed;
				//rnd.init(ROTL(sd, 20));
				msinkwrap.finishReadOne(
					&msobj.shs,           // seed results for mate 1
					exhaustive[mate],     // exhausted seed hits for mate 1?
					msobj.rnd,            // pseudo-random generator
					msconsts->sc,                   // scoring scheme
					true,                 // suppress seed summaries?
					false,                // suppress alignments?
					scUnMapped,           // Consider soft-clipped bases unmapped when calculating TLEN
					xeq);

			} // if (!done_reading[mate])
		   } // for mate
		} // while(true)

		// Merge in the metrics
		// Must be done sequentially
		for (size_t mate=0; mate<num_parallel_tasks; mate++) {
			msink.mergeMetricsUnsafe(g_msinkwrap[mate].rpm);
		}

		v_msinkwrap.clear();
		delete rp;

		delete[] exhaustive;
		delete[] filt;

		delete[] rds;
		delete[] minsc;
		delete[] nceil;

		delete[] g_msobjs;

		// do not delete objects managed by the allocator
		// accept memory leak
	}

	delete msconsts;

	return;
}

#ifdef SUPPORT_PAIRED
// NOTE: Unsupported, likely does not work
static void multiseedSearchWorkerPaired(const size_t num_parallel_tasks) {
	int tid = 0;
	assert(multiseed_ebwtFw != NULL);
	assert(multiseedMms == 0 || multiseed_ebwtBw != NULL);
	PatternSourceReadAheadFactory& readahead_factory =  *multiseed_readahead_factory;
	const Ebwt&             ebwtFw   = *multiseed_ebwtFw;
	const Ebwt*             ebwtBw   = multiseed_ebwtBw;
	const Scoring&          sc       = *multiseed_sc;
	const BitPairReference& ref      = *multiseed_refs;
	AlnSink&                msink    = *multiseed_msink;

	{
		// Sinks: these are so that we can print tables encoding counts for
		// events of interest on a per-read, per-seed, per-join, or per-SW
		// level.  These in turn can be used to diagnose performance
		// problems, or generally characterize performance.

		//const BitPairReference& refs   = *multiseed_refs;

		AlignmentCache scCurrent(seedCacheCurrentMB * 1024 * 1024);

		// Interfaces for alignment and seed caches
		AlignmentCacheIface ca(&scCurrent);

		// Instantiate an object for holding reporting-related parameters.
		ReportingParams rp(
			(allHits ? std::numeric_limits<THitInt>::max() : khits), // -k
			mhits,             // -m/-M
			0,                 // penalty gap (not used now)
			msample,           // true -> -M was specified, otherwise assume -m
			gReportDiscordant, // report discordang paired-end alignments?
			gReportMixed);     // report unpaired alignments for paired reads?

		// Instantiate a mapping quality calculator
		unique_ptr<Mapq> bmapq(new_mapq(mapqv, scoreMin, sc));

		// Make a per-thread wrapper for the global MHitSink object.
		AlnSinkWrap msinkwrap(
			msink,         // global sink
			rp,            // reporting parameters
			*bmapq,        // MAPQ calculator
			(size_t)tid);  // thread id

		SeedAligner al;
		SwDriver sd(exactCacheCurrentMB * 1024 * 1024);
		SwAligner sw(NULL), osw(NULL);
		SeedResults shs[2];
		ReportingMetrics rpm;
		RandomSource rnd;

		int pepolFlag;
		if(gMate1fw && gMate2fw) {
			pepolFlag = PE_POLICY_FF;
		} else if(gMate1fw && !gMate2fw) {
			pepolFlag = PE_POLICY_FR;
		} else if(!gMate1fw && gMate2fw) {
			pepolFlag = PE_POLICY_RF;
		} else {
			pepolFlag = PE_POLICY_RR;
		}
		assert_geq(gMaxInsert, gMinInsert);
		assert_geq(gMinInsert, 0);
		PairedEndPolicy pepol(
			pepolFlag,
			gMaxInsert,
			gMinInsert,
			localAlign,
			gFlippedMatesOK,
			gDovetailMatesOK,
			gContainMatesOK,
			gOlapMatesOK,
			gExpandToFrag);

		EList<Seed> seeds1, seeds2;
		EList<Seed> *seeds[2] = { &seeds1, &seeds2 };

		PerReadMetrics prm;

		// Used by thread with threadid == 1 to measure time elapsed
		time_t iTime = time(0);

		// Keep track of whether last search was exhaustive for mates 1 and 2
		bool exhaustive[2] = { false, false };
		// Keep track of whether mates 1/2 were filtered out last time through
		bool filt[2]    = { true, true };
		// Keep track of whether mates 1/2 were filtered out due Ns last time
		bool nfilt[2]   = { true, true };
		// Keep track of whether mates 1/2 were filtered out due to not having
		// enough characters to rise about the score threshold.
		bool scfilt[2]  = { true, true };
		// Keep track of whether mates 1/2 were filtered out due to not having
		// more characters than the number of mismatches permitted in a seed.
		bool lenfilt[2] = { true, true };
		// Keep track of whether mates 1/2 were filtered out by upstream qc
		bool qcfilt[2]  = { true, true };

                std::unique_ptr<PatternSourceReadAhead> g_psrah(new PatternSourceReadAhead(readahead_factory));
                do { //while have_next_read(g_psrah)
                        pair<bool, bool> ret = g_psrah.get()->nextReadPair();
			{
			  bool success = ret.first;
			  bool done = ret.second;
			  if(!success && done) {
				break;
			  } else if(!success) {
				continue;
			  }
			}

			PatternSourcePerThread* const ps = g_psrah.get()->ptr();
			TReadId rdid = ps->read_a().rdid;

				// Align this read/pair
				//
				// Check if there is metrics reporting for us to do.
				//
				prm.reset(); // per-read metrics
				prm.doFmString = false;

					ca.nextRead(); // clear the cache
					assert(!ca.aligning());
					bool paired = !ps->read_b().empty();
					const size_t rdlen1 = ps->read_a().length();
					const size_t rdlen2 = paired ? ps->read_b().length() : 0;
					msinkwrap.nextRead(
						&ps->read_a(),
						paired ? &ps->read_b() : NULL,
						rdid,
						sc.qualitiesMatter());
					assert(msinkwrap.inited());
					size_t rdlens[2] = { rdlen1, rdlen2 };
					size_t rdrows[2] = { rdlen1, rdlen2 };
					// Calculate the minimum valid score threshold for the read
					TAlScore minsc[2];
					minsc[0] = minsc[1] = std::numeric_limits<TAlScore>::max();
					{
						minsc[0] = scoreMin.f<TAlScore>(rdlens[0]);
						if(paired) minsc[1] = scoreMin.f<TAlScore>(rdlens[1]);
						{
							if(minsc[0] > 0) {
								if(!gQuiet) printEEScoreMsg(*ps, paired, true);
								minsc[0] = 0;
							}
							if(paired && minsc[1] > 0) {
								if(!gQuiet) printEEScoreMsg(*ps, paired, false);
								minsc[1] = 0;
							}
						}
					}
					// N filter; does the read have too many Ns?
					size_t readns[2] = {0, 0};
					sc.nFilterPair(
						&ps->read_a().patFw,
						paired ? &ps->read_b().patFw : NULL,
						readns[0],
						readns[1],
						nfilt[0],
						nfilt[1]);
					// Score filter; does the read enough character to rise above
					// the score threshold?
					scfilt[0] = sc.scoreFilter(minsc[0], rdlens[0]);
					scfilt[1] = sc.scoreFilter(minsc[1], rdlens[1]);
					lenfilt[0] = lenfilt[1] = true;
					if(rdlens[0] <= (size_t)multiseedMms || rdlens[0] < 2) {
						if(!gQuiet) printMmsSkipMsg(*ps, paired, true, multiseedMms);
						lenfilt[0] = false;
					}
					if((rdlens[1] <= (size_t)multiseedMms || rdlens[1] < 2) && paired) {
						if(!gQuiet) printMmsSkipMsg(*ps, paired, false, multiseedMms);
						lenfilt[1] = false;
					}
					if(rdlens[0] < 2) {
						if(!gQuiet) printLenSkipMsg(*ps, paired, true);
						lenfilt[0] = false;
					}
					if(rdlens[1] < 2 && paired) {
						if(!gQuiet) printLenSkipMsg(*ps, paired, false);
						lenfilt[1] = false;
					}
					qcfilt[0] = qcfilt[1] = true;
					if(qcFilter) {
						qcfilt[0] = (ps->read_a().filter != '0');
						qcfilt[1] = (ps->read_b().filter != '0');
					}
					filt[0] = (nfilt[0] && scfilt[0] && lenfilt[0] && qcfilt[0]);
					filt[1] = (nfilt[1] && scfilt[1] && lenfilt[1] && qcfilt[1]);
					prm.nFilt += (filt[0] ? 0 : 1) + (filt[1] ? 0 : 1);
					Read* rds[2] = { &ps->read_a(), &ps->read_b() };
					// For each mate...
					assert(msinkwrap.empty());
					sd.nextRead(paired, rdrows[0], rdrows[1]); // SwDriver
					size_t minedfw[2] = { 0, 0 };
					size_t minedrc[2] = { 0, 0 };
					// Calcualte nofw / no rc
					bool nofw[2] = { false, false };
					bool norc[2] = { false, false };
					nofw[0] = paired ? (gMate1fw ? gNofw : gNorc) : gNofw;
					norc[0] = paired ? (gMate1fw ? gNorc : gNofw) : gNorc;
					nofw[1] = paired ? (gMate2fw ? gNofw : gNorc) : gNofw;
					norc[1] = paired ? (gMate2fw ? gNorc : gNofw) : gNorc;
					// Calculate nceil
					int nceil[2] = { 0, 0 };
					nceil[0] = nCeil.f<int>((double)rdlens[0]);
					nceil[0] = min(nceil[0], (int)rdlens[0]);
					if(paired) {
						nceil[1] = nCeil.f<int>((double)rdlens[1]);
						nceil[1] = min(nceil[1], (int)rdlens[1]);
					}
					exhaustive[0] = exhaustive[1] = false;
					size_t matemap[2] = { 0, 1 };
					bool pairPostFilt = filt[0] && filt[1];
					if(pairPostFilt) {
						rnd.init(ps->read_a().seed ^ ps->read_b().seed);
					} else {
						rnd.init(ps->read_a().seed);
					}
					// Calculate interval length for both mates
					int interval[2] = { 0, 0 };
					for(size_t mate = 0; mate < (paired ? 2:1); mate++) {
						interval[mate] = msIval.f<int>((double)rdlens[mate]);
						if(filt[0] && filt[1]) {
							// Boost interval length by 20% for paired-end reads
							interval[mate] = (int)(interval[mate] * 1.2 + 0.5);
						}
						interval[mate] = max(interval[mate], 1);
					}
					// Calculate streak length
					size_t streak[2]    = { maxDpStreak,   maxDpStreak };
					size_t mtStreak[2]  = { maxMateStreak, maxMateStreak };
					size_t mxDp[2]      = { maxDp,         maxDp       };
					size_t mxUg[2]      = { maxUg,         maxUg       };
					size_t mxIter[2]    = { maxIters,      maxIters    };
					if(allHits) {
						streak[0]   = streak[1]   = std::numeric_limits<size_t>::max();
						mtStreak[0] = mtStreak[1] = std::numeric_limits<size_t>::max();
						mxDp[0]     = mxDp[1]     = std::numeric_limits<size_t>::max();
						mxUg[0]     = mxUg[1]     = std::numeric_limits<size_t>::max();
						mxIter[0]   = mxIter[1]   = std::numeric_limits<size_t>::max();
					} else if(khits > 1) {
						for(size_t mate = 0; mate < 2; mate++) {
							streak[mate]   += (khits-1) * maxStreakIncr;
							mtStreak[mate] += (khits-1) * maxStreakIncr;
							mxDp[mate]     += (khits-1) * maxItersIncr;
							mxUg[mate]     += (khits-1) * maxItersIncr;
							mxIter[mate]   += (khits-1) * maxItersIncr;
						}
					}
					if(filt[0] && filt[1]) {
						streak[0] = (size_t)ceil((double)streak[0] / 2.0);
						streak[1] = (size_t)ceil((double)streak[1] / 2.0);
						assert_gt(streak[1], 0);
					}
					prm.maxDPFails = streak[0];
					assert_gt(streak[0], 0);
					// Calculate # seed rounds for each mate
					size_t nrounds[2] = { nSeedRounds, nSeedRounds };
					if(filt[0] && filt[1]) {
						nrounds[0] = (size_t)ceil((double)nrounds[0] / 2.0);
						nrounds[1] = (size_t)ceil((double)nrounds[1] / 2.0);
						assert_gt(nrounds[1], 0);
					}
					assert_gt(nrounds[0], 0);
					// Increment counters according to what got filtered
					for(size_t mate = 0; mate < (paired ? 2:1); mate++) {
						if(!filt[mate]) {
							// Mate was rejected by N filter
						} else {
							shs[mate].clear();
							shs[mate].nextRead(mate == 0 ? ps->read_a() : ps->read_b());
							assert(shs[mate].empty());
						}
					}
					const size_t eePeEeltLimit = std::numeric_limits<size_t>::max();
					// Whether we're done with mate1 / mate2
					bool done[2] = { !filt[0], !filt[1] };
					size_t nelt[2] = {0, 0};

					// Find end-to-end exact alignments for each read
					{
						for(size_t matei = 0; matei < (paired ? 2:1); matei++) {
							size_t mate = matemap[matei];
							if(!filt[mate] || done[mate] || msinkwrap.state().doneWithMate(mate == 0)) {
								continue;
							}
							nelt[mate] = al.exactSweep(
								ebwtFw,        // index
								*rds[mate],    // read
								sc,            // scoring scheme
								nofw[mate],    // nofw?
								norc[mate],    // norc?
								2,             // max # edits we care about
								minedfw[mate], // minimum # edits for fw mate
								minedrc[mate], // minimum # edits for rc mate
								true,          // report 0mm hits
								shs[mate]);     // put end-to-end results here
						}
						matemap[0] = 0; matemap[1] = 1;
						if(nelt[0] > 0 && nelt[1] > 0 && nelt[0] > nelt[1]) {
							// Do the mate with fewer exact hits first
							// TODO: Consider mates & orientations separately?
							matemap[0] = 1; matemap[1] = 0;
						}
						for(size_t matei = 0; matei < 2; matei++) {
							size_t mate = matemap[matei];
							if(nelt[mate] == 0 || nelt[mate] > eePeEeltLimit) {
								shs[mate].clearExactE2eHits();
								continue;
							}
							if(msinkwrap.state().doneWithMate(mate == 0)) {
								shs[mate].clearExactE2eHits();
								done[mate] = true;
								continue;
							}
							assert(filt[mate]);
							assert(matei == 0 || paired);
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							int ret = 0;
							if(paired) {
								// Paired-end dynamic programming driver
								ret = sd.extendSeedsPaired(
									*rds[mate],     // mate to align as anchor
									*rds[mate ^ 1], // mate to align as opp.
									mate == 0,      // anchor is mate 1?
									!filt[mate ^ 1],// opposite mate filtered out?
									shs[mate],      // seed hits for anchor
									ebwtFw,         // bowtie index
									ebwtBw,         // rev bowtie index
									ref,            // packed reference strings
									sw,             // dyn prog aligner, anchor
									osw,            // dyn prog aligner, opposite
									sc,             // scoring scheme
									pepol,          // paired-end policy
									-1,             // # mms allowed in a seed
									0,              // length of a seed
									0,              // interval between seeds
									minsc[mate],    // min score for anchor
									minsc[mate^1],  // min score for opp.
									nceil[mate],    // N ceil for anchor
									nceil[mate^1],  // N ceil for opp.
									nofw[mate],     // don't align forward read
									norc[mate],     // don't align revcomp read
									maxhalf,        // max width on one DP side
									mxIter[mate],   // max extend loop iters
									mxUg[mate],     // max # ungapped extends
									mxDp[mate],     // max # DPs
									streak[mate],   // stop after streak of this many end-to-end fails
									streak[mate],   // stop after streak of this many ungap fails
									streak[mate],   // stop after streak of this many dp fails
									mtStreak[mate], // max mate fails per seed range
									doExtend,       // extend seed hits
									enable8,        // use 8-bit SSE where possible
									cminlen,        // checkpoint if read is longer
									cpow2,          // checkpointer interval, log2
									doTri,          // triangular mini-fills?
									tighten,        // -M score tightening mode
									ca,             // seed alignment cache
									rnd,            // pseudo-random source
									prm,            // per-read metrics
									&msinkwrap,     // for organizing hits
									true,           // seek mate immediately
									true,           // report hits once found
									gReportDiscordant,// look for discordant alns?
									gReportMixed,   // look for unpaired alns?
									exhaustive[mate]);
								// Might be done, but just with this mate
							} else {
								// Unpaired dynamic programming driver
								ret = sd.extendSeeds(
									*rds[mate],     // read
									mate == 0,      // mate #1?
									shs[mate],      // seed hits
									ebwtFw,         // bowtie index
									ebwtBw,         // rev bowtie index
									ref,            // packed reference strings
									sw,             // dynamic prog aligner
									sc,             // scoring scheme
									-1,             // # mms allowed in a seed
									0,              // length of a seed
									0,              // interval between seeds
									minsc[mate],    // minimum score for valid
									nceil[mate],    // N ceil for anchor
									maxhalf,        // max width on one DP side
									mxIter[mate],   // max extend loop iters
									mxUg[mate],     // max # ungapped extends
									mxDp[mate],     // max # DPs
									streak[mate],   // stop after streak of this many end-to-end fails
									streak[mate],   // stop after streak of this many ungap fails
									doExtend,       // extend seed hits
									enable8,        // use 8-bit SSE where possible
									cminlen,        // checkpoint if read is longer
									cpow2,          // checkpointer interval, log2
									doTri,          // triangular mini-fills
									tighten,        // -M score tightening mode
									ca,             // seed alignment cache
									rnd,            // pseudo-random source
									prm,            // per-read metrics
									&msinkwrap,     // for organizing hits
									true,           // report hits once found
									exhaustive[mate]);
							}
							assert_gt(ret, 0);
							// Clear out the exact hits so that we don't try to
							// extend them again later!
							shs[mate].clearExactE2eHits();
							if(ret == EXTEND_EXHAUSTED_CANDIDATES) {
								// Not done yet
							} else if(ret == EXTEND_POLICY_FULFILLED) {
								// Policy is satisfied for this mate at least
								if(msinkwrap.state().doneWithMate(mate == 0)) {
									done[mate] = true;
								}
								if(msinkwrap.state().doneWithMate(mate == 1)) {
									done[mate^1] = true;
								}
							} else if(ret == EXTEND_PERFECT_SCORE) {
								// We exhausted this mode at least
								done[mate] = true;
							} else if(ret == EXTEND_EXCEEDED_HARD_LIMIT) {
								// We exceeded a per-read limit
								done[mate] = true;
							} else if(ret == EXTEND_EXCEEDED_SOFT_LIMIT) {
								// Not done yet
							} else {
								//
								cerr << "Bad return value: " << ret << endl;
								// We do not have clean exception handling, so just report to the user 
								done[mate] = true;
							}
							if(!done[mate]) {
								TAlScore perfectScore = sc.perfectScore(rdlens[mate]);
								if(!done[mate] && minsc[mate] == perfectScore) {
									done[mate] = true;
								}
							}
						}
					}

					// 1-mismatch
					{
						for(size_t matei = 0; matei < (paired ? 2:1); matei++) {
							size_t mate = matemap[matei];
							if(!filt[mate] || done[mate] || nelt[mate] > eePeEeltLimit) {
								// Done with this mate
								shs[mate].clear1mmE2eHits();
								nelt[mate] = 0;
								continue;
							}
							nelt[mate] = 0;
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							//rnd.init(ROTL(rds[mate]->seed, 10));
							assert(shs[mate].empty());
							assert(shs[mate].repOk(&ca.current()));
							bool yfw = minedfw[mate] <= 1 && !nofw[mate];
							bool yrc = minedrc[mate] <= 1 && !norc[mate];
							if(yfw || yrc) {
								// Clear out the exact hits
								al.oneMmSearch(
									&ebwtFw,        // BWT index
									ebwtBw,         // BWT' index
									*rds[mate],     // read
									sc,             // scoring scheme
									minsc[mate],    // minimum score
									!yfw,           // don't align forward read
									!yrc,           // don't align revcomp read
									false,          // do exact match
									true,           // do 1mm
									shs[mate]);     // seed hits (hits installed here)
								nelt[mate] = shs[mate].num1mmE2eHits();
							}
						}
						// Possibly reorder the mates
						matemap[0] = 0; matemap[1] = 1;
						if(nelt[0] > 0 && nelt[1] > 0 && nelt[0] > nelt[1]) {
							// Do the mate with fewer exact hits first
							// TODO: Consider mates & orientations separately?
							matemap[0] = 1; matemap[1] = 0;
						}
						for(size_t matei = 0; matei < 2; matei++) {
							size_t mate = matemap[matei];
							if(nelt[mate] == 0 || nelt[mate] > eePeEeltLimit) {
								continue;
							}
							if(msinkwrap.state().doneWithMate(mate == 0)) {
								done[mate] = true;
								continue;
							}
							int ret = 0;
							if(paired) {
								// Paired-end dynamic programming driver
								ret = sd.extendSeedsPaired(
									*rds[mate],     // mate to align as anchor
									*rds[mate ^ 1], // mate to align as opp.
									mate == 0,      // anchor is mate 1?
									!filt[mate ^ 1],// opposite mate filtered out?
									shs[mate],      // seed hits for anchor
									ebwtFw,         // bowtie index
									ebwtBw,         // rev bowtie index
									ref,            // packed reference strings
									sw,             // dyn prog aligner, anchor
									osw,            // dyn prog aligner, opposite
									sc,             // scoring scheme
									pepol,          // paired-end policy
									-1,             // # mms allowed in a seed
									0,              // length of a seed
									0,              // interval between seeds
									minsc[mate],    // min score for anchor
									minsc[mate^1],  // min score for opp.
									nceil[mate],    // N ceil for anchor
									nceil[mate^1],  // N ceil for opp.
									nofw[mate],     // don't align forward read
									norc[mate],     // don't align revcomp read
									maxhalf,        // max width on one DP side
									mxIter[mate],   // max extend loop iters
									mxUg[mate],     // max # ungapped extends
									mxDp[mate],     // max # DPs
									streak[mate],   // stop after streak of this many end-to-end fails
									streak[mate],   // stop after streak of this many ungap fails
									streak[mate],   // stop after streak of this many dp fails
									mtStreak[mate], // max mate fails per seed range
									doExtend,       // extend seed hits
									enable8,        // use 8-bit SSE where possible
									cminlen,        // checkpoint if read is longer
									cpow2,          // checkpointer interval, log2
									doTri,          // triangular mini-fills?
									tighten,        // -M score tightening mode
									ca,             // seed alignment cache
									rnd,            // pseudo-random source
									prm,            // per-read metrics
									&msinkwrap,     // for organizing hits
									true,           // seek mate immediately
									true,           // report hits once found
									gReportDiscordant,// look for discordant alns?
									gReportMixed,   // look for unpaired alns?
									exhaustive[mate]);
								// Might be done, but just with this mate
							} else {
								// Unpaired dynamic programming driver
								ret = sd.extendSeeds(
									*rds[mate],     // read
									mate == 0,      // mate #1?
									shs[mate],      // seed hits
									ebwtFw,         // bowtie index
									ebwtBw,         // rev bowtie index
									ref,            // packed reference strings
									sw,             // dynamic prog aligner
									sc,             // scoring scheme
									-1,             // # mms allowed in a seed
									0,              // length of a seed
									0,              // interval between seeds
									minsc[mate],    // minimum score for valid
									nceil[mate],    // N ceil for anchor
									maxhalf,        // max width on one DP side
									mxIter[mate],   // max extend loop iters
									mxUg[mate],     // max # ungapped extends
									mxDp[mate],     // max # DPs
									streak[mate],   // stop after streak of this many end-to-end fails
									streak[mate],   // stop after streak of this many ungap fails
									doExtend,       // extend seed hits
									enable8,        // use 8-bit SSE where possible
									cminlen,        // checkpoint if read is longer
									cpow2,          // checkpointer interval, log2
									doTri,          // triangular mini-fills?
									tighten,        // -M score tightening mode
									ca,             // seed alignment cache
									rnd,            // pseudo-random source
									prm,            // per-read metrics
									&msinkwrap,     // for organizing hits
									true,           // report hits once found
									exhaustive[mate]);
							}
							assert_gt(ret, 0);
							// Clear out the 1mm hits so that we don't try to
							// extend them again later!
							shs[mate].clear1mmE2eHits();
							if(ret == EXTEND_EXHAUSTED_CANDIDATES) {
								// Not done yet
							} else if(ret == EXTEND_POLICY_FULFILLED) {
								// Policy is satisfied for this mate at least
								if(msinkwrap.state().doneWithMate(mate == 0)) {
									done[mate] = true;
								}
								if(msinkwrap.state().doneWithMate(mate == 1)) {
									done[mate^1] = true;
								}
							} else if(ret == EXTEND_PERFECT_SCORE) {
								// We exhausted this mode at least
								done[mate] = true;
							} else if(ret == EXTEND_EXCEEDED_HARD_LIMIT) {
								// We exceeded a per-read limit
								done[mate] = true;
							} else if(ret == EXTEND_EXCEEDED_SOFT_LIMIT) {
								// Not done yet
							} else {
								//
								cerr << "Bad return value: " << ret << endl;
								// We do not have clean exception handling, so just report to the user 
								done[mate] = true;
							}
							if(!done[mate]) {
								TAlScore perfectScore = sc.perfectScore(rdlens[mate]);
								if(!done[mate] && minsc[mate] == perfectScore) {
									done[mate] = true;
								}
							}
						}
					}
					int seedlens[2] = { multiseedLen, multiseedLen };
					nrounds[0] = min<size_t>(nrounds[0], interval[0]);
					nrounds[1] = min<size_t>(nrounds[1], interval[1]);
					Constraint gc = Constraint::penaltyFuncBased(scoreMin);
					size_t seedsTried = 0;
					size_t seedsTriedMS[] = {0, 0, 0, 0};
					size_t nUniqueSeeds = 0, nRepeatSeeds = 0, seedHitTot = 0;
					size_t nUniqueSeedsMS[] = {0, 0, 0, 0};
					size_t nRepeatSeedsMS[] = {0, 0, 0, 0};
					size_t seedHitTotMS[] = {0, 0, 0, 0};
					for(size_t roundi = 0; roundi < nSeedRounds; roundi++) {
						ca.nextRead(); // Clear cache in preparation for new search
						shs[0].clearSeeds();
						shs[1].clearSeeds();
						assert(shs[0].empty());
						assert(shs[1].empty());
						assert(shs[0].repOk(&ca.current()));
						assert(shs[1].repOk(&ca.current()));
						//if(roundi > 0) {
						//	if(seedlens[0] > 8) seedlens[0]--;
						//	if(seedlens[1] > 8) seedlens[1]--;
						//}
						for(size_t matei = 0; matei < (paired ? 2:1); matei++) {
							size_t mate = matemap[matei];
							if(done[mate] || msinkwrap.state().doneWithMate(mate == 0)) {
								// Done with this mate
								done[mate] = true;
								continue;
							}
							if(roundi >= nrounds[mate]) {
								// Not doing this round for this mate
								continue;
							}
							// Figure out the seed offset
							if(interval[mate] <= (int)roundi) {
								// Can't do this round, seeds already packed as
								// tight as possible
								continue;
							}
							size_t offset = (interval[mate] * roundi) / nrounds[mate];
							assert(roundi == 0 || offset > 0);
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							//rnd.init(ROTL(rds[mate]->seed, 10));
							assert(shs[mate].repOk(&ca.current()));
							// Set up seeds
							seeds[mate]->clear();
							Seed::mmSeeds(
								multiseedMms,    // max # mms per seed
								seedlens[mate],  // length of a multiseed seed
								*seeds[mate],    // seeds
								gc);             // global constraint
							// Check whether the offset would drive the first seed
							// off the end
							if(offset > 0 && (*seeds[mate])[0].len + offset > rds[mate]->length()) {
								continue;
							}
							// Instantiate the seeds
							std::pair<int, int> instFw, instRc;
							std::pair<int, int> inst = al.instantiateSeeds(
								*seeds[mate],   // search seeds
								offset,         // offset to begin extracting
								interval[mate], // interval between seeds
								*rds[mate],     // read to align
								sc,             // scoring scheme
								nofw[mate],     // don't align forward read
								norc[mate],     // don't align revcomp read
								ca,             // holds some seed hits from previous reads
								shs[mate],      // holds all the seed hits
								instFw,
								instRc);
							assert(shs[mate].repOk(&ca.current()));
							if(inst.first + inst.second == 0) {
								// No seed hits!  Done with this mate.
								assert(shs[mate].empty());
								done[mate] = true;
								break;
							}
							seedsTried += (inst.first + inst.second);
							seedsTriedMS[mate * 2 + 0] = instFw.first + instFw.second;
							seedsTriedMS[mate * 2 + 1] = instRc.first + instRc.second;
							// Align seeds
							al.searchAllSeeds(
								*seeds[mate],     // search seeds
								&ebwtFw,          // BWT index
								ebwtBw,           // BWT' index
								*rds[mate],       // read
								sc,               // scoring scheme
								ca,               // alignment cache
								shs[mate],        // store seed hits here
								prm);             // per-read metrics
							assert(shs[mate].repOk(&ca.current()));
							if(shs[mate].empty()) {
								// No seed alignments!  Done with this mate.
								done[mate] = true;
								break;
							}
						}
						// shs contain what we need to know to update our seed
						// summaries for this seeding
						for(size_t mate = 0; mate < 2; mate++) {
							if(!shs[mate].empty()) {
								nUniqueSeeds += shs[mate].numUniqueSeeds();
								nUniqueSeedsMS[mate * 2 + 0] += shs[mate].numUniqueSeedsStrand(true);
								nUniqueSeedsMS[mate * 2 + 1] += shs[mate].numUniqueSeedsStrand(false);
								nRepeatSeeds += shs[mate].numRepeatSeeds();
								nRepeatSeedsMS[mate * 2 + 0] += shs[mate].numRepeatSeedsStrand(true);
								nRepeatSeedsMS[mate * 2 + 1] += shs[mate].numRepeatSeedsStrand(false);
								seedHitTot += shs[mate].numElts();
								seedHitTotMS[mate * 2 + 0] += shs[mate].numEltsFw();
								seedHitTotMS[mate * 2 + 1] += shs[mate].numEltsRc();
							}
						}
						double uniqFactor[2] = { 0.0f, 0.0f };
						for(size_t i = 0; i < 2; i++) {
							if(!shs[i].empty()) {
								uniqFactor[i] = shs[i].uniquenessFactor();
							}
						}
						// Possibly reorder the mates
						matemap[0] = 0; matemap[1] = 1;
						if(!shs[0].empty() && !shs[1].empty() && uniqFactor[1] > uniqFactor[0]) {
							// Do the mate with fewer exact hits first
							// TODO: Consider mates & orientations separately?
							matemap[0] = 1; matemap[1] = 0;
						}
						for(size_t matei = 0; matei < (paired ? 2:1); matei++) {
							size_t mate = matemap[matei];
							if(done[mate] || msinkwrap.state().doneWithMate(mate == 0)) {
								// Done with this mate
								done[mate] = true;
								continue;
							}
							assert(!msinkwrap.maxed());
							assert(msinkwrap.repOk());
							//rnd.init(ROTL(rds[mate]->seed, 10));
							assert(shs[mate].repOk(&ca.current()));
							{
								// If there aren't any seed hits...
								if(shs[mate].empty()) {
									continue; // on to the next mate
								}
								// Sort seed hits into ranks
								shs[mate].rankSeedHits(rnd, msinkwrap.allHits());
								int ret = 0;
								if(paired) {
									// Paired-end dynamic programming driver
									ret = sd.extendSeedsPaired(
										*rds[mate],     // mate to align as anchor
										*rds[mate ^ 1], // mate to align as opp.
										mate == 0,      // anchor is mate 1?
										!filt[mate ^ 1],// opposite mate filtered out?
										shs[mate],      // seed hits for anchor
										ebwtFw,         // bowtie index
										ebwtBw,         // rev bowtie index
										ref,            // packed reference strings
										sw,             // dyn prog aligner, anchor
										osw,            // dyn prog aligner, opposite
										sc,             // scoring scheme
										pepol,          // paired-end policy
										multiseedMms,   // # mms allowed in a seed
										seedlens[mate], // length of a seed
										interval[mate], // interval between seeds
										minsc[mate],    // min score for anchor
										minsc[mate^1],  // min score for opp.
										nceil[mate],    // N ceil for anchor
										nceil[mate^1],  // N ceil for opp.
										nofw[mate],     // don't align forward read
										norc[mate],     // don't align revcomp read
										maxhalf,        // max width on one DP side
										mxIter[mate],   // max extend loop iters
										mxUg[mate],     // max # ungapped extends
										mxDp[mate],     // max # DPs
										streak[mate],   // stop after streak of this many end-to-end fails
										streak[mate],   // stop after streak of this many ungap fails
										streak[mate],   // stop after streak of this many dp fails
										mtStreak[mate], // max mate fails per seed range
										doExtend,       // extend seed hits
										enable8,        // use 8-bit SSE where possible
										cminlen,        // checkpoint if read is longer
										cpow2,          // checkpointer interval, log2
										doTri,          // triangular mini-fills?
										tighten,        // -M score tightening mode
										ca,             // seed alignment cache
										rnd,            // pseudo-random source
										prm,            // per-read metrics
										&msinkwrap,     // for organizing hits
										true,           // seek mate immediately
										true,           // report hits once found
										gReportDiscordant,// look for discordant alns?
										gReportMixed,   // look for unpaired alns?
										exhaustive[mate]);
									// Might be done, but just with this mate
								} else {
									// Unpaired dynamic programming driver
									ret = sd.extendSeeds(
										*rds[mate],     // read
										mate == 0,      // mate #1?
										shs[mate],      // seed hits
										ebwtFw,         // bowtie index
										ebwtBw,         // rev bowtie index
										ref,            // packed reference strings
										sw,             // dynamic prog aligner
										sc,             // scoring scheme
										multiseedMms,   // # mms allowed in a seed
										seedlens[mate], // length of a seed
										interval[mate], // interval between seeds
										minsc[mate],    // minimum score for valid
										nceil[mate],    // N ceil for anchor
										maxhalf,        // max width on one DP side
										mxIter[mate],   // max extend loop iters
										mxUg[mate],     // max # ungapped extends
										mxDp[mate],     // max # DPs
										streak[mate],   // stop after streak of this many end-to-end fails
										streak[mate],   // stop after streak of this many ungap fails
										doExtend,       // extend seed hits
										enable8,        // use 8-bit SSE where possible
										cminlen,        // checkpoint if read is longer
										cpow2,          // checkpointer interval, log2
										doTri,          // triangular mini-fills?
										tighten,        // -M score tightening mode
										ca,             // seed alignment cache
										rnd,            // pseudo-random source
										prm,            // per-read metrics
										&msinkwrap,     // for organizing hits
										true,           // report hits once found
										exhaustive[mate]);
								}
								assert_gt(ret, 0);
								if(ret == EXTEND_EXHAUSTED_CANDIDATES) {
									// Not done yet
								} else if(ret == EXTEND_POLICY_FULFILLED) {
									// Policy is satisfied for this mate at least
									if(msinkwrap.state().doneWithMate(mate == 0)) {
										done[mate] = true;
									}
									if(msinkwrap.state().doneWithMate(mate == 1)) {
										done[mate^1] = true;
									}
								} else if(ret == EXTEND_PERFECT_SCORE) {
									// We exhausted this made at least
									done[mate] = true;
								} else if(ret == EXTEND_EXCEEDED_HARD_LIMIT) {
									// We exceeded a per-read limit
									done[mate] = true;
								} else if(ret == EXTEND_EXCEEDED_SOFT_LIMIT) {
									// Not done yet
								} else {
									//
									cerr << "Bad return value: " << ret << endl;
									// We do not have clean exception handling, so just report to the user 
									done[mate] = true;
								}
							}
						} // for(size_t matei = 0; matei < 2; matei++)

						// We don't necessarily have to continue investigating both
						// mates.  We continue on a mate only if its average
						// interval length is high (> 1000)
						for(size_t mate = 0; mate < 2; mate++) {
							if(!done[mate] && shs[mate].averageHitsPerSeed() < seedBoostThresh) {
								done[mate] = true;
							}
						}
					} // end loop over reseeding rounds
					if(seedsTried > 0) {
						prm.seedPctUnique = (float)nUniqueSeeds / seedsTried;
						prm.seedPctRep = (float)nRepeatSeeds / seedsTried;
						prm.seedHitAvg = (float)seedHitTot / seedsTried;
					} else {
						prm.seedPctUnique = -1.0f;
						prm.seedPctRep = -1.0f;
						prm.seedHitAvg = -1.0f;
					}
					for(int i = 0; i < 4; i++) {
						if(seedsTriedMS[i] > 0) {
							prm.seedPctUniqueMS[i] = (float)nUniqueSeedsMS[i] / seedsTriedMS[i];
							prm.seedPctRepMS[i] = (float)nRepeatSeedsMS[i] / seedsTriedMS[i];
							prm.seedHitAvgMS[i] = (float)seedHitTotMS[i] / seedsTriedMS[i];
						} else {
							prm.seedPctUniqueMS[i] = -1.0f;
							prm.seedPctRepMS[i] = -1.0f;
							prm.seedHitAvgMS[i] = -1.0f;
						}
					}
					size_t totnucs = 0;
					for(size_t mate = 0; mate < (paired ? 2:1); mate++) {
						if(filt[mate]) {
							size_t len = rdlens[mate];
							if(!nofw[mate] && !norc[mate]) {
								len *= 2;
							}
							totnucs += len;
						}
					}
					prm.seedsPerNuc = totnucs > 0 ? ((float)seedsTried / totnucs) : -1;
					for(int i = 0; i < 4; i++) {
						prm.seedsPerNucMS[i] = totnucs > 0 ? ((float)seedsTriedMS[i] / totnucs) : -1;
					}
					for(size_t i = 0; i < 2; i++) {
						assert_leq(prm.nExIters, mxIter[i]);
						assert_leq(prm.nExDps,   mxDp[i]);
						assert_leq(prm.nMateDps, mxDp[i]);
						assert_leq(prm.nExUgs,   mxUg[i]);
						assert_leq(prm.nMateUgs, mxUg[i]);
						assert_leq(prm.nDpFail,  streak[i]);
						assert_leq(prm.nUgFail,  streak[i]);
						assert_leq(prm.nEeFail,  streak[i]);
					}

				// Commit and report paired-end/unpaired alignments
				//uint32_t sd = rds[0]->seed ^ rds[1]->seed;
				//rnd.init(ROTL(sd, 20));
				msinkwrap.finishRead(
					&shs[0],              // seed results for mate 1
					&shs[1],              // seed results for mate 2
					exhaustive[0],        // exhausted seed hits for mate 1?
					exhaustive[1],        // exhausted seed hits for mate 2?
					nfilt[0],
					nfilt[1],
					scfilt[0],
					scfilt[1],
					lenfilt[0],
					lenfilt[1],
					qcfilt[0],
					qcfilt[1],
					rnd,                  // pseudo-random generator
					rpm,                  // reporting metrics
					prm,                  // per-read metrics
					sc,                   // scoring scheme
					true,                 // suppress seed summaries?
					false,                // suppress alignments?
					scUnMapped,           // Consider soft-clipped bases unmapped when calculating TLEN
					xeq);

		} while (have_next_read(g_psrah)); // must read the whole cached buffer

		// Merge in the metrics
		msink.mergeMetricsUnsafe(rpm);
	}

	return;
}
#endif

#ifdef ENABLE_2P5

// NOTE: Unsupported, likely does not work

//void multiseedSearchWorker_2p5::operator()() const {
static void multiseedSearchWorker_2p5(const size_t num_parallel_tasks) {
	int tid = 0;
	assert(multiseed_ebwtFw != NULL);
	assert(multiseedMms == 0 || multiseed_ebwtBw != NULL);
	PatternSourceReadAheadFactory& readahead_factory =  *multiseed_readahead_factory;
	const Ebwt&             ebwtFw   = *multiseed_ebwtFw;
	const Ebwt&             ebwtBw   = *multiseed_ebwtBw;
	const Scoring&          sc       = *multiseed_sc;
	const BitPairReference& ref      = *multiseed_refs;
	AlnSink&                msink    = *multiseed_msink;

	// Sinks: these are so that we can print tables encoding counts for
	// events of interest on a per-read, per-seed, per-join, or per-SW
	// level.  These in turn can be used to diagnose performance
	// problems, or generally characterize performance.

	ThreadCounter tc;

	// Instantiate an object for holding reporting-related parameters.
	ReportingParams rp(
		(allHits ? std::numeric_limits<THitInt>::max() : khits), // -k
		mhits,             // -m/-M
		0,                 // penalty gap (not used now)
		msample,           // true -> -M was specified, otherwise assume -m
		gReportDiscordant, // report discordang paired-end alignments?
		gReportMixed);     // report unpaired alignments for paired reads?

	// Instantiate a mapping quality calculator
	unique_ptr<Mapq> bmapq(new_mapq(mapqv, scoreMin, sc));

	// Make a per-thread wrapper for the global MHitSink object.
	AlnSinkWrap msinkwrap(
		msink,         // global sink
		rp,            // reporting parameters
		*bmapq,        // MAPQ calculator
		(size_t)tid);  // thread id

	DescentMetrics descm;
	ReportingMetrics rpm;
	RandomSource rnd;

	int pepolFlag;
	if(gMate1fw && gMate2fw) {
		pepolFlag = PE_POLICY_FF;
	} else if(gMate1fw && !gMate2fw) {
		pepolFlag = PE_POLICY_FR;
	} else if(!gMate1fw && gMate2fw) {
		pepolFlag = PE_POLICY_RF;
	} else {
		pepolFlag = PE_POLICY_RR;
	}
	assert_geq(gMaxInsert, gMinInsert);
	assert_geq(gMinInsert, 0);
	PairedEndPolicy pepol(
		pepolFlag,
		gMaxInsert,
		gMinInsert,
		localAlign,
		gFlippedMatesOK,
		gDovetailMatesOK,
		gContainMatesOK,
		gOlapMatesOK,
		gExpandToFrag);

	AlignerDriver ald(
		descConsExp,         // exponent for interpolating maximum penalty
		descPrioritizeRoots, // whether to select roots with scores and weights
		msIval,              // interval length, as function of read length
		descLanding,         // landing length
		gVerbose,            // verbose?
		descentTotSz,        // limit on total bytes of best-first search data
		descentTotFmops);    // limit on total number of FM index ops in BFS

	PerReadMetrics prm;

	// Used by thread with threadid == 1 to measure time elapsed
	time_t iTime = time(0);

	// Keep track of whether last search was exhaustive for mates 1 and 2
	bool exhaustive[2] = { false, false };
	// Keep track of whether mates 1/2 were filtered out last time through
	bool filt[2]    = { true, true };
	// Keep track of whether mates 1/2 were filtered out due Ns last time
	bool nfilt[2]   = { true, true };
	// Keep track of whether mates 1/2 were filtered out due to not having
	// enough characters to rise about the score threshold.
	bool scfilt[2]  = { true, true };
	// Keep track of whether mates 1/2 were filtered out due to not having
	// more characters than the number of mismatches permitted in a seed.
	bool lenfilt[2] = { true, true };
	// Keep track of whether mates 1/2 were filtered out by upstream qc
	bool qcfilt[2]  = { true, true };

	while(true) {
	   PatternSourceReadAhead psrah(readahead_factory);
	   PatternSourcePerThread* const ps = psrah.ptr();
	   do {
		pair<bool, bool> ret = ps->nextReadPair();
		bool success = ret.first;
		bool done = ret.second;
		if(!success && done) {
			break;
		} else if(!success) {
			continue;
		}
		TReadId rdid = ps->read_a().rdid;

			//
			// Check if there is metrics reporting for us to do.
			//
			prm.reset(); // per-read metrics
			prm.doFmString = sam_print_zm;
			// If we're reporting how long each read takes, get the initial time
			// measurement here
			// Try to align this read
			bool paired = !ps->read_b().empty();
			const size_t rdlen1 = ps->read_a().length();
			const size_t rdlen2 = paired ? ps->read_b().length() : 0;
			// Check if read is identical to previous read
			rnd.init(ROTL(ps->read_a().seed, 5));
			msinkwrap.nextRead(
				&ps->read_a(),
				paired ? &ps->read_b() : NULL,
				rdid,
				sc.qualitiesMatter());
			assert(msinkwrap.inited());
			size_t rdlens[2] = { rdlen1, rdlen2 };
			// Calculate the minimum valid score threshold for the read
			TAlScore minsc[2], maxpen[2];
			minsc[0] = minsc[1] = std::numeric_limits<TAlScore>::max();
			setupMinScores(*ps, paired, sc, rdlens, minsc, maxpen);
			// N filter; does the read have too many Ns?
			size_t readns[2] = {0, 0};
			sc.nFilterPair(
				&ps->read_a().patFw,
				paired ? &ps->read_b().patFw : NULL,
				readns[0],
				readns[1],
				nfilt[0],
				nfilt[1]);
			// Score filter; does the read enough character to rise above
			// the score threshold?
			scfilt[0] = sc.scoreFilter(minsc[0], rdlens[0]);
			scfilt[1] = sc.scoreFilter(minsc[1], rdlens[1]);
			lenfilt[0] = lenfilt[1] = true;
			if(rdlens[0] <= (size_t)multiseedMms || rdlens[0] < 2) {
				if(!gQuiet) printMmsSkipMsg(*ps, paired, true, multiseedMms);
				lenfilt[0] = false;
			}
			if((rdlens[1] <= (size_t)multiseedMms || rdlens[1] < 2) && paired) {
				if(!gQuiet) printMmsSkipMsg(*ps, paired, false, multiseedMms);
				lenfilt[1] = false;
			}
			if(rdlens[0] < 2) {
				if(!gQuiet) printLenSkipMsg(*ps, paired, true);
				lenfilt[0] = false;
			}
			if(rdlens[1] < 2 && paired) {
				if(!gQuiet) printLenSkipMsg(*ps, paired, false);
				lenfilt[1] = false;
			}
			qcfilt[0] = qcfilt[1] = true;
			if(qcFilter) {
				qcfilt[0] = (ps->read_a().filter != '0');
				qcfilt[1] = (ps->read_b().filter != '0');
			}
			filt[0] = (nfilt[0] && scfilt[0] && lenfilt[0] && qcfilt[0]);
			filt[1] = (nfilt[1] && scfilt[1] && lenfilt[1] && qcfilt[1]);
			prm.nFilt += (filt[0] ? 0 : 1) + (filt[1] ? 0 : 1);
			Read* rds[2] = { &ps->read_a(), &ps->read_b() };
			assert(msinkwrap.empty());
			// Calcualte nofw / no rc
			bool nofw[2] = { false, false };
			bool norc[2] = { false, false };
			nofw[0] = paired ? (gMate1fw ? gNofw : gNorc) : gNofw;
			norc[0] = paired ? (gMate1fw ? gNorc : gNofw) : gNorc;
			nofw[1] = paired ? (gMate2fw ? gNofw : gNorc) : gNofw;
			norc[1] = paired ? (gMate2fw ? gNorc : gNofw) : gNorc;
			// Calculate nceil
			int nceil[2] = { 0, 0 };
			nceil[0] = nCeil.f<int>((double)rdlens[0]);
			nceil[0] = min(nceil[0], (int)rdlens[0]);
			if(paired) {
				nceil[1] = nCeil.f<int>((double)rdlens[1]);
				nceil[1] = min(nceil[1], (int)rdlens[1]);
			}
			exhaustive[0] = exhaustive[1] = false;
			bool pairPostFilt = filt[0] && filt[1];
			if(pairPostFilt) {
				rnd.init(ROTL((rds[0]->seed ^ rds[1]->seed), 10));
			}
			// Calculate streak length
			size_t streak[2]    = { maxDpStreak,   maxDpStreak };
			size_t mtStreak[2]  = { maxMateStreak, maxMateStreak };
			size_t mxDp[2]      = { maxDp,         maxDp       };
			size_t mxUg[2]      = { maxUg,         maxUg       };
			size_t mxIter[2]    = { maxIters,      maxIters    };
			if(allHits) {
				streak[0]   = streak[1]   = std::numeric_limits<size_t>::max();
				mtStreak[0] = mtStreak[1] = std::numeric_limits<size_t>::max();
				mxDp[0]     = mxDp[1]     = std::numeric_limits<size_t>::max();
				mxUg[0]     = mxUg[1]     = std::numeric_limits<size_t>::max();
				mxIter[0]   = mxIter[1]   = std::numeric_limits<size_t>::max();
			} else if(khits > 1) {
				for(size_t mate = 0; mate < 2; mate++) {
					streak[mate]   += (khits-1) * maxStreakIncr;
					mtStreak[mate] += (khits-1) * maxStreakIncr;
					mxDp[mate]     += (khits-1) * maxItersIncr;
					mxUg[mate]     += (khits-1) * maxItersIncr;
					mxIter[mate]   += (khits-1) * maxItersIncr;
				}
			}
			// If paired-end and neither mate filtered...
			if(filt[0] && filt[1]) {
				// Reduce streaks for either mate
				streak[0] = (size_t)ceil((double)streak[0] / 2.0);
				streak[1] = (size_t)ceil((double)streak[1] / 2.0);
				assert_gt(streak[1], 0);
			}
			assert_gt(streak[0], 0);
			prm.maxDPFails = streak[0];
			if(filt[0]) {
				ald.initRead(ps->read_a(), nofw[0], norc[0], minsc[0], maxpen[0], filt[1] ? &ps->read_b() : NULL);
			} else if(filt[1]) {
				ald.initRead(ps->read_b(), nofw[1], norc[1], minsc[1], maxpen[1], NULL);
			}
			if(filt[0] || filt[1]) {
				ald.go(sc, ebwtFw, ebwtBw, ref, descm, prm, rnd, msinkwrap);
			}
			// Commit and report paired-end/unpaired alignments
			uint32_t sd = rds[0]->seed ^ rds[1]->seed;
			rnd.init(ROTL(sd, 20));
			msinkwrap.finishRead(
				NULL,                 // seed results for mate 1
				NULL,                 // seed results for mate 2
				exhaustive[0],        // exhausted seed results for 1?
				exhaustive[1],        // exhausted seed results for 2?
				nfilt[0],
				nfilt[1],
				scfilt[0],
				scfilt[1],
				lenfilt[0],
				lenfilt[1],
				qcfilt[0],
				qcfilt[1],
				rnd,                  // pseudo-random generator
				rpm,                  // reporting metrics
				prm,                  // per-read metrics
				sc,                   // scoring scheme
				true,                 // suppress seed summaries?
				false,                // suppress alignments?
				scUnMapped,           // Consider soft-clipped bases unmapped when calculating TLEN
				xeq);

	   } while (ps->nextReadPairReady()); // must read the whole cached buffer
	} // while(true)

	// Merge in the metrics
	msink.mergeMetricsUnsafe(rpm);

	return;
}

#endif

#ifndef _WIN32
/**
 * Print friendly-ish message pertaining to failed system call.
 */
static void errno_message() {
	int errnum = errno;
	cerr << "errno is " << errnum << endl;
	perror("perror error: ");
}

#if 0
/**
 * Delete PID file.  Raise error if the file doesn't exist or if
 * we fail to delete it.
 */
void del_pid(const char* dirname,int pid) {
	char* fname = (char*)calloc(FNAME_SIZE, sizeof(char));
	if(fname == NULL) {
		errno_message();
		cerr << "del_pid: could not allocate buffer" << endl;
		throw 1;
	}
	snprintf(fname, FNAME_SIZE, "%s/%d", dirname, pid);
	if(unlink(fname) != 0) {
		if(errno != ENOENT) {
			errno_message();
			cerr << "del_pid: could not delete PID file " << fname << endl;
			free(fname);
			throw 1;
		} else {
			// Probably just a race between processes
		}
	}
	free(fname);
}

/**
 * Write PID file.
 */
static void write_pid(const char* dirname,int pid) {
	struct stat dinfo;
	if(stat(dirname, &dinfo) != 0) {
		if(mkdir(dirname, 0755) != 0) {
			if(errno != EEXIST) {
				errno_message();
				cerr << "write_pid: could not create PID directory " << dirname << endl;
				throw 1;
			}
		}
	}
	char* fname = (char*)calloc(FNAME_SIZE, sizeof(char));
	if(fname == NULL) {
		errno_message();
		cerr << "write_pid: could not allocate buffer" << endl;
		throw 1;
	}
	snprintf(fname, FNAME_SIZE, "%s/%d", dirname, pid);
	FILE *f = fopen(fname, "w");
	if(f == NULL) {
		errno_message();
		cerr << "write_pid: could not open PID file " << fname << endl;
		throw 1;
	}
	if(fclose(f) != 0) {
		errno_message();
		cerr << "write_pid: could not close PID file " << fname << endl;
		throw 1;
	}
	free(fname);
}

/**
 * Read all the PID files in the given PID directory.  If the
 * process corresponding to a PID file seems to have expired,
 * delete the PID file.  Return the lowest PID encountered for
 * a still-valid process.
 */
static int read_dir(const char* dirname, int* num_pids) {
	DIR *dir;
	struct dirent *ent;
	char* fname = (char*)calloc(FNAME_SIZE, sizeof(char));
	if(fname == NULL) {
		errno_message();
		cerr << "read_dir: could not allocate buffer" << endl;
		throw 1;
	}
	dir = opendir(dirname);
	if(dir == NULL) {
		errno_message();
		cerr << "read_dir: could not open directory " << dirname << endl;
		free(fname);
		throw 1;
	}
	int lowest_pid = -1;
	while(true) {
		errno = 0;
		ent = readdir(dir);
		if(ent == NULL) {
			if(errno != 0) {
				errno_message();
				cerr << "read_dir: could not read directory " << dirname << endl;
				free(fname);
				throw 1;
			}
			break;
		}
		if(ent->d_name[0] == '.') {
			continue;
		}
		int pid = atoi(ent->d_name);
		if(kill(pid, 0) != 0) {
			if(errno == ESRCH) {
				del_pid(dirname, pid);
				continue;
			} else {
				errno_message();
				cerr << "read_dir: could not interrogate pid " << pid << endl;
				free(fname);
				throw 1;
			}
		}
		(*num_pids)++;
		if(pid < lowest_pid || lowest_pid == -1) {
			lowest_pid = pid;
		}
	}
	if(closedir(dir) != 0) {
		errno_message();
		cerr << "read_dir: could not close directory " << dir << endl;
		free(fname);
		throw 1;
	}
	free(fname);
	return lowest_pid;
}
#endif


#endif
/**
 * Called once per alignment job.  Sets up global pointers to the
 * shared global data structures, creates per-thread structures, then
 * enters the search loop.
 */
static void multiseedSearch(
	Scoring& sc,
	const PatternParams& pp,
	PatternComposer& patsrc,      // pattern source
	AlnSink& msink,               // hit sink
	Ebwt& ebwtFw,                 // index of original text
	Ebwt* ebwtBw)                 // index of mirror text
{
	multiseed_msink  = &msink;
	multiseed_ebwtFw = &ebwtFw;
	multiseed_ebwtBw = ebwtBw;
	multiseed_sc     = &sc;
	Timer *_t = new Timer(cerr, "Time loading reference: ", timing);
	unique_ptr<BitPairReference> refs(
		new BitPairReference(
			adjIdxBase,
			false,
			sanityCheck,
			NULL,
			NULL,
			false,
			useMm,
			useShmem,
			mmSweep,
			gVerbose,
			startVerbose)
		);
	delete _t;
	if(!refs->loaded()) throw 1;
	multiseed_refs = refs.get();
#ifndef _WIN32
	sigset_t set;
	sigemptyset(&set);
	sigaddset(&set, SIGPIPE);
	pthread_sigmask(SIG_BLOCK, &set, NULL);
#endif
	{
		// Load the other half of the index into memory
		assert(!ebwtFw.isInMemory());
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(
			0,  // colorspace?
			-1, // not the reverse index
			true,         // load SA samp? (yes, need forward index's SA samp)
			true,         // load ftab (in forward index)
			true,         // load rstarts (in forward index)
			!noRefNames,  // load names?
			startVerbose);
	}
	{
		// Load the other half of the index into memory
		assert(!ebwtBw->isInMemory());
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw->loadIntoMemory(
			0, // colorspace?
			// It's bidirectional search, so we need the reverse to be
			// constructed as the reverse of the concatenated strings.
			1,
			false,        // don't load SA samp in reverse index
			true,         // yes, need ftab in reverse index
			false,        // don't load rstarts in reverse index
			!noRefNames,  // load names?
			startVerbose);
	}

	// Important: Need at least nthreads+1 elements, more is OK
	PatternSourceReadAheadFactory readahead_factory(patsrc,pp,2*nthreads+1);
	multiseed_readahead_factory = &readahead_factory;

	// Start the metrics thread

	std::atomic<int> all_threads_done;
	all_threads_done = 0;
	{
		Timer _t(cerr, "Multiseed full-index search: ", timing);

		{

#ifdef ENABLE_2P5
			if(bowtie2p5) { // WARNING: generally unsupported
				multiseedSearchWorker_2p5(nthreads);
			} else {
#endif
#ifdef ENABLE_PAIRED
			if(paired) { // WARNING: generally unsupported
				multiseedSearchWorkerPaired(nthreads);
			} else {
#endif
				multiseedSearchWorker(nthreads);
#ifdef ENABLE_2P5
			}
#endif
		}

	}
}

static string argstr;

template<typename TStr>
static void driver(
	const char * type,
	const string& bt2indexBase,
	const string& outfile)
{
	if(gVerbose || startVerbose)  {
		cerr << "Entered driver(): "; logTime(cerr, true);
	}
	// Vector of the reference sequences; used for sanity-checking
	EList<SString<char> > names, os;
	EList<size_t> nameLens, seqLens;
	// Read reference sequences from the command-line or from a FASTA file
	if(!origString.empty()) {
		// Read fasta file(s)
		EList<string> origFiles;
		tokenize(origString, ",", origFiles);
		parseFastas(origFiles, names, nameLens, os, seqLens);
	}
	PatternParams pp(
		format,        // file format
		interleaved,   // some or all of the reads are interleaved
		fileParallel,  // true -> wrap files with separate PairedPatternSources
		seed,          // pseudo-random seed
		readsPerBatch, // # reads in a light parsing batch
		solexaQuals,   // true -> qualities are on solexa64 scale
		phred64Quals,  // true -> qualities are on phred64 scale
		integerQuals,  // true -> qualities are space-separated numbers
		gTrim5,        // amt to hard clip from 5' end
		gTrim3,        // amt to hard clip from 3' end
		trimTo,        // trim reads exceeding given length from either 3' or 5'-end
		fastaContLen,  // length of sampled reads for FastaContinuous...
		fastaContFreq, // frequency of sampled reads for FastaContinuous...
		nthreads,      //number of threads for locking
		outType != OUTPUT_SAM, // whether to fix mate names
		preserve_tags, // keep existing tags when aligning BAM files
		align_paired_reads // Align only the paired reads in BAM file
		);
	if(gVerbose || startVerbose) {
		cerr << "Creating PatternSource: "; logTime(cerr, true);
	}
	PatternComposer *patsrc = PatternComposer::setupPatternComposer(
		queries,     // singles, from argv
		mates1,      // mate1's, from -1 arg
		mates2,      // mate2's, from -2 arg
		mates12,     // both mates on each line, from --12 arg
		qualities,   // qualities associated with singles
		qualities1,  // qualities associated with m1
		qualities2,  // qualities associated with m2
#ifdef USE_SRA
		sra_accs,    // SRA accessions
#endif
		pp,          // read read-in parameters
		gVerbose || startVerbose); // be talkative
	// Open hit output file
	if(gVerbose || startVerbose) {
		cerr << "Opening hit output file: "; logTime(cerr, true);
	}
	OutFileBuf *fout;
	if(!outfile.empty()) {
		fout = new OutFileBuf(outfile.c_str(), false);
	} else {
		fout = new OutFileBuf();
	}
	// Initialize Ebwt object and read in header
	if(gVerbose || startVerbose) {
		cerr << "About to initialize fw Ebwt: "; logTime(cerr, true);
	}
	adjIdxBase = adjustEbwtBase(argv0, bt2indexBase, gVerbose);

	// cannot use stack objects on GPUs
	std::unique_ptr<Ebwt> pebwt( new Ebwt(
		adjIdxBase,
		0,        // index is colorspace
		-1,       // fw index
		true,     // index is for the forward direction
		/* overriding: */ offRate,
		0, // amount to add to index offrate or <= 0 to do nothing
		useMm,    // whether to use memory-mapped files
		useShmem, // whether to use shared memory
		mmSweep,  // sweep memory-mapped files
		!noRefNames, // load names?
		true,        // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
		gVerbose, // whether to be talkative
		startVerbose, // talkative during initialization
		false /*passMemExc*/,
		sanityCheck));
	Ebwt& ebwt = *(pebwt.get());

	if(sanityCheck && !os.empty()) {
		// Sanity check number of patterns and pattern lengths in Ebwt
		// against original strings
		assert_eq(os.size(), ebwt.nPat());
		for(size_t i = 0; i < os.size(); i++) {
			assert_eq(os[i].length(), ebwt.plen()[i]);
		}
	}
	// Sanity-check the restored version of the Ebwt
	if(sanityCheck && !os.empty()) {
		ebwt.loadIntoMemory(
			0,
			-1, // fw index
			true, // load SA sample
			true, // load ftab
			true, // load rstarts
			!noRefNames,
			startVerbose);
		ebwt.checkOrigs(os, false, false);
		ebwt.evictFromMemory();
	}
	OutputQueue oq(
		*fout,                           // out file buffer
		reorder && (nthreads > 1), // whether to reorder
		nthreads,                        // # threads
		nthreads > 1 , // whether to be thread-safe
		readsPerBatch,                   // size of output buffer of reads
		0);                      // first read will have this rdid
	{
		Timer _t(cerr, "Time searching: ", timing);
		// Set up pexnalities
		if(bonusMatch > 0 && !localAlign) {
			cerr << "Warning: Match bonus always = 0 in --end-to-end mode; ignoring user setting" << endl;
			bonusMatch = 0;
		}
		// cannot use stack objects on GPUs
		std::unique_ptr<Scoring> psc( new Scoring(
			bonusMatch,     // constant reward for match
			penMmcType,     // how to penalize mismatches
			penMmcMax,      // max mm pelanty
			penMmcMin,      // min mm pelanty
			scoreMin,       // min score as function of read len
			nCeil,          // max # Ns as function of read len
			penNType,       // how to penalize Ns in the read
			penN,           // constant if N pelanty is a constant
			penNCatPair,    // whether to concat mates before N filtering
			penRdGapConst,  // constant coeff for read gap cost
			penRfGapConst,  // constant coeff for ref gap cost
			penRdGapLinear, // linear coeff for read gap cost
			penRfGapLinear, // linear coeff for ref gap cost
			gGapBarrier));  // # rows at top/bot only entered diagonally
		Scoring& sc = *(psc.get());
		EList<size_t> reflens(size_t(32*1024),0);
		for(size_t i = 0; i < ebwt.nPat(); i++) {
			reflens.push_back(ebwt.plen()[i]);
		}
		EList<string> refnames(size_t(32*1024),0);
		readEbwtRefnames(adjIdxBase, refnames);
		SamConfig samc(
			refnames,               // reference sequence names
			reflens,                // reference sequence lengths
			samTruncQname,          // whether to truncate QNAME to 255 chars
			samAppendComment,	// append FASTA/FASTQ comment to SAM record
			samOmitSecSeqQual,      // omit SEQ/QUAL for 2ndary alignments?
			samNoUnal,              // omit unaligned-read records?
			string("bowtie2"),      // program id
			string("bowtie2"),      // program name
			string(BOWTIE2_VERSION), // program version
			argstr,                 // command-line
			rgs_optflag,            // read-group string
			sam_print_as,
			sam_print_xs,
			sam_print_xss,
			sam_print_yn,
			sam_print_xn,
			sam_print_x0,
			sam_print_x1,
			sam_print_xm,
			sam_print_xo,
			sam_print_xg,
			sam_print_nm,
			sam_print_md,
			sam_print_yf,
			sam_print_yi,
			sam_print_ym,
			sam_print_yp,
			sam_print_yt,
			sam_print_ys,
			sam_print_zs,
			sam_print_xr,
			sam_print_xt,
			sam_print_xd,
			sam_print_xu,
			sam_print_yl,
			sam_print_ye,
			sam_print_yu,
			sam_print_xp,
			sam_print_yr,
			sam_print_zb,
			sam_print_zr,
			sam_print_zf,
			sam_print_zm,
			sam_print_zi,
			sam_print_zp,
			sam_print_zu,
			sam_print_zt);
		// Set up hit sink; if sanityCheck && !os.empty() is true,
		// then instruct the sink to "retain" hits in a vector in
		// memory so that we can easily sanity check them later on
		AlnSink *mssink = NULL;
		switch(outType) {
		case OUTPUT_SAM: {
			mssink = new AlnSinkSam(
				oq,           // output queue
				samc,         // settings & routines for SAM output
				refnames,     // reference names
				gQuiet);      // don't print alignment summary at end
			if(!samNoHead) {
				bool printHd = true, printSq = true;
				BTString buf;
				samc.printHeader(buf, rgid, rgs, printHd, !samNoSQ, printSq);
				fout->writeString(buf);
			}
			break;
		}
		default:
			cerr << "Invalid output type: " << outType << endl;
			throw 1;
		}
		if(gVerbose || startVerbose) {
			cerr << "Dispatching to search driver: "; logTime(cerr, true);
		}
		// Do the search for all input reads
		assert(patsrc != NULL);
		assert(mssink != NULL);
                {
			if(gVerbose || startVerbose) {
				cerr << "About to initialize rev Ebwt: "; logTime(cerr, true);
			}

			// We need the mirror index if mismatches are allowed
			// cannot use stack objects on GPUs
			std::unique_ptr<Ebwt> pebwtBw( new Ebwt(
				adjIdxBase + ".rev",
				0,       // index is colorspace
				1,       // TODO: maybe not
				false, // index is for the reverse direction
				/* overriding: */ offRate,
				0, // amount to add to index offrate or <= 0 to do nothing
				useMm,    // whether to use memory-mapped files
				useShmem, // whether to use shared memory
				mmSweep,  // sweep memory-mapped files
				!noRefNames, // load names?
				true,        // load SA sample?
				true,        // load ftab?
				true,        // load rstarts?
				gVerbose,    // whether to be talkative
				startVerbose, // talkative during initialization
				false /*passMemExc*/,
				sanityCheck));
			Ebwt& ebwtBw = *(pebwtBw.get());

			multiseedSearch(
				sc,      // scoring scheme
				pp,      // pattern params
				*patsrc, // pattern source
				*mssink, // hit sink
				ebwt,    // BWT
				&ebwtBw); // BWT'
		}

		// Evict any loaded indexes from memory
		if(ebwt.isInMemory()) {
			ebwt.evictFromMemory();
		}

		if(!gQuiet) {
			size_t repThresh = mhits;
			if(repThresh == 0) {
				repThresh = std::numeric_limits<size_t>::max();
			}
			mssink->finish(
				repThresh,
				gReportDiscordant,
				gReportMixed,
				hadoopOut);
		}
		oq.flush(true);
		assert_eq(oq.numStarted(), oq.numFinished());
		assert_eq(oq.numStarted(), oq.numFlushed());
		delete patsrc;
		delete mssink;
		if(fout != NULL) {
			delete fout;
		}
	}
}

// C++ name mangling is disabled for the bowtie() function to make it
// easier to use Bowtie as a library.
extern "C" {

/**
 * Main bowtie entry function.  Parses argc/argv style command-line
 * options, sets global configuration variables, and calls the driver()
 * function.
 */
int bowtie(int argc, const char **argv) {
	try {
#ifdef WITH_AFFINITY
		//CWILKS: adjust this depending on # of hyperthreads per core
		pinning_observer pinner( 2 /* the number of hyper threads on each core */ );
        	pinner.observe( true );
#endif
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();
		for(int i = 0; i < argc; i++) {
			argstr += argv[i];
			if(i < argc-1) argstr += " ";
		}
		if(startVerbose) { cerr << "Entered main(): "; logTime(cerr, true); }
		parseOptions(argc, argv);
		argv0 = argv[0];
		if(showVersion) {
			cout << argv0 << " version " << BOWTIE2_VERSION << endl;
			if(sizeof(void*) == 4) {
				cout << "32-bit" << endl;
			} else if(sizeof(void*) == 8) {
				cout << "64-bit" << endl;
			} else {
				cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
			}
			cout << "Built on " << BUILD_HOST << endl;
			cout << BUILD_TIME << endl;
			cout << "Compiler: " << COMPILER_VERSION << endl;
			cout << "Options: " << COMPILER_OPTIONS << endl;
			cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
			     << sizeof(int)
			     << ", " << sizeof(long) << ", " << sizeof(long long)
			     << ", " << sizeof(void *) << ", " << sizeof(size_t)
			     << ", " << sizeof(off_t) << "}" << endl;
			return 0;
		}
		{
			Timer _t(cerr, "Overall time: ", timing);
			if(startVerbose) {
				cerr << "Parsing index and read arguments: "; logTime(cerr, true);
			}

			// Get index basename (but only if it wasn't specified via --index)
			if(bt2index.empty()) {
				cerr << "No index, query, or output file specified!" << endl;
				printUsage(cerr);
				return 1;
			}

			// Get query filename
			bool got_reads = !queries.empty() || !mates1.empty() || !mates12.empty();
#ifdef USE_SRA
			got_reads = got_reads || !sra_accs.empty();
#endif
			if(optind >= argc) {
				if(!got_reads) {
					printUsage(cerr);
					cerr << "***" << endl
#ifdef USE_SRA
					     << "Error: Must specify at least one read input with -U/-1/-2/--sra-acc" << endl;
#else
					<< "Error: Must specify at least one read input with -U/-1/-2" << endl;
#endif
					return 1;
				}
			} else if(!got_reads) {
				// Tokenize the list of query files
				tokenize(argv[optind++], ",", queries);
				if(queries.empty()) {
					cerr << "Tokenized query file list was empty!" << endl;
					printUsage(cerr);
					return 1;
				}
			}

			// Get output filename
			if(optind < argc && outfile.empty()) {
				outfile = argv[optind++];
				cerr << "Warning: Output file '" << outfile.c_str()
				     << "' was specified without -S.  This will not work in "
				     << "future Bowtie 2 versions.  Please use -S instead."
				     << endl;
			}

			// Extra parametesr?
			if(optind < argc) {
				cerr << "Extra parameter(s) specified: ";
				for(int i = optind; i < argc; i++) {
					cerr << "\"" << argv[i] << "\"";
					if(i < argc-1) cerr << ", ";
				}
				cerr << endl;
				if(mates1.size() > 0) {
					cerr << "Note that if <mates> files are specified using -1/-2, a <singles> file cannot" << endl
					     << "also be specified.  Please run bowtie separately for mates and singles." << endl;
				}
				throw 1;
			}

			// Optionally summarize
			if(gVerbose) {
				cout << "Input " + gEbwt_ext +" file: \"" << bt2index.c_str() << "\"" << endl;
				cout << "Query inputs (DNA, " << file_format_names[format].c_str() << "):" << endl;
				for(size_t i = 0; i < queries.size(); i++) {
					cout << "  " << queries[i].c_str() << endl;
				}
				cout << "Quality inputs:" << endl;
				for(size_t i = 0; i < qualities.size(); i++) {
					cout << "  " << qualities[i].c_str() << endl;
				}
				cout << "Output file: \"" << outfile.c_str() << "\"" << endl;
				cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
				cout << "Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
#ifdef NDEBUG
				cout << "Assertions: disabled" << endl;
#else
				cout << "Assertions: enabled" << endl;
#endif
			}
			if(ipause) {
				cout << "Press key to continue..." << endl;
				getchar();
			}
			driver<SString<char> >("DNA", bt2index, outfile);
		}
#ifdef WITH_AFFINITY
		// Always disable observation before observers destruction
    		//tracker.observe( false );
    		pinner.observe( false );
#endif
		return 0;
	} catch(std::exception& e) {
		cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Error: Encountered internal Bowtie 2 exception (#" << e << ")" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
} // bowtie()
} // extern "C"
