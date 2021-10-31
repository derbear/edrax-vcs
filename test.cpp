#include "vcs.h"

#include "benchmark.h"

#include <iostream>
#include <random>
#include <thread>
#include <vector>
#include <sstream>
#include <map>

#include "test_point.hpp"
#include "bn.h"

#include <gmp.h>
#include <gmpxx.h>

using namespace std;
using namespace bn;

int main(int argc, char** argv){
	// init
	int L = atoi(argv[1]);

	bn::CurveParam cp = bn::CurveFp254BNb;
	Param::init(cp);
	const Point& pt = selectPoint(cp);
	const Ec2 g2(
		Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
		Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
	);
	const Ec1 g1(pt.g1.a, pt.g1.b);
	
	mpz_class p;
	p.set_str("16798108731015832284940804142231733909759579603404752749028378864165570215949",10);

	// from cppreference.com
	random_device rd("/dev/urandom");  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	// seed_seq rseed{0, 0, 0, 0, 0, 0, 0, 0};
	// mt19937 gen(rseed); //Standard mersenne_twister_engine seeded with rd()

	uniform_int_distribution<> distrib(0, numeric_limits<int>::max());

	// cout << distrib(gen) << endl;
	// cout << distrib(gen) << endl;
	// cout << distrib(gen) << endl;

	// srand(time(NULL));
	gmp_randstate_t r_state;
	unsigned long int seed = distrib(gen);
	gmp_randinit_default (r_state);
	gmp_randseed_ui(r_state, seed);
	
	vcs a(L,p,g1,g2);
	
	auto t2 = chrono::steady_clock::now();
	
	vector<vector<Ec1> > prk;
	vector<Ec2> vrk;

	a.keygen(prk, vrk);
	a.load_key(prk,vrk);
	
	auto t3 = chrono::steady_clock::now();
	// auto t4 = t3 - t2;
	auto t4 = chrono::duration<double, micro>(t3 - t2);
	// cout << "keygen time: " << chrono::duration<double, milli>(t4).count()/1000 << "s" << endl;
	auto t1 = t4;
	auto tmin = t4;
	auto tmax = t4;
	
	// start
	int errs = 0;
	// cout << "B," << a.N << endl;
	auto iters = 1000;
	auto warm = 100;
	auto tot_iters = warm + iters;
	// cout << "iters," << iters << endl;

	// commit
	Ec1 digest = g1*0;
	vector<mpz_class> vals(a.N);

	vector<chrono::duration<double, micro>> commit_m(tot_iters);
	t1 = chrono::duration<double, micro>::zero();
	// for (int k = 0; k < tot_iters; k++) {
	for (int k = 0; k < 1; k++) { // cut to save time (not measuring commit)
	  for (int i = 0; i < a.N; i++) {
	    vals[i] = distrib(gen);
	  }

	  t2 = chrono::steady_clock::now();
	  {
	    for (int i = 0; i < a.N; i++) {
	      auto upk_i = a.calc_update_key(i, prk);
	      benchmark::DoNotOptimize(digest = a.update_digest(digest, i, vals[i], upk_i));
	      benchmark::ClobberMemory();
	    }
	  }
	  t3 = chrono::steady_clock::now();
	  t4 = chrono::duration<double, micro>(t3 - t2);
	  commit_m[k] = t4;
	}

	t1 = commit_m[warm];
	tmin = t1;
	tmax = t1;
	for (int k = warm+1; k < tot_iters; k++) {
	  t1 += commit_m[k];
	  if (commit_m[k] > tmax) {
	    tmax = commit_m[k];
	  }
	  if (commit_m[k] < tmin) {
	    tmin = commit_m[k];
	  }
	}
	// cout << "commit," << int(t1.count()/iters) << "," << int(tmin.count()) << "," << int(tmax.count()) << endl;

	// open
	vector<int> open_indexes(tot_iters);
	vector<vector<Ec1> > proofs(tot_iters);
	for (int j = 0; j < tot_iters; j++) {
	  open_indexes[j] = distrib(gen) % a.N;
	}

	{
	  vector<chrono::duration<double, micro>> open_m(tot_iters);
	  for (int j = 0; j < tot_iters; j++) {
	    auto i = open_indexes[j];
	    t2 = chrono::steady_clock::now();
	    benchmark::DoNotOptimize(proofs[j] = a.prove(i, vals, prk));
	    benchmark::ClobberMemory();
	    t3 = chrono::steady_clock::now();
	    t4 = chrono::duration<double, micro>(t3 - t2);
	    open_m[j] = t4;
	  }
	  t1 = open_m[warm];
	  tmin = t1;
	  tmax = t1;
	  for (int k = warm+1; k < tot_iters; k++) {
	    t1 += open_m[k];
	    if (open_m[k] > tmax) {
	      tmax = open_m[k];
	    }
	    if (open_m[k] < tmin) {
	      tmin = open_m[k];
	    }
	  }
	  // cout << "open," << int(t1.count()/iters) << "," << int(tmin.count()) << "," << int(tmax.count()) << endl;
	}

	// verify
	{
	  vector<chrono::duration<double, micro>> verify_m(tot_iters);
	  for (int j = 0; j < tot_iters; j++) {
	    auto i = open_indexes[j];

	    t2 = chrono::steady_clock::now();
	    bool ok;
	    benchmark::DoNotOptimize(ok = a.verify(digest, i, vals[i], proofs[j], vrk));
	    if (!ok) {
	      errs += 1;
	    }
	    benchmark::ClobberMemory();
	    t3 = chrono::steady_clock::now();
	    t4 = chrono::duration<double, micro>(t3 - t2);
	    verify_m[j] = t4;
	  }
	  cout << "verify,";
	  for (int j = warm; j < tot_iters; j++) {
	    cout << int(verify_m[j].count()) << ",";
	  }
	  cout << endl;
	  t1 = verify_m[warm];
	  tmin = t1;
	  tmax = t1;
	  for (int k = warm+1; k < tot_iters; k++) {
	    t1 += verify_m[k];
	    if (verify_m[k] > tmax) {
	      tmax = verify_m[k];
	    }
	    if (verify_m[k] < tmin) {
	      tmin = verify_m[k];
	    }
	  }
	  // cout << "verify," << int(t1.count()/iters) << "," << int(tmin.count()) << "," << int(tmax.count()) << endl;
	}
	
	// update commit
	vector<long long int> update_indexes(tot_iters);
	vector<int> update_vals(tot_iters);
	for (int j = 0; j < tot_iters; j++) {
	  update_indexes[j] = distrib(gen) % a.N;
	  update_vals[j] = distrib(gen);
	}
	auto upk_is = a.calc_update_key_batch(update_indexes, prk);

	{
	  vector<chrono::duration<double, micro>> cupdate_m(tot_iters);
	  for (int j = 0; j < tot_iters; j++) {
	    auto i = update_indexes[j];
	    auto upk_i = upk_is[j];

	    t2 = chrono::steady_clock::now();
	    benchmark::DoNotOptimize(digest = a.update_digest(digest, i, update_vals[j], upk_i));
	    benchmark::ClobberMemory();
	    t3 = chrono::steady_clock::now();
	    t4 = chrono::duration<double, micro>(t3 - t2);;
	    cupdate_m[j] = t4;
	  }
	  cout << "commit_update,";
	  for (int j = warm; j < tot_iters; j++) {
	    cout << int(cupdate_m[j].count()) << ",";
	  }
	  cout << endl;
	  t1 = cupdate_m[warm];
	  tmin = t1;
	  tmax = t1;
	  for (int k = warm+1; k < tot_iters; k++) {
	    t1 += cupdate_m[k];
	    if (cupdate_m[k] > tmax) {
	      tmax = cupdate_m[k];
	    }
	    if (cupdate_m[k] < tmin) {
	      tmin = cupdate_m[k];
	    }
	  }
	  // cout << "commit_update," << int(t1.count()/iters) << "," << int(tmin.count()) << "," << int(tmax.count()) << endl;
	}
	
	// update proof
	{
	  vector<chrono::duration<double, micro>> pupdate_m(tot_iters);
	  for (int j = 0; j < tot_iters; j++) {
	    auto i = update_indexes[j];
	    auto upk_i = upk_is[j];

	    auto l = distrib(gen) % a.N;
	    auto proofl = a.prove(l, vals, prk);

	    t2 = chrono::steady_clock::now();
	    benchmark::DoNotOptimize(proofl = a.update_proof(proofl, i, l, update_vals[j], upk_i));
	    benchmark::ClobberMemory();
	    t3 = chrono::steady_clock::now();
	    t4 = chrono::duration<double, micro>(t3 - t2);;
	    pupdate_m[j] = t4;
	  }
	  cout << "proof_update,";
	  for (int j = warm; j < tot_iters; j++) {
	    cout << int(pupdate_m[j].count()) << ",";
	  }
	  cout << endl;
	  t1 = pupdate_m[warm];
	  tmin = t1;
	  tmax = t1;
	  for (int k = warm+1; k < tot_iters; k++) {
	    t1 += pupdate_m[k];
	    if (pupdate_m[k] > tmax) {
	      tmax = pupdate_m[k];
	    }
	    if (pupdate_m[k] < tmin) {
	      tmin = pupdate_m[k];
	    }
	  }
	  // cout << "proof_update," << int(t1.count()/iters) << "," << int(tmin.count()) << "," << int(tmax.count()) << endl;
	}

	// 100 micros sanity check
	{
	  vector<chrono::duration<double, micro>> sanity_m(tot_iters);
	  for (int j = 0; j < tot_iters; j++) {
	    auto i = update_indexes[j];
	    auto upk_i = upk_is[j];

	    auto l = distrib(gen) % a.N;
	    auto proofl = a.prove(l, vals, prk);

	    t2 = chrono::steady_clock::now();
	    this_thread::sleep_for(chrono::duration<double, micro>(100));
	    benchmark::ClobberMemory();
	    t3 = chrono::steady_clock::now();
	    t4 = chrono::duration<double, micro>(t3 - t2);;
	    sanity_m[j] = t4;
	  }
	  cout << "100micros,";
	  for (int j = warm; j < tot_iters; j++) {
	    cout << int(sanity_m[j].count()) << ",";
	  }
	  cout << endl;
	  t1 = sanity_m[warm];
	  tmin = t1;
	  tmax = t1;
	  for (int k = warm+1; k < tot_iters; k++) {
	    t1 += sanity_m[k];
	    if (sanity_m[k] > tmax) {
	      tmax = sanity_m[k];
	    }
	    if (sanity_m[k] < tmin) {
	      tmin = sanity_m[k];
	    }
	  }
	  // cout << "100micros," << int(t1.count()/iters) << "," << int(tmin.count()) << "," << int(tmax.count()) << endl;
	}
	
	if (errs != 0) {
	  cout << "errs," << errs << endl;
	  return 1;
	}
	return 0;
}
